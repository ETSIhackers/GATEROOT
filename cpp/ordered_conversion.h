#pragma once
#include "event_batch.h"
#include <filesystem>
#include <future>
#include <queue>
#include <thread>
#include <limits>
#include <chrono>
#include <tuple>
#include <TFile.h>
#include <TTree.h>

struct ConversionOptions {
  std::string execution = "auto", clock = "auto", delayeds = "auto";
  size_t workers = 0, batch_size = 500000;
  double duration_ms = 10000, start_ms = 0, stop_ms = std::numeric_limits<double>::infinity();
  std::filesystem::path temporary;
};

struct FileSpan {
  double low = std::numeric_limits<double>::infinity();
  double high = -std::numeric_limits<double>::infinity();
  bool delayed = false;
};

inline FileSpan file_span(const std::string& filename, bool include_delayeds) {
  std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("Cannot read ROOT file " + filename);
  FileSpan span;
  span.delayed = file->Get("delay") != nullptr;
  for (const auto* name : {"Coincidences", "delay"}) {
    if (std::string(name) == "delay" && (!include_delayeds || !span.delayed)) continue;
    auto* tree = dynamic_cast<TTree*>(file->Get(name));
    if (!tree) throw std::runtime_error("Missing coincidence tree in " + filename);
    double clock;
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("time1",1);
    if (tree->SetBranchAddress("time1",&clock) < 0) throw std::runtime_error("Missing time1 branch");
    for (Long64_t i=0;i<tree->GetEntries();++i) {
      tree->GetEntry(i);
      if (!std::isfinite(clock) || clock < 0) throw std::runtime_error("Invalid ROOT timestamp");
      span.low = std::min(span.low,clock*1000.);
      span.high = std::max(span.high,clock*1000.);
    }
    tree->ResetBranchAddresses();
  }
  return span;
}

struct SpoolRecord {
  double time_ms;
  uint64_t entry;
  uint32_t file, kind, a, b, tof;
};

inline bool record_less(const SpoolRecord& a,const SpoolRecord& b) {
  return std::tie(a.time_ms,a.file,a.kind,a.entry) < std::tie(b.time_ms,b.file,b.kind,b.entry);
}

class SpoolMerge {
  struct Cursor { SpoolRecord row; size_t index; };
  struct Greater { bool operator()(const Cursor& a,const Cursor& b) const {
    return record_less(b.row,a.row);
  }};
  std::vector<std::ifstream> readers;
  std::priority_queue<Cursor,std::vector<Cursor>,Greater> queue;
public:
  explicit SpoolMerge(const std::vector<std::filesystem::path>& paths) {
    readers.reserve(paths.size());
    for (const auto& path : paths) {
      readers.emplace_back(path,std::ios::binary);
      if (!readers.back()) throw std::runtime_error("Cannot open event spool");
      SpoolRecord row;
      if (readers.back().read(reinterpret_cast<char*>(&row),sizeof row)) queue.push({row,readers.size()-1});
    }
  }
  bool next(SpoolRecord& out) {
    if (queue.empty()) return false;
    auto head = queue.top(); queue.pop(); out = head.row;
    auto& stream = readers[head.index];
    if (stream.read(reinterpret_cast<char*>(&head.row),sizeof head.row)) queue.push(head);
    else if (!stream.eof()) throw std::runtime_error("Cannot read event spool");
    return true;
  }
};

inline std::filesystem::path merge_file_runs(std::vector<std::filesystem::path> runs,
    const std::filesystem::path& directory, size_t file) {
  size_t pass = 0;
  while (runs.size() > 1) {
    std::vector<std::filesystem::path> next;
    for (size_t begin=0;begin<runs.size();begin+=32) {
      auto end = std::min(runs.size(),begin+32);
      std::vector<std::filesystem::path> group(runs.begin()+begin,runs.begin()+end);
      auto path = directory/("merge_"+std::to_string(file)+"_"+std::to_string(pass)+"_"+std::to_string(begin));
      {
        std::ofstream output(path,std::ios::binary);
        if (!output) throw std::runtime_error("Cannot create merged event spool");
        SpoolMerge merge(group); SpoolRecord row;
        while (merge.next(row)) output.write(reinterpret_cast<const char*>(&row),sizeof row);
        if (!output) throw std::runtime_error("Event spool disk is full");
      }
      for (const auto& used : group) std::filesystem::remove(used);
      next.push_back(path);
    }
    runs = std::move(next); ++pass;
  }
  return runs.empty() ? std::filesystem::path{} : runs[0];
}

struct FileResult {
  std::filesystem::path spool;
  uint64_t prompts = 0, delayeds = 0;
};

inline FileResult spool_root_file(const std::string& prefix, int file_index, double offset,
    bool delayeds, double delayed_offset_ps, const ScannerGeometry& geometry,
    const petsird::ScannerInformation& scanner, const std::array<bool,num_categories>& selected,
    const ConversionOptions& options, bool gpu, const std::filesystem::path& directory) {
  FileResult result;
  std::vector<std::filesystem::path> runs;
  const auto& axis = scanner.tof_bin_edges[0][0].edges;
  std::vector<float> edges(axis.begin(),axis.end());
  const auto& energy = scanner.event_energy_bin_edges[0].edges;
  for (int kind=0;kind<(delayeds?2:1);++kind) {
    GateCoincidences tree(kind?"delay":"Coincidences",prefix,1,file_index);
    const auto total = tree.chain.GetEntries();
    for (Long64_t begin=0;begin<total;begin+=options.batch_size) {
      std::vector<BatchInput> raw;
      const auto end = std::min(total,begin+static_cast<Long64_t>(options.batch_size));
      raw.reserve(end-begin);
      for (auto entry=begin;entry<end;++entry) {
        if (tree.chain.GetEntry(entry) <= 0) throw std::runtime_error("ROOT entry could not be read");
        if (!std::isfinite(tree.time1) || !std::isfinite(tree.time2) || tree.time1<0)
          throw std::runtime_error("Invalid ROOT coincidence timing");
        uint32_t ids[2];
        petsird::ExpandedDetectionBin expanded;
        expanded.energy_index = 0;
        expanded.module_index = calculate_module_index(tree.gantryID1,tree.rsectorID1,geometry);
        expanded.element_index = calculate_element_index(tree.moduleID1,tree.submoduleID1,tree.crystalID1,tree.layerID1,geometry);
        ids[0] = petsird_helpers::make_detection_bin(scanner,petsird::TypeOfModule{0},expanded);
        expanded.module_index = calculate_module_index(tree.gantryID2,tree.rsectorID2,geometry);
        expanded.element_index = calculate_element_index(tree.moduleID2,tree.submoduleID2,tree.crystalID2,tree.layerID2,geometry);
        ids[1] = petsird_helpers::make_detection_bin(scanner,petsird::TypeOfModule{0},expanded);
        raw.push_back({tree.time1*1000.+offset,(tree.time1-tree.time2)*1e12+(kind?delayed_offset_ps:0.),
                       tree.energy1*1000.,tree.energy2*1000.,static_cast<uint64_t>(entry),ids[0],ids[1],
                       kind || selected[tree.category()]});
      }
      auto transformed = gpu ? gpu_transform(raw,edges,energy[0],energy[1],options.start_ms,options.stop_ms)
                             : cpu_transform(raw,edges,energy[0],energy[1],options.start_ms,options.stop_ms);
      if (transformed.empty()) continue;
      auto path = directory/("run_"+std::to_string(file_index)+"_"+std::to_string(kind)+"_"+std::to_string(begin));
      std::ofstream stream(path,std::ios::binary);
      for (const auto& e : transformed) {
        SpoolRecord record{e.time_ms,e.entry,static_cast<uint32_t>(file_index-1),
                           static_cast<uint32_t>(kind),e.a,e.b,e.tof};
        stream.write(reinterpret_cast<const char*>(&record),sizeof record);
      }
      if (!stream) throw std::runtime_error("Cannot write event spool; check free disk space");
      runs.push_back(path);
      (kind?result.delayeds:result.prompts) += transformed.size();
    }
  }
  result.spool = merge_file_runs(std::move(runs),directory,file_index);
  return result;
}

inline int ordered_conversion(const std::string& prefix,int count,const std::string& output,
    const ScannerGeometry& geometry,petsird::Header header,const std::array<bool,num_categories>& selected,
    double delayed_offset_ps,ConversionOptions options) {
  const auto started=std::chrono::steady_clock::now();
  if (count<0 || options.batch_size==0 || !std::isfinite(options.duration_ms) || options.duration_ms<=0
      || !std::isfinite(options.start_ms) || options.start_ms<0 || !(options.stop_ms>options.start_ms))
    throw std::runtime_error("Invalid file count, batch size, duration or acquisition time window");
  const bool gpu = options.execution != "cpu" && options.execution != "none" && gpu_available();
  if (options.execution == "gpu" && !gpu) throw std::runtime_error("Requested GPU mode is unavailable; build with CUDA and expose a GPU");
  const size_t workers = options.execution=="none" ? 1 : options.workers ? options.workers : std::min(8u,std::max(1u,std::thread::hardware_concurrency()));
  ROOT::EnableThreadSafety();
  std::cout << "Execution: " << (gpu?"gpu":options.execution=="none"?"none":"cpu")
            << ", " << workers << " independent ROOT readers" << std::endl;
  std::vector<FileSpan> spans;
  std::vector<std::future<FileSpan>> probes;
  for (int i=1;i<=count;++i) {
    probes.push_back(std::async(options.execution=="none"?std::launch::deferred:std::launch::async,[&,i] {
      return file_span(prefix+std::to_string(i)+".root",options.delayeds!="off");
    }));
    if (probes.size()>=workers || i==count) {
      for (auto& probe:probes) spans.push_back(probe.get());
      probes.clear();
    }
  }
  std::vector<double> offsets;
  double shift = 0, previous = -std::numeric_limits<double>::infinity();
  bool any_delays = false;
  for (size_t i=0;i<spans.size();++i) {
    const auto& span = spans[i];
    if (options.delayeds=="required" && !span.delayed) throw std::runtime_error("Required delay tree is missing");
    any_delays |= span.delayed && options.delayeds!="off";
    if (options.clock=="per-file") {
      shift=i*options.duration_ms;
      if (span.high>=options.duration_ms) throw std::runtime_error("Per-file timestamps exceed file duration");
    } else if (options.clock=="auto" && span.low+shift<previous) {
      shift=std::ceil(previous/options.duration_ms)*options.duration_ms;
    }
    offsets.push_back(shift);
    if (std::isfinite(span.high)) previous=std::max(previous,span.high+shift);
  }
  header.scanner.delayed_event_policy=any_delays?petsird::CoincidencePolicy::kOther:petsird::CoincidencePolicy::kNone;
  auto destination=std::filesystem::path(output), partial=std::filesystem::path(output+".partial");
  auto manifest=std::filesystem::path(output+".json");
  if (std::filesystem::exists(destination) || std::filesystem::exists(partial) || std::filesystem::exists(manifest))
    throw std::runtime_error("Output exists; select a new PETSIRD filename");
  if (!destination.parent_path().empty()) std::filesystem::create_directories(destination.parent_path());
  auto directory=options.temporary.empty()?std::filesystem::temp_directory_path():options.temporary;
  directory/=std::string("gateroot-")+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
  std::filesystem::create_directories(directory);
  std::vector<std::future<FileResult>> pending;
  std::vector<std::filesystem::path> spools;
  nlohmann::json report={{"converter_version","1.0.0"},{"sdk","0.13.3"},
      {"execution_requested",options.execution},{"execution_used",gpu?"gpu":options.execution=="none"?"none":"cpu"},
      {"workers",workers},{"clock_offsets_ms",offsets},{"start_time_s",options.start_ms/1000.},
      {"files",nlohmann::json::array()}};
  report["end_time_s"]=std::isfinite(options.stop_ms)?nlohmann::json(options.stop_ms/1000.):nlohmann::json(nullptr);
  try {
    for (int i=1;i<=count;++i) {
      pending.push_back(std::async(options.execution=="none"?std::launch::deferred:std::launch::async,[&,i] {
        if (spans[i-1].high+offsets[i-1]<options.start_ms || spans[i-1].low+offsets[i-1]>=options.stop_ms)
          return FileResult{};
        return spool_root_file(prefix,i,offsets[i-1],spans[i-1].delayed && options.delayeds!="off",
                               delayed_offset_ps,geometry,header.scanner,selected,options,gpu,directory);
      }));
      if (pending.size()>=workers || i==count) {
        const size_t begin=static_cast<size_t>(i)-pending.size();
        for (size_t j=0;j<pending.size();++j) {
          auto result=pending[j].get();
          if (!result.spool.empty()) spools.push_back(result.spool);
          report["files"].push_back({{"file",prefix+std::to_string(begin+j+1)+".root"},
              {"prompts",result.prompts},{"delayeds",result.delayeds},{"offset_ms",offsets[begin+j]}});
        }
        pending.clear();
      }
    }
    report["read_transform_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    petsird::binary::PETSIRDWriter writer(partial.string()); writer.WriteHeader(header);
    SpoolMerge merge(spools); SpoolRecord row;
    bool has=merge.next(row);
    while (has) {
      const double tag=std::floor(row.time_ms/geometry.LM_time_block_duration_ms)*geometry.LM_time_block_duration_ms;
      if (tag+geometry.LM_time_block_duration_ms>std::numeric_limits<uint32_t>::max()) throw std::runtime_error("PETSIRD timestamp overflow");
      auto block=create_new_time_block(header.scanner,any_delays);
      block.time_interval.start=static_cast<uint32_t>(tag);
      block.time_interval.stop=static_cast<uint32_t>(tag+geometry.LM_time_block_duration_ms);
      do {
        petsird::CoincidenceEvent event; event.detection_bins={row.a,row.b}; event.tof_idx=row.tof;
        if (row.kind) block.delayed_events[0][0].push_back(event);
        else block.prompt_events[0][0].push_back(event);
        has=merge.next(row);
      } while (has && row.time_ms<tag+geometry.LM_time_block_duration_ms);
      writer.WriteTimeBlocks(block);
    }
    writer.EndTimeBlocks(); writer.Close();
    std::filesystem::rename(partial,destination);
    report["total_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::ofstream report_file(manifest); report_file << report.dump(2) << '\n';
  } catch (...) {
    for (auto& future:pending) if (future.valid()) future.wait();
    std::filesystem::remove(partial);
    std::filesystem::remove_all(directory);
    throw;
  }
  std::filesystem::remove_all(directory);
  std::cout << report.dump(2) << std::endl;
  return 0;
}
