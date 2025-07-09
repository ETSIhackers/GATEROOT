//  *********************************************************************
//  * DISCLAIMER                                                        *
//  *                                                                   *
//  * Neither the authors of this software system, nor their employing  *
//  * institutes, nor the agencies providing financial support for this *
//  * work  make  any representation or  warranty, express or implied,  *
//  * regarding  this  software system or assume any liability for its  *
//  * use.                                                              *
//  *                                                                   *
//  * This  code  implementation is the  intellectual property  of the  *
//  * OpenGATE collaboration.                                           *
//  * By copying,  distributing  or modifying the Program (or any work  *
//  * based  on  the Program)  you indicate  your  acceptance of  this  *
//  * statement, and all its terms.                                     *
//  *********************************************************************
//
//######################################################################################
//# Authors   : Nicolas A Karakatsanis,                                                #
//#                                                                                    #
//# Date  03-NOV-2024                                                                  #
//#                                                                                    #
//# Objective : To read the coincidences TTree from the .root file, and generates the  #
//#             corresponding Michelogram and Projection files.                        #
//#                                                                                    #
//# Input     : Monte Carlo data from GATE                                             #
//#                                                                                    #
//# Output    : 1 Michelogram  files according to various binning definitions          #
//#           : 2 Projection  files according to various binning definitions           #
//#                                                                                    #
//######################################################################################
//#                                                                                    #
//# This file is last modified on Nov 03, 2024 by: N. Karakatsanis                     #
//#                                                                                    #
//# The data are input from a root file produced by Gate simulating extended FOV of    #
//# mCT scanner. This scanner will have 5xFOV thus the root file contains information  #
//# on every gantry. In this case there are 5 gantries with gantryID (0 -> 4).         #
//# The ring numbers are defined based on the gantryID.                                #
//#                                                                                    #
//# The virtual rings between the blocks are taken into consideration here.            #
//#                                                                                    #
//# The central FOV is taken into consideration                                        #
//# The maximum and minimum rings should be specified if the user wishes to change     #
//# the number or the order of gantries.                                               #
//#                                                                                    #
//#                                                                                    #
//#                                                                                    #
//#                                                                                    #
//#  NEW WAY TO RUN:                                                                   #
//#                                                                                    #
//# HOW TO COMPILE:                                                                    #
//# 1) Compile using this command line in the terminal:                                #
//#       g++ (name of file) `root-config --cflags --libs`                             #
//#                                                                                    #
//# HOW TO RUN:                                                                        #
//# 1) After compiling the code type the following command line:                       #
//#       ./a.out 'directory name of root files' 'output file name' 'yes/no'           #
//#                                                                                    #
//# NOTE: 1) To drop the odd ring put 'yes' as an argument when running, otherwise the #
//# rings will not be droped.                                                          #
//#       2) You will be prompted to enter the minimum ring number                     #
//######################################################################################



#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "TROOT.h"
#include "TChain.h"
#include <cmath>

#include <nlohmann/json.hpp>

// PETSIRD Includes
#include "protocols.h"
#include "types.h"
#include "binary/protocols.h"
#include "petsird_helpers/create.h"
#include "petsird_helpers/geometry.h"

using namespace std ;

//#include "petsird_helpers.h"

struct ScannerGeometry
{
  int n_rings;
  int n_det;
  int s_width;
  int n_rsec;
  int n_rsec_xy;
  int n_rsec_z;
  int n_module;
  int n_mod_xy;
  int n_mod_z;
  int n_submod;
  int n_smod_xy;
  int n_smod_z;
  int n_crystal;
  int n_cry_xy;
  int n_cry_z;
  int n_cry_layers;
  int cry_ax_gap;
  int cry_tx_gap;
  int smod_ax_gap;
  int smod_tx_gap;
  int mod_ax_gap;
  int mod_tx_gap;
  int rsec_ax_gap;
  int rsec_tx_gap;
  int max_d_ring;
  int number_of_tof_bins;
  int number_of_energy_bins;
  float radius;
  int tx_virtual_crystal_num;
  int ax_virtual_crystal_num;
  int tx_phys_crystal_num;
  int ax_phys_crystal_num;
  float detector_x_dim, detector_y_dim, detector_z_dim;
  float energy_LLD, energy_ULD;
  float EnergyResolutionAt511;
  float TOF_resolution;
  float LM_TimeBlockDuration;
  float ArcLength;
  float TxFOV;
  float TxFOV_TOF;
  float module_axial_pitch;
};

void WriteScannerGeometry(const ScannerGeometry& scanner_geometry, const std::string& filename)
{
  nlohmann::json j;
  j["n_rings"] = scanner_geometry.n_rings;
  j["n_det"] = scanner_geometry.n_det;
  j["s_width"] = scanner_geometry.s_width;
  j["n_rsec"] = scanner_geometry.n_rsec;
  j["n_rsec_xy"] = scanner_geometry.n_rsec_xy;
  j["n_rsec_z"] = scanner_geometry.n_rsec_z;
  j["n_module"] = scanner_geometry.n_module;
  j["n_mod_xy"] = scanner_geometry.n_mod_xy;
  j["n_mod_z"] = scanner_geometry.n_mod_z;
  j["n_submod"] = scanner_geometry.n_submod;
  j["n_smod_xy"] = scanner_geometry.n_smod_xy;
  j["n_smod_z"] = scanner_geometry.n_smod_z;
  j["n_crystal"] = scanner_geometry.n_crystal;
  j["n_cry_xy"] = scanner_geometry.n_cry_xy;
  j["n_cry_z"] = scanner_geometry.n_cry_z;
  j["n_cry_layers"] = scanner_geometry.n_cry_layers;
  j["cry_ax_gap"] = scanner_geometry.cry_ax_gap;
  j["cry_tx_gap"] = scanner_geometry.cry_tx_gap;
  j["smod_ax_gap"] = scanner_geometry.smod_ax_gap;
  j["smod_tx_gap"] = scanner_geometry.smod_tx_gap;
  j["mod_ax_gap"] = scanner_geometry.mod_ax_gap;
  j["mod_tx_gap"] = scanner_geometry.mod_tx_gap;
  j["rsec_ax_gap"] = scanner_geometry.rsec_ax_gap;
  j["rsec_tx_gap"] = scanner_geometry.rsec_tx_gap;
  j["max_d_ring"] = scanner_geometry.max_d_ring;
  j["number_of_tof_bins"] = scanner_geometry.number_of_tof_bins;
  j["number_of_energy_bins"] = scanner_geometry.number_of_energy_bins;
  j["radius"] = scanner_geometry.radius;
  j["tx_virtual_crystal_num"] = scanner_geometry.tx_virtual_crystal_num;
  j["ax_virtual_crystal_num"] = scanner_geometry.ax_virtual_crystal_num;
  j["tx_phys_crystal_num"] = scanner_geometry.tx_phys_crystal_num;
  j["ax_phys_crystal_num"] = scanner_geometry.ax_phys_crystal_num;
  j["detector_x_dim"] = scanner_geometry.detector_x_dim;
  j["detector_y_dim"] = scanner_geometry.detector_y_dim;
  j["detector_z_dim"] = scanner_geometry.detector_z_dim;
  j["energy_LLD"] = scanner_geometry.energy_LLD;
  j["energy_ULD"] = scanner_geometry.energy_ULD;
  j["EnergyResolutionAt511"] = scanner_geometry.EnergyResolutionAt511;
  j["TOF_resolution"] = scanner_geometry.TOF_resolution;
  j["LM_TimeBlockDuration"] = scanner_geometry.LM_TimeBlockDuration;

  j["ArcLength"] = scanner_geometry.s_width * scanner_geometry.detector_y_dim / 2.0f;
  j["TxFOV"] = 2 * scanner_geometry.radius * sin (scanner_geometry.ArcLength / (2 * scanner_geometry.radius) );
  j["TxFOV_TOF"] = scanner_geometry.TxFOV + 0.3*scanner_geometry.TOF_resolution;
  j["module_axial_pitch"] = scanner_geometry.n_cry_z * scanner_geometry.detector_z_dim + (scanner_geometry.n_cry_z - 1) * scanner_geometry.cry_ax_gap;


  std::ofstream o(filename);
  o << std::setw(4) << j << std::endl;
}

// Function for reading json scanner geometry
ScannerGeometry ReadScannerGeometry(const std::string& filename)
{
  std::ifstream i(filename);
  nlohmann::json j;
  i >> j;

  ScannerGeometry scanner_geometry;
  scanner_geometry.n_rings = j["n_rings"];
  scanner_geometry.n_det = j["n_det"];
  scanner_geometry.s_width = j["s_width"];
  scanner_geometry.n_rsec = j["n_rsec"];
  scanner_geometry.n_rsec_xy = j["n_rsec_xy"];
  scanner_geometry.n_rsec_z = j["n_rsec_z"];
  scanner_geometry.n_module = j["n_module"];
  scanner_geometry.n_mod_xy = j["n_mod_xy"];
  scanner_geometry.n_mod_z = j["n_mod_z"];
  scanner_geometry.n_submod = j["n_submod"];
  scanner_geometry.n_smod_xy = j["n_smod_xy"];
  scanner_geometry.n_smod_z = j["n_smod_z"];
  scanner_geometry.n_crystal = j["n_crystal"];
  scanner_geometry.n_cry_xy = j["n_cry_xy"];
  scanner_geometry.n_cry_z = j["n_cry_z"];
  scanner_geometry.n_cry_layers = j["n_cry_layers"];
  scanner_geometry.cry_ax_gap = j["cry_ax_gap"];
  scanner_geometry.cry_tx_gap = j["cry_tx_gap"];
  scanner_geometry.smod_ax_gap = j["smod_ax_gap"];
  scanner_geometry.smod_tx_gap = j["smod_tx_gap"];
  scanner_geometry.mod_ax_gap = j["mod_ax_gap"];
  scanner_geometry.mod_tx_gap = j["mod_tx_gap"];
  scanner_geometry.rsec_ax_gap = j["rsec_ax_gap"];
  scanner_geometry.rsec_tx_gap = j["rsec_tx_gap"];
  scanner_geometry.max_d_ring = j["max_d_ring"];
  scanner_geometry.number_of_tof_bins = j["number_of_tof_bins"];
  scanner_geometry.number_of_energy_bins = j["number_of_energy_bins"];
  scanner_geometry.radius = j["radius"];
  scanner_geometry.tx_virtual_crystal_num = j["tx_virtual_crystal_num"];
  scanner_geometry.ax_virtual_crystal_num = j["ax_virtual_crystal_num"];
  scanner_geometry.tx_phys_crystal_num = j["tx_phys_crystal_num"];
  scanner_geometry.ax_phys_crystal_num = j["ax_phys_crystal_num"];
  scanner_geometry.detector_x_dim = j["detector_x_dim"];
  scanner_geometry.detector_y_dim = j["detector_y_dim"];
  scanner_geometry.detector_z_dim = j["detector_z_dim"];
  scanner_geometry.energy_LLD = j["energy_LLD"];
  scanner_geometry.energy_ULD = j["energy_ULD"];
  scanner_geometry.EnergyResolutionAt511 = j["EnergyResolutionAt511"];
  scanner_geometry.TOF_resolution = j["TOF_resolution"];
  scanner_geometry.LM_TimeBlockDuration = j["LM_TimeBlockDuration"];

  scanner_geometry.ArcLength = scanner_geometry.s_width * scanner_geometry.detector_y_dim / 2.0f;
  scanner_geometry.TxFOV = 2 * scanner_geometry.radius * sin (scanner_geometry.ArcLength / (2 * scanner_geometry.radius) );
  scanner_geometry.TxFOV_TOF = scanner_geometry.TxFOV + 0.3*scanner_geometry.TOF_resolution;
  scanner_geometry.module_axial_pitch = scanner_geometry.n_cry_z * scanner_geometry.detector_z_dim + (scanner_geometry.n_cry_z - 1) * scanner_geometry.cry_ax_gap;

  return scanner_geometry;
}

void usage()
{
  std::cout << "Usage: root_to_petsird [options]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -r, --root-prefix <root_prefix>             Prefix of root files" << std::endl;
  std::cout << "  -s, --scanner-geometry-file <filename>      Scanner geometry file" << std::endl;
  std::cout << "  -p, --petsird-file <filename>               PETSiRD file" << std::endl;
  std::cout << "  -n, --number-of-root-files <number>         Number of root files" << std::endl;
  std::cout << "  -v, --verbose                               Verbose output" << std::endl;
  std::cout << "  -h, --help                                  Print this help message" << std::endl;
}

int calculate_element_index(int module_id, int submodule_id, int crystal_id, const ScannerGeometry& scannerGeometry)
{
  //int N_DET = scannerGeometry.n_det;
  //int N_MOD_xy = scannerGeometry.n_mod_xy;
  //int N_MOD_z = scannerGeometry.n_mod_z;
  int N_SMOD_xy = scannerGeometry.n_smod_xy;
  int N_SMOD_z = scannerGeometry.n_smod_z;
  int N_CRY_xy = scannerGeometry.n_cry_xy;
  int N_CRY_z = scannerGeometry.n_cry_z;
  // int N_CRY_layers = scannerGeometry.n_cry_layers;

  return (Int_t)(((module_id*N_SMOD_xy*N_SMOD_z)
                  + submodule_id)*N_CRY_xy*N_CRY_z
                 + crystal_id);
}

int calculate_module_index(int gantry_id, int rsector_id, const ScannerGeometry& scannerGeometry)
{
  const int N_RSEC_xy = scannerGeometry.n_rsec_xy;
  const int N_RSEC_z = scannerGeometry.n_rsec_z;

  return (Int_t)(gantry_id)*N_RSEC_z*N_RSEC_xy + (Int_t)rsector_id;
}

//! return a cuboid volume
petsird::BoxSolidVolume
get_crystal(ScannerGeometry& scannerGeometry)
{
  using petsird::Coordinate;
  petsird::BoxShape crystal_shape{ Coordinate{ { 0, 0, 0 } },
                                   Coordinate{ { 0, 0, scannerGeometry.detector_z_dim } },
                                   Coordinate{ { 0, scannerGeometry.detector_y_dim, scannerGeometry.detector_z_dim } },
                                   Coordinate{ { 0, scannerGeometry.detector_y_dim, 0 } },
                                   Coordinate{ { scannerGeometry.detector_x_dim, 0, 0 } },
                                   Coordinate{ { scannerGeometry.detector_x_dim, 0, scannerGeometry.detector_z_dim } },
                                   Coordinate{ { scannerGeometry.detector_x_dim, scannerGeometry.detector_y_dim, scannerGeometry.detector_z_dim } },
                                   Coordinate{ { scannerGeometry.detector_x_dim, scannerGeometry.detector_y_dim, 0 } } };

  petsird::BoxSolidVolume crystal{ crystal_shape, /* material_id */ 1 };
  return crystal;
}

//! return a module of NUM_CRYSTALS_PER_MODULE cuboids
petsird::DetectorModule
get_detector_module(ScannerGeometry& scannerGeometry)
{
  petsird::ReplicatedBoxSolidVolume rep_volume;
  {
    rep_volume.object = get_crystal(scannerGeometry);
    for (int rep_mod_z = 0; rep_mod_z < scannerGeometry.n_mod_z; ++rep_mod_z)
      for (int rep_mod_xy = 0; rep_mod_xy < scannerGeometry.n_mod_xy; ++rep_mod_xy)
        for (int rep_smod_z = 0; rep_smod_z < scannerGeometry.n_smod_z; ++rep_smod_z)
          for (int rep_smod_xy = 0; rep_smod_xy < scannerGeometry.n_smod_xy; ++rep_smod_xy)
            for (int rep_cry_z = 0; rep_cry_z < scannerGeometry.n_cry_z; ++rep_cry_z)
              for (int rep_cry_xy = 0; rep_cry_xy < scannerGeometry.n_cry_xy; ++rep_cry_xy)
                  for (int rep_cry_layer = 0; rep_cry_layer < scannerGeometry.n_cry_layers; ++rep_cry_layer)
                    {
                      petsird::RigidTransformation transform{ { { 1.0, 0.0, 0.0, scannerGeometry.radius + rep_cry_layer * scannerGeometry.detector_x_dim },
                                                                { 0.0, 1.0, 0.0, (rep_mod_xy - scannerGeometry.n_mod_xy / 2) * scannerGeometry.n_smod_xy * scannerGeometry.n_cry_xy * scannerGeometry.detector_y_dim
                                                                               + (rep_smod_xy - scannerGeometry.n_smod_xy / 2) * scannerGeometry.n_cry_xy * scannerGeometry.detector_y_dim
                                                                               + (rep_cry_xy - scannerGeometry.n_cry_xy / 2) * scannerGeometry.detector_y_dim  },
                                                                { 0.0, 0.0, 1.0, (rep_mod_z - scannerGeometry.n_mod_z / 2) * scannerGeometry.n_smod_z * scannerGeometry.n_cry_z * scannerGeometry.detector_z_dim
                                                                               + (rep_smod_z - scannerGeometry.n_smod_z / 2) * scannerGeometry.n_cry_z * scannerGeometry.detector_z_dim
                                                                               + (rep_cry_z - scannerGeometry.n_cry_z / 2) * scannerGeometry.detector_z_dim } } };
                      rep_volume.transforms.push_back(transform);
                    }
  }

  petsird::DetectorModule detector_module;
  detector_module.detecting_elements = rep_volume;

  return detector_module;
}


//! return scanner build by rotating a module around the (0,0,1) axis
petsird::ScannerGeometry
get_scanner_geometry(ScannerGeometry& scannerGeometry)
{
  petsird::ReplicatedDetectorModule rep_module;
  {
    rep_module.object = get_detector_module(scannerGeometry);
    std::vector<float> angles;
    for (int i = 0; i < scannerGeometry.n_rsec_xy; ++i)
      {
        angles.push_back(static_cast<float>((2 * M_PI * i) / scannerGeometry.n_rsec_xy));
      }

    for (auto angle : angles)
      for (int rep_rsec_z = 0; rep_rsec_z < scannerGeometry.n_rsec_z; ++rep_rsec_z)
        {
          petsird::RigidTransformation transform{ { { std::cos(angle), std::sin(angle), 0.F, 0.F },
                                                    { -std::sin(angle), std::cos(angle), 0.F, 0.F},
                                                    { 0.F, 0.F, 1.F, (rep_rsec_z - scannerGeometry.n_rsec_z / 2) * scannerGeometry.n_mod_z * scannerGeometry.n_smod_z * scannerGeometry.n_cry_z * scannerGeometry.detector_z_dim} } };

          rep_module.transforms.push_back(transform);
        }
  }
  petsird::ScannerGeometry scanner_geometry;
  scanner_geometry.replicated_modules.push_back(rep_module);
  return scanner_geometry;
}

// single ring as example
petsird::ScannerInformation
get_scanner_info(ScannerGeometry& scannerGeometry)
{
  // float radius = scannerGeometry.radius;
  // int n_detectors = scannerGeometry.n_det;
  // int n_rings = scannerGeometry.n_rings;
  unsigned long NUMBER_OF_TOF_BINS = static_cast<unsigned long>(scannerGeometry.number_of_tof_bins);
  unsigned long NUMBER_OF_EVENT_ENERGY_BINS = static_cast<unsigned long>(scannerGeometry.number_of_energy_bins);
  float energy_LLD = scannerGeometry.energy_LLD;
  float energy_ULD =scannerGeometry.energy_ULD;

  petsird::ScannerInformation scanner_info;
  scanner_info.model_name = "PETSIRD_GATEROOT"; // TODO

  const auto num_types_of_modules = 1;
  // Pre-allocate various structures to have the correct size for num_types_of_modules
  // (We will still have to set descent values into each of these.)
  petsird_helpers::create::initialize_scanner_information_dimensions(scanner_info, num_types_of_modules,
                                                                     /* allocate_detection_bin_efficiencies = */ false,
                                                                     /* allocate_module_pair_efficiencies = */ false);

  // TODO scanner_info.bulk_materials

  // geometry
  scanner_info.scanner_geometry = get_scanner_geometry(scannerGeometry);

  // TOF and energy information
  {
    auto& all_tof_bin_edges = scanner_info.tof_bin_edges;
    auto& all_tof_resolutions = scanner_info.tof_resolution;
    auto& all_event_energy_bin_edges = scanner_info.event_energy_bin_edges;
    auto& all_event_energy_resolutions = scanner_info.energy_resolution_at_511;

    // only 1 type of module in the current scanner
    assert(num_types_of_modules == 1);
    const petsird::TypeOfModule type_of_module{ 0 };

    typedef yardl::NDArray<float, 1> FArray1D;
    // TOF info (in mm)
    FArray1D tof_bin_edges_arr;
    yardl::resize(tof_bin_edges_arr, { NUMBER_OF_TOF_BINS + 1 });
    for (std::size_t i = 0; i < tof_bin_edges_arr.size(); ++i)
      tof_bin_edges_arr[i] = (i - NUMBER_OF_TOF_BINS / 2.F) / NUMBER_OF_TOF_BINS * 2 * scannerGeometry.TxFOV_TOF;
    const petsird::BinEdges tof_bin_edges{ tof_bin_edges_arr };
    all_tof_bin_edges[type_of_module][type_of_module] = tof_bin_edges;

    // TODO use speed-of-light here
    all_tof_resolutions[type_of_module][type_of_module] = scannerGeometry.TOF_resolution*0.3; // conversion from psec to mm (e.g. 200ps TOF is equivalent to 60mm uncertainty)

    FArray1D event_energy_bin_edges_arr;
    yardl::resize(event_energy_bin_edges_arr, { NUMBER_OF_EVENT_ENERGY_BINS + 1 });
    for (std::size_t i = 0; i < event_energy_bin_edges_arr.size(); ++i)
      event_energy_bin_edges_arr[i] = energy_LLD + i * (energy_ULD - energy_LLD) / NUMBER_OF_EVENT_ENERGY_BINS;
    petsird::BinEdges event_energy_bin_edges{ event_energy_bin_edges_arr };
    all_event_energy_bin_edges[type_of_module] = event_energy_bin_edges;
    all_event_energy_resolutions[type_of_module] = scannerGeometry.EnergyResolutionAt511;    // as fraction of 511 (e.g. 0.11F)
  }

  // TODO scanner_info.coincidence_policy = petsird::CoincidencePolicy::kRejectMultiples;
  scanner_info.delayed_coincidences_are_stored = false;
  scanner_info.triple_events_are_stored = false;
  return scanner_info;
}



uint32_t tofToIdx(double delta_time_psec, const petsird::ScannerInformation& scanner_info)
{
  constexpr petsird::TypeOfModule type_of_module{ 0 };
  const auto& tof_bin_edges = scanner_info.tof_bin_edges[type_of_module][type_of_module].edges;

  float tofPos_mm = delta_time_psec * 0.15; //conversion from time difference (in psec) to spatial position in LOR (in mm) DT*C/2
  for (size_t i = 0; i < tof_bin_edges.size() - 1; ++i)
  {
    if (tofPos_mm >= tof_bin_edges[i] && tofPos_mm < tof_bin_edges[i+1])
    {
      return static_cast<uint32_t>(i);
    }
  }
  std::cout << "WARNING: TOF out of range: " << tofPos_mm << std::endl;
  //std::stringstream ss;
  //ss << "WARNING: TOF out of range: " << tofPos_mm;
  throw std::runtime_error("TOF out of range");
}

uint32_t energyToIdx(float energy, const petsird::ScannerInformation& scanner_info)
{
  constexpr petsird::TypeOfModule type_of_module{ 0 };
  const auto& energy_bin_edges = scanner_info.event_energy_bin_edges[type_of_module].edges;
  for (size_t i = 0; i < energy_bin_edges.size() - 1; ++i)
  {
    if (energy >= energy_bin_edges[i] && energy < energy_bin_edges[i+1])
    {
      return static_cast<uint32_t>(i);
    }
  }
  std::stringstream ss;
  ss << "WARNING: Energy out of range: " << energy;
  throw std::runtime_error(ss.str());
}

petsird::Coordinate mean_position(const petsird::BoxShape& box_shape)
{
  petsird::Coordinate mean;
  mean.c = {0, 0, 0};
  for (auto& corner : box_shape.corners)
    {
      mean.c += corner.c;
    }
  mean.c /= box_shape.corners.size();
  return mean;
}

auto create_new_time_block(const petsird::ScannerInformation& scanner)
{
  const auto num_module_types = scanner.scanner_geometry.replicated_modules.size();
  petsird::EventTimeBlock time_block;
  // allocate lists for prompt_events
  // warning: would have to do the same for delayeds once we store them
  time_block.prompt_events = petsird_helpers::create::construct_2D_nested_vector<petsird::ListOfCoincidenceEvents>(num_module_types, num_module_types);
  return time_block;
}

int main(int argc, char** argv)
{
  std::string root_prefix  = std::string{};
  std::string scanner_geometry_file = std::string{};
  std::string petsird_file = std::string{};
  int number_of_root_files = 2;
  bool verbose = false;

  // Parse command line args:
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-r" || arg == "--root-prefix") {
      root_prefix = argv[++i];
    } else if (arg == "-s" || arg == "--scanner-geometry-file") {
      scanner_geometry_file = argv[++i];
    } else if (arg == "-p" || arg == "--petsird-file") {
      petsird_file = argv[++i];
    } else if (arg == "-n" || arg == "--number-of-root-files") {
      number_of_root_files = atoi(argv[++i]);
    } else if (arg == "-v" || arg == "--verbose") {
      verbose = true;
    } else if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    } else {
      std::cerr << "Unknown argument: " << arg << std::endl;
      return 1;
    }
  }

  if (root_prefix.empty()) {
    std::cerr << "Missing root prefix" << std::endl;
    usage();
    return 1;
  }

  if (petsird_file.empty()) {
    std::cerr << "Missing petsird file" << std::endl;
    usage();
    return 1;
  }

  // Print arguments and exit
  std::cout << "root_prefix: " << root_prefix << std::endl;
  std::cout << "scanner_geometry_file: " << scanner_geometry_file << std::endl;
  std::cout << "petsird_file: " << petsird_file << std::endl;

  // Read scanner geometry
  ScannerGeometry scannerGeometry;
  if (scanner_geometry_file.empty()) {
    std::cout << "Using default scanner geometry" << std::endl;
    return 1;
  } else {
    scannerGeometry = ReadScannerGeometry(scanner_geometry_file);
  }

  string filedir, inputfilename;
  Int_t   Trues = 0, Scatters = 0, Randoms = 0;

  //####################################################################
  //#             Declaration of leaves types - TTree Coincidences     #
  //####################################################################
  Float_t         			axialPos, rotationAngle, sinogramS, sinogramTheta;
  Char_t          			comptVolName1[255], comptVolName2[255];
  Int_t           			compton1, compton2, gantryID1, gantryID2;
  Int_t           			runID, sourceID1, sourceID2, eventID1, eventID2;
  Int_t           			layerID1, layerID2, crystalID1, crystalID2;
  Int_t           			submoduleID1, submoduleID2, moduleID1, moduleID2, rsectorID1, rsectorID2;
  Int_t           			comptonPhantom1, comptonPhantom2;
  Float_t         			energy1, energy2; //in MeV
  Float_t         			globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2; //in mm
  Float_t         			sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2; //in mm
  Double_t        			time1, time2; //in sec
  unsigned long long int   	nentries;

  //######################################################################################
  //#                        Set branch addresses - TTree Coincidences                   #
  //######################################################################################

  filedir = root_prefix;
  TChain *Coincidences = new TChain("Coincidences");

  for (int i = 0; i < number_of_root_files; i++) {
    std::ostringstream fileNumber;
    fileNumber << i + 1;
    inputfilename = filedir + fileNumber.str() + ".root";
    std::cout << "Input file name is " << inputfilename << std::endl;
    Coincidences->Add(inputfilename.c_str());
  }


  Coincidences->SetBranchStatus("*",0);
  Coincidences->SetBranchAddress("axialPos",&axialPos);
  Coincidences->SetBranchAddress("comptVolName1",&comptVolName1);
  Coincidences->SetBranchAddress("comptVolName2",&comptVolName2);
  Coincidences->SetBranchAddress("comptonCrystal1",&compton1);
  Coincidences->SetBranchAddress("comptonCrystal2",&compton2);
  Coincidences->SetBranchAddress("crystalID1",&crystalID1);
  Coincidences->SetBranchAddress("crystalID2",&crystalID2);
  Coincidences->SetBranchAddress("comptonPhantom1",&comptonPhantom1);
  Coincidences->SetBranchAddress("comptonPhantom2",&comptonPhantom2);
  Coincidences->SetBranchAddress("energy1",&energy1);
  Coincidences->SetBranchAddress("energy2",&energy2);
  Coincidences->SetBranchAddress("eventID1",&eventID1);
  Coincidences->SetBranchAddress("eventID2",&eventID2);
  Coincidences->SetBranchAddress("globalPosX1",&globalPosX1);
  Coincidences->SetBranchAddress("globalPosX2",&globalPosX2);
  Coincidences->SetBranchAddress("globalPosY1",&globalPosY1);
  Coincidences->SetBranchAddress("globalPosY2",&globalPosY2);
  Coincidences->SetBranchAddress("globalPosZ1",&globalPosZ1);
  Coincidences->SetBranchAddress("globalPosZ2",&globalPosZ2);
  Coincidences->SetBranchAddress("layerID1",&layerID1);
  Coincidences->SetBranchAddress("layerID2",&layerID2);
  Coincidences->SetBranchAddress("moduleID1",&moduleID1);
  Coincidences->SetBranchAddress("moduleID2",&moduleID2);
  Coincidences->SetBranchAddress("rotationAngle",&rotationAngle);
  Coincidences->SetBranchAddress("rsectorID1",&rsectorID1);
  Coincidences->SetBranchAddress("rsectorID2",&rsectorID2);
  Coincidences->SetBranchAddress("runID",&runID);
  Coincidences->SetBranchAddress("sinogramS",&sinogramS);
  Coincidences->SetBranchAddress("sinogramTheta",&sinogramTheta);
  Coincidences->SetBranchAddress("sourceID1",&sourceID1);
  Coincidences->SetBranchAddress("sourceID2",&sourceID2);
  Coincidences->SetBranchAddress("sourcePosX1",&sourcePosX1);
  Coincidences->SetBranchAddress("sourcePosX2",&sourcePosX2);
  Coincidences->SetBranchAddress("sourcePosY1",&sourcePosY1);
  Coincidences->SetBranchAddress("sourcePosY2",&sourcePosY2);
  Coincidences->SetBranchAddress("sourcePosZ1",&sourcePosZ1);
  Coincidences->SetBranchAddress("sourcePosZ2",&sourcePosZ2);
  Coincidences->SetBranchAddress("submoduleID1",&submoduleID1);
  Coincidences->SetBranchAddress("submoduleID2",&submoduleID2);
  Coincidences->SetBranchAddress("time1",&time1);
  Coincidences->SetBranchAddress("time2",&time2);
  Coincidences->SetBranchAddress("gantryID1",&gantryID1);
  Coincidences->SetBranchAddress("gantryID2",&gantryID2);

  nentries = (unsigned long long int)(Coincidences->GetEntries());

  printf("Total Number of Coincidence Events in the ROOT file:= %llu \n",nentries );

  // Output PETSIRD
  petsird::Header header;
  header.scanner = get_scanner_info(scannerGeometry);
  auto& scanner = header.scanner;

  if (verbose) {
    // Print scanner information
    std::cout << "Scanner information:" << std::endl;
    constexpr petsird::TypeOfModule type_of_module{ 0 };
    const auto& tof_bin_edges = header.scanner.tof_bin_edges[type_of_module][type_of_module];
    const auto num_tof_bins = tof_bin_edges.NumberOfBins();
    std::cout << "Number of TOF bins: " << num_tof_bins << std::endl;
    std::cout << "TOF bin edges: " << tof_bin_edges.edges << std::endl;
    const auto& event_energy_bin_edges = header.scanner.event_energy_bin_edges[type_of_module];
    const auto num_event_energy_bins = event_energy_bin_edges.NumberOfBins();
    std::cout << "Number of energy bins: " << num_event_energy_bins << std::endl;
    std::cout << "Event energy bin edges: " << event_energy_bin_edges.edges << std::endl;
  }

  // Write PETSiRD file
  petsird::binary::PETSIRDWriter writer(petsird_file);
  writer.WriteHeader(header);

  long current_time_block = -1;
  auto time_block = create_new_time_block(header.scanner);

  const auto event_time_block_duration = scannerGeometry.LM_TimeBlockDuration; // ms
  const petsird::TypeOfModule type_of_module{ 0 };
  unsigned long Counts_binned = 0;
  for (unsigned long long int i = 0 ; i < nentries ; i++)
  {
    if (i % 1000000 == 0) {
      printf("Processing event %llu of %llu, (%f percent)\n", i, nentries, 100.0f*i/nentries);
    }

    Coincidences->GetEntry(i);
    if (eventID1 == eventID2)
    {
	    if (comptonPhantom1 == 0 && comptonPhantom2 == 0) {
        petsird::CoincidenceEvent event;
        petsird::ExpandedDetectionBin expanded_detection_bin;
        expanded_detection_bin.module_index = calculate_module_index(gantryID1, rsectorID1, scannerGeometry);
        expanded_detection_bin.element_index = calculate_element_index(moduleID1, submoduleID1, crystalID1, scannerGeometry);
        expanded_detection_bin.energy_index = static_cast<uint32_t>(energyToIdx(1.0e3*energy1, scanner));
        event.detection_bins[0] = petsird_helpers::make_detection_bin(header.scanner, type_of_module, expanded_detection_bin);
        expanded_detection_bin.module_index = calculate_module_index(gantryID2, rsectorID2, scannerGeometry);
        expanded_detection_bin.element_index = calculate_element_index(moduleID2, submoduleID2, crystalID2, scannerGeometry);
        expanded_detection_bin.energy_index = static_cast<uint32_t>(energyToIdx(1.0e3*energy2, scanner));
        event.detection_bins[1] = petsird_helpers::make_detection_bin(header.scanner, type_of_module, expanded_detection_bin);
        double dt_psec = 1.0e12f*(time1 - time2); //in psec
        if (abs(dt_psec) > scannerGeometry.TxFOV_TOF/0.3f) {
          continue;
        }
        event.tof_idx = static_cast<uint32_t>(tofToIdx(dt_psec, scanner));

        if (verbose && i%100000 == 0) {
          std::cout << "Event " << i << std::endl;
          const auto expanded_detection_bin0
            = petsird_helpers::expand_detection_bin(header.scanner, type_of_module, event.detection_bins[0]);
          const auto expanded_detection_bin1
            = petsird_helpers::expand_detection_bin(header.scanner, type_of_module, event.detection_bins[1]);
          std::cout << "    "
                    << "[ExpandedDetectionBin(module=" << expanded_detection_bin0.module_index << ", "
                    << "el=" << expanded_detection_bin0.element_index << ", "
                    << "energy_index=" << expanded_detection_bin0.energy_index
                    << "), ExpandedDetectionBin(module=" << expanded_detection_bin1.module_index << ", "
                    << "el=" << expanded_detection_bin1.element_index << ", "
                    << "energy_index=" << expanded_detection_bin1.energy_index << ")]\n";
          std::cout << "  tof_idx: " << event.tof_idx << std::endl;
          std::cout << "gantryID1=" << gantryID1 << ", rsectorID1=" << rsectorID1
                    << ", moduleID1=" << moduleID1 << ", submoduleID1=" << submoduleID1 << ", crystalID1=" << crystalID1 << '\n';
          std::cout << "gantryID2=" << gantryID2 << ", rsectorID2=" << rsectorID2
                    << ", moduleID2=" << moduleID2 << ", submoduleID2=" << submoduleID2 << ", crystalID2=" << crystalID2 << '\n';
          const auto box_shape0 = petsird_helpers::geometry::get_detecting_box(header.scanner, type_of_module, expanded_detection_bin0);
          const auto mean_pos0 = mean_position(box_shape0);
          const auto box_shape1 = petsird_helpers::geometry::get_detecting_box(header.scanner, type_of_module, expanded_detection_bin1);
          const auto mean_pos1 = mean_position(box_shape1);
          std::cout << "  pos 1          : " << mean_pos0.c[0] << ", " << mean_pos0.c[1] << ", " << mean_pos0.c[2] << "\n";
          std::cout << "  pos 2          : " << mean_pos1.c[0] << ", " << mean_pos1.c[1] << ", " << mean_pos1.c[2] << "\n";

          std::cout << "  GlobalPosition 1: " << globalPosX1 << ", " << globalPosY1 << ", " << globalPosZ1 << std::endl;
          float distance_1 = std::sqrt(std::pow(mean_pos0.c[0]-globalPosX1, 2) + std::pow(mean_pos0.c[1]+globalPosY1, 2) + std::pow(mean_pos0.c[2]-globalPosZ1, 2));
          std::cout << "  Distance 1: " << distance_1 << std::endl;
          std::cout << "  GlobalPosition 2: " << globalPosX2 << ", " << globalPosY2 << ", " << globalPosZ2 << std::endl;
          float distance_2 = std::sqrt(std::pow(mean_pos1.c[0]-globalPosX2, 2) + std::pow(mean_pos1.c[1]+globalPosY2, 2) + std::pow(mean_pos1.c[2]-globalPosZ2, 2));
          std::cout << "  Distance 2: " << distance_2 << std::endl;
        }
        long this_time_block = static_cast<long>(time1*1.0e3 / event_time_block_duration);
        if (this_time_block != current_time_block) {
          if (current_time_block != -1) {
            writer.WriteTimeBlocks(time_block);
          }
          current_time_block = this_time_block;
          time_block = create_new_time_block(header.scanner);
          time_block.time_interval.start = time1*1.0e3;
          time_block.time_interval.stop = time1*1.0e3 + event_time_block_duration;
        }
        time_block.prompt_events[type_of_module][type_of_module].push_back(event);
        Counts_binned++;
        Trues++;
      } else {
        Scatters++;
      }
	  } else {
      Randoms++;
    }
  }
  writer.WriteTimeBlocks(time_block);
  writer.EndTimeBlocks();
  writer.Close();

  printf("Total Number of Coincidence Events in the ROOT file:= %llu ...\n",nentries );
  printf("Total Number of Coincidence Events registered in list-mode or sinogram format:= %lu ...\n", Counts_binned);  return(0);
}
