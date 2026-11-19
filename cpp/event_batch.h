#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

struct BatchInput {
  double time_ms, dt_ps;
  double energy1, energy2;
  uint64_t entry;
  uint32_t a, b;
  bool selected;
};

struct BatchEvent {
  double time_ms;
  uint64_t entry;
  uint32_t a, b, tof;
  bool valid;
};

bool gpu_available();
std::vector<BatchEvent> gpu_transform(const std::vector<BatchInput>& input,
    const std::vector<float>& edges, double energy_low, double energy_high,
    double start_ms, double end_ms);

inline std::vector<BatchEvent> cpu_transform(const std::vector<BatchInput>& input,
    const std::vector<float>& edges, double energy_low, double energy_high,
    double start_ms, double end_ms) {
  std::vector<BatchEvent> out;
  out.reserve(input.size());
  for (const auto& raw : input) {
    double position = raw.dt_ps * .299792458 / 2;
    if (raw.a < raw.b) position = -position;
    auto k = std::upper_bound(edges.begin(), edges.end(), position) - edges.begin() - 1;
    if (!raw.selected || raw.a == raw.b || raw.time_ms < start_ms || raw.time_ms >= end_ms
        || raw.energy1 < energy_low || raw.energy1 >= energy_high
        || raw.energy2 < energy_low || raw.energy2 >= energy_high
        || k < 0 || k >= static_cast<decltype(k)>(edges.size()-1)) continue;
    out.push_back({raw.time_ms, raw.entry, std::max(raw.a,raw.b), std::min(raw.a,raw.b),
                   static_cast<uint32_t>(k), true});
  }
  std::stable_sort(out.begin(), out.end(), [](const auto& a,const auto& b) {
    return a.time_ms < b.time_ms;
  });
  return out;
}
