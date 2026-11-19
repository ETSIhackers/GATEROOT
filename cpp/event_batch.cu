#include "event_batch.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <mutex>
#include <stdexcept>

bool gpu_available() {
  int count = 0;
  return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

__global__ void transform_kernel(const BatchInput* input, BatchEvent* output, size_t count,
    const float* edges, int bins, double low, double high, double start, double stop) {
  const size_t i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i >= count) return;
  const auto raw = input[i];
  double pos = raw.dt_ps * .299792458 / 2;
  if (raw.a < raw.b) pos = -pos;
  int left = 0, right = bins+1;
  while (left < right) {
    int mid = (left+right)/2;
    if (pos < static_cast<double>(edges[mid])) right = mid;
    else left = mid+1;
  }
  const int k = left-1;
  const bool valid = raw.selected && raw.a != raw.b && raw.time_ms >= start && raw.time_ms < stop
      && raw.energy1 >= low && raw.energy1 < high && raw.energy2 >= low && raw.energy2 < high
      && k >= 0 && k < bins;
  output[i] = {raw.time_ms,raw.entry,raw.a>raw.b?raw.a:raw.b,raw.a>raw.b?raw.b:raw.a,
                static_cast<uint32_t>(k),valid};
}

struct EventLess {
  __host__ __device__ bool operator()(const BatchEvent& a,const BatchEvent& b) const {
    if (a.valid != b.valid) return a.valid;
    return a.time_ms < b.time_ms;
  }
};

std::vector<BatchEvent> gpu_transform(const std::vector<BatchInput>& input,
    const std::vector<float>& edges, double low, double high, double start, double stop) {
  // Independent CPU readers overlap; serialize device work to bound device memory.
  static std::mutex gpu_mutex;
  const std::lock_guard<std::mutex> lock(gpu_mutex);
  if (input.empty()) return {};
  thrust::device_vector<BatchInput> raw(input.begin(),input.end());
  thrust::device_vector<float> axis(edges.begin(),edges.end());
  thrust::device_vector<BatchEvent> transformed(input.size());
  transform_kernel<<<(input.size()+255)/256,256>>>(thrust::raw_pointer_cast(raw.data()),
      thrust::raw_pointer_cast(transformed.data()), input.size(), thrust::raw_pointer_cast(axis.data()),
      static_cast<int>(edges.size()-1),low,high,start,stop);
  if (cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess)
    throw std::runtime_error("CUDA event transformation failed");
  thrust::stable_sort(transformed.begin(), transformed.end(), EventLess{});
  std::vector<BatchEvent> result(input.size());
  thrust::copy(transformed.begin(), transformed.end(),result.begin());
  result.erase(std::find_if(result.begin(), result.end(), [](const auto& e){return !e.valid;}), result.end());
  return result;
}
