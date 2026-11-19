#include "event_batch.h"
#include <stdexcept>
bool gpu_available() { return false; }
std::vector<BatchEvent> gpu_transform(const std::vector<BatchInput>&,
    const std::vector<float>&, double, double, double, double) {
  throw std::runtime_error("GPU mode requires a CUDA-enabled build and a visible NVIDIA GPU");
}
