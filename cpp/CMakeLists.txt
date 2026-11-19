cmake_minimum_required(VERSION 3.20.0) # older would work, but could give warnings on policy CMP0074
project(root_to_petsird VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)

#Set the build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif ()

if(WIN32)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/W3> $<$<COMPILE_LANGUAGE:CXX>:/WX>)
else()
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wall> $<$<COMPILE_LANGUAGE:CXX>:-Wextra> $<$<COMPILE_LANGUAGE:CXX>:-pedantic>)
endif()

find_package(PETSIRD 0.13.3 REQUIRED CONFIG)
find_package(ROOT REQUIRED CONFIG COMPONENTS Tree)

find_package(Threads REQUIRED)
include(CheckLanguage)
set(GATEROOT_CUDA "AUTO" CACHE STRING "CUDA event processing: AUTO, ON, OFF")
set_property(CACHE GATEROOT_CUDA PROPERTY STRINGS AUTO ON OFF)
set(EVENT_BATCH_SOURCE event_batch_cpu.cpp)
if(NOT GATEROOT_CUDA STREQUAL "OFF")
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(EVENT_BATCH_SOURCE event_batch.cu)
  elseif(GATEROOT_CUDA STREQUAL "ON")
    message(FATAL_ERROR "CUDA requested but nvcc was not found")
  endif()
endif()
add_executable(root_to_petsird main.cpp ${EVENT_BATCH_SOURCE})
set_target_properties(root_to_petsird PROPERTIES CUDA_STANDARD 17)
if(CMAKE_CUDA_COMPILER AND NOT GATEROOT_CUDA STREQUAL "OFF")
  target_compile_options(root_to_petsird PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>)
endif()
target_link_libraries(root_to_petsird PRIVATE Threads::Threads)
# needed for helpers
target_link_libraries(root_to_petsird PUBLIC PETSIRD::helpers)
target_link_libraries(root_to_petsird PUBLIC ROOT::Tree)
install(TARGETS root_to_petsird DESTINATION bin)
