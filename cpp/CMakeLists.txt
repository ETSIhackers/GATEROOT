cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(PETSIRDUseCaseTemplate VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

#Set the build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif ()

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

add_subdirectory(../PETSIRD/cpp/generated PETSIRD_generated)

find_package(ROOT REQUIRED COMPONENTS Core Imt RIO Net Hist Graf Graf3d Gpad Tree TreePlayer Rint Postscript Matrix Physics MathCore Thread MultiProc ROOTDataFrame)
add_executable(root_to_petsird main.cpp)
include_directories(root_to_petsird PUBLIC ../PETSIRD/cpp/generated)
target_link_libraries(root_to_petsird PUBLIC prd_generated)
target_link_libraries(root_to_petsird PUBLIC ROOT::Core ROOT::Imt ROOT::RIO ROOT::Net ROOT::Hist ROOT::Graf ROOT::Graf3d ROOT::Gpad ROOT::Tree ROOT::TreePlayer ROOT::Rint ROOT::Postscript ROOT::Matrix ROOT::Physics ROOT::MathCore ROOT::Thread ROOT::MultiProc ROOT::ROOTDataFrame)
install(TARGETS root_to_petsird DESTINATION bin)
