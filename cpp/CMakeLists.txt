cmake_minimum_required(VERSION 3.20.0) # older would work, but could give warnings on policy CMP0074
project(root_to_petsird VERSION 0.5.0)

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

find_package(PETSIRD 0.13 REQUIRED CONFIG)
find_package(ROOT REQUIRED CONFIG COMPONENTS Tree)

add_executable(root_to_petsird main.cpp)
# needed for helpers
target_link_libraries(root_to_petsird PUBLIC PETSIRD::helpers)
target_link_libraries(root_to_petsird PUBLIC ROOT::Tree)
install(TARGETS root_to_petsird DESTINATION bin)
