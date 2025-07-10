cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(root_to_petsird VERSION 0.3.0)

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

set(PETSIRD_dir ../PETSIRD/cpp/generated)
add_subdirectory(${PETSIRD_dir} PETSIRD_generated)

find_package(ROOT REQUIRED COMPONENTS Tree)

add_executable(root_to_petsird main.cpp)
target_include_directories(root_to_petsird PUBLIC ${PETSIRD_dir})
# needed for helpers
target_include_directories(root_to_petsird PUBLIC ${PETSIRD_dir}/..)
target_include_directories(root_to_petsird PUBLIC ${PETSIRD_dir}/../helpers/include)
target_link_libraries(root_to_petsird PUBLIC petsird_generated)
target_link_libraries(root_to_petsird PUBLIC ROOT::Tree)
install(TARGETS root_to_petsird DESTINATION bin)
