cmake_minimum_required (VERSION 2.8.8)

# Fallback for using newer policies on CMake < 2.8
if(${CMAKE_VERSION} VERSION_LESS 2.8)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# Project name
project(aa2cg)

# Add the libraries and packages
add_subdirectory(lib/fmt-6.2.1)
find_package(OpenMP REQUIRED)
find_package(Eigen3 3.3 NO_MODULE)

# Set the standard version of C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wall -fconcepts")

# Add the execuatbles,path to source codes, and dependencies
set(COMMON_SRC src/inputs.cpp src/mapping_matrix.cpp src/lammps_data.cpp
    src/lammps_dump.cpp src/blocks.cpp src/reorder.cpp 
    src/cg_data.cpp src/cg_dump.cpp src/outputs.cpp)
add_executable(aa2cg src/main.cpp)
add_library(common ${COMMON_SRC})

# Directories to the required files and libraries
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/lib/fmt-6.2.1/include)
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(/packages/7x/eigen/3.3.7-gcc-stock/include/eigen3)
include_directories(/usr/include/eigen3)
include_directories(/usr/local/include/eigen3)

# Linking the libraries
target_link_libraries(aa2cg common fmt)
