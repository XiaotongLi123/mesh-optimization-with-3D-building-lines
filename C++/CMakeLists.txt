cmake_minimum_required(VERSION 3.10)

project(mesh_optimization)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include header file directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Find Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find PCL
find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
add_definitions(${PCL_DEFINITIONS})

# Collect source files
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/include/*.h")

add_executable(mesh_optimization ${SOURCES} ${HEADERS})

# Link PCL and Eigen3
target_link_libraries(mesh_optimization ${PCL_LIBRARIES} Eigen3::Eigen)
