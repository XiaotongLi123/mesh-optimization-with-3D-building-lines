cmake_minimum_required(VERSION 3.10)

project(mesh_optimization)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 包含头文件目录
include_directories(${CMAKE_SOURCE_DIR}/include)

# 查找 Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# 查找 zlib
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

# 查找 PCL
find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
add_definitions(${PCL_DEFINITIONS})

# 收集源文件
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/include/*.h")

add_executable(mesh_optimization ${SOURCES} ${HEADERS})

# 链接 PCL、Eigen3 和 zlib
target_link_libraries(mesh_optimization ${PCL_LIBRARIES} Eigen3::Eigen ${ZLIB_LIBRARIES})
