cmake_minimum_required(VERSION 2.8.3)
project(mpl_external_planner)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-declarations")

set(MOTION_PRIMITIVE_LIBRARY_INCLUDE_DIRS "../motion_primitive_library/include")
set(DECOMP_UTIL_INCLUDE_DIRS "../DecompROS/DecompUtil/include")

set(MOTION_PRIMITIVE_LIBRARY_LIBRARIES
  "/home/gaussian/cmu_ri_phd/phd_research/csbpl_ros/devel/lib/libpoly_solver.so"
  "/home/gaussian/cmu_ri_phd/phd_research/csbpl_ros/devel/lib/libmap_planner.so")

find_package(catkin_simple REQUIRED)
find_package(cmake_modules REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(motion_primitive_library REQUIRED)
#find_package(decomp_util REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR} ${MOTION_PRIMITIVE_LIBRARY_INCLUDE_DIRS} ${DECOMP_UTIL_INCLUDE_DIRS})

catkin_simple(ALL_DEPS_REQUIRED)

cs_add_executable(test_primitive_collide.cpp src/test_primitive_collide.cpp)

cs_install()

cs_export()
