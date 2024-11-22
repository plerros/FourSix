cmake_minimum_required(VERSION 2.8.10)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_C_FLAGS_RELEASE "-O2" CACHE STRING "Flags used by the compiler during release builds." FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O2" CACHE STRING "Flags used by the compiler during release builds." FORCE)
set(CMAKE_C_FLAGS_MINSIZEREL "-Os" CACHE STRING "Flags used by the compiler during release builds." FORCE)
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os" CACHE STRING "Flags used by the compiler during release builds." FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "Flags used by the compiler during release builds." FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "Flags used by the compiler during release builds." FORCE)

set(CMAKE_BUILD_PARALLEL_LEVEL $(nproc))

project(FourSix)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

add_executable(FourSix
	src/data/data_in.cpp
	src/data/data.cpp
	src/data/data_out.cpp
	src/json/json.cpp
	src/helper/helper.cpp
	src/main.cpp
	src/obtuse_polygon/obtuse_polygon.cpp
	src/optimization_methods/ant_colony.cpp
	src/optimization_methods/local_search.cpp
	src/optimization_methods/mixed_recursive.cpp
	src/optimization_methods/simulated_annealing.cpp
	src/triangulation/triangulation.cpp
	)

target_include_directories(FourSix PRIVATE
	src
	src/data
	src/json
	src/helper
	src/triangulation
	src/obtuse_polygon
	)

target_link_libraries(FourSix
	CGAL::CGAL
	CGAL::CGAL_Basic_viewer
	)
