cmake_minimum_required(VERSION 2.8.10)

project(a_cute_triangle)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

add_executable(a_cute_triangle
	src/data/data_in.cpp
	src/data/data.cpp
	src/data/data_out.cpp
	src/json/json.cpp
	src/main.cpp
	src/triangulation/triangulation.cpp
	)

target_include_directories(a_cute_triangle PRIVATE
	src
	src/data
	src/json
	src/triangulation
	)

target_link_libraries(a_cute_triangle
	CGAL::CGAL
	CGAL::CGAL_Basic_viewer
	)
