cmake_minimum_required(VERSION 2.8.10)

project(a_cute_triangle)

find_package(CGAL REQUIRED)

add_executable(a_cute_triangle
	src/main.cpp
	src/json/json.cpp
	src/data/data_in.cpp
	)

target_include_directories(a_cute_triangle PRIVATE
	src/json
	src/data
	)

target_link_libraries(a_cute_triangle
	CGAL::CGAL
	)
