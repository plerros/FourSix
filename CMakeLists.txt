cmake_minimum_required(VERSION 2.8.10)

project(HalfTurn)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

add_executable(HalfTurn
	src/data/data_in.cpp
	src/data/data.cpp
	src/data/data_out.cpp
	src/json/json.cpp
	src/main.cpp
	src/triangulation/triangulation.cpp
	)

target_include_directories(HalfTurn PRIVATE
	src
	src/data
	src/json
	src/triangulation
	)

target_link_libraries(HalfTurn
	CGAL::CGAL
	CGAL::CGAL_Basic_viewer
	)
