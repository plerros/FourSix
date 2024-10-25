cmake_minimum_required(VERSION 2.8.10)

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
