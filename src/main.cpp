#include <boost/json.hpp>
#include <boost/json/src.hpp>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <iostream>

#include "json.hpp"
#include "data_in.hpp"

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cerr <<
			"Usage: pretty <filename>"
			<< std::endl;
		return EXIT_FAILURE;
	}

	try {
		// Parse the file as JSON
		auto const jv = parse_file(argv[1]);

		// Now pretty-print the value
		//pretty_print(std::cout, jv);


		data_in input{jv};
		//input.print();

		auto points = input.get_points();
		auto constraints = input.get_constraints();
		CDT cdt;
		for (auto it = points.begin(); it != points.end(); it++)
			cdt.insert(Point(it[0].first, it[0].second));
		
		for (auto it = constraints.begin(); it != constraints.end(); it++) {
			std::array<std::array<uint64_t, 2>, 2> arr;
			arr[0][0] = points[it[0].first].first;
			arr[0][1] = points[it[0].first].second;
			arr[1][0] = points[it[0].second].first;
			arr[1][1] = points[it[0].second].second;

			cdt.insert_constraint(Point(arr[0][0], arr[0][1]), Point(arr[1][0], arr[1][1]));
		}
		CGAL::draw(cdt);
	}
	catch(std::exception const& e) {
		std::cerr <<
			"Caught exception: "
			<< e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
