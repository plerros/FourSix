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

		auto positions = input.get_points();
		std::vector<Point> points;
		for (const auto& pos : positions)
			points.push_back(Point(pos.first, pos.second));

		auto constraints = input.get_constraints();
		CDT cdt;
		for (const auto& point : points)
			cdt.insert(point);
		for (const auto& constraint : constraints)
			cdt.insert_constraint(points[constraint.first], points[constraint.second]);

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
