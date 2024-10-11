#include <boost/json.hpp>
#include <boost/json/src.hpp>
#include <iostream>

#include "configuration.hpp"

#include "json.hpp"
#include "data_in.hpp"
#include "data.hpp"

size_t count_obtuse(CDT *cdt)
{
	size_t ret = 0;
	for (auto it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {
		auto triangle = cdt->triangle(it);
		std::array <CGAL::Vector_2<K>, 3> edge;

		for (size_t i = 0; i < 3; i++) {
			size_t j = i + 1;
			if (j == 3)
				j = 0;
			
			edge[i] = CGAL::Vector_2<K>(triangle.vertex(i), triangle.vertex(j));
		}

		std::array <int, 3> angle_type;

		for (size_t i = 0; i < 3; i++) {
			size_t j = i + 1;
			if (j == 3)
				j = 0;
			
			angle_type[i] = CGAL::angle(-edge[i], edge[j]);

			switch (angle_type[i]) {
				case CGAL::OBTUSE:
					std::cout << "Obtuse" << std::endl;
					std::cout << i << "[" << edge[i] << "]" << std::endl;
					std::cout << j << "[" << edge[j] << "]" << std::endl;
					ret++;
					break;
				default:
					break;
			}
		}
	}
	return ret;
}

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

		CDT cdt;
		data_t data{input};
		for (const auto& point : data.get_points())
			cdt.insert(point);
		for (const auto& edge : data.get_boundary())
			cdt.insert_constraint(edge.first, edge.second);
		for (const auto& constraint : data.get_constraints())
			cdt.insert_constraint(constraint.first, constraint.second);

		//CGAL::make_conforming_Delaunay_2(cdt);
		//CGAL::make_conforming_Gabriel_2(cdt);

		std::cout << count_obtuse(&cdt) << std::endl;
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
