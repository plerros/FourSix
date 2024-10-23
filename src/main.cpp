#include <boost/json.hpp>
#include <boost/json/src.hpp>
#include <chrono>
#include <iostream>

#include "configuration.hpp"

#include "json.hpp"
#include "data_in.hpp"
#include "data.hpp"
#include "data_out.hpp"
#include "triangulation.hpp"

int main(int argc, char** argv) {
	if (argc != 2 && argc != 3) {
		std::cerr << "Usage: FourSix <input json>" << std::endl;
		std::cerr << "       FourSix <input json> <output json>" << std::endl;
		return EXIT_FAILURE;
	}

	try {
		// Parse the file as JSON
		auto const jv = parse_file(argv[1]);

		// Now pretty-print the value
		//pretty_print(std::cout, jv);

		data_in input{jv};
		//input.print();

		data_t data{input};
		struct triangulation_t triangulation(&data);
		std::cout << "Input Obtuse: " << triangulation.get_obtuse() << std::endl;

		//CGAL::make_conforming_Delaunay_2(cdt);
		//CGAL::make_conforming_Gabriel_2(cdt);
		
		// Runtime start
		auto t1 = std::chrono::high_resolution_clock::now();


		//steiner_mixed(&triangulation, &data, 5);
		unsigned int depth = triangulation.get_obtuse() / 3;
		if (triangulation.get_obtuse() > 0 && depth < 5)
			depth = 5;

		triangulation.steiner_mixed_recursive(depth);
		triangulation.set_progression_check(progression_less_equal);
		triangulation.steiner_mixed_recursive(depth);


		// Runtime end
		auto t2 = std::chrono::high_resolution_clock::now();
		auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		std::chrono::duration<double, std::milli> ms_double = t2 - t1;
		std::cout << ms_int.count() << "ms" << std::endl;
		std::cout << ms_double.count() << "ms" << std::endl;


		std::cout << "obtuse " << triangulation.get_obtuse() << std::endl;

		data_out output(&data, &triangulation);

		CGAL::draw(triangulation.get_cdt());

		if (argc == 3) {
			std::ofstream outfile;
			outfile.open(argv[2]);
			pretty_print(outfile, output.get_jsonvalue());
			outfile.close();
		} else {
			pretty_print(std::cout, output.get_jsonvalue());
		}
	}
	catch(std::exception const& e) {
		std::cerr <<
			"Caught exception: "
			<< e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
