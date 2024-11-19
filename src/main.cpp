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
	if (argc != 3 && argc != 5) {
		std::cerr << "Usage: FourSix -i <input json>" << std::endl;
		std::cerr << "       FourSix -i <input json> -o <output json>" << std::endl;
		return EXIT_FAILURE;
	}
	if (argc == 5 && OUTPUT_TRIANGULATION == false) {
		std::cerr << "Usage: FourSix -i <input json>" << std::endl;
		std::cout << std::endl;
		std::cout << "Unavailable with current configuration.h:" << std::endl;
		std::cerr << "       FourSix -i <input json> -o <output json> | (needs OUTPUT_TRIANGULATION == true)" << std::endl;
		return EXIT_FAILURE;
	}

	try {
		// Parse the file as JSON
		auto const jv = parse_file(argv[2]);

		// Now pretty-print the value
		//pretty_print(std::cout, jv);

		data_in input{jv};
		//input.print();

		data_t data{input};
		triangulation_t triangulation(&data);
		std::cout << "Input Obtuse: " << triangulation.size_obtuse() << std::endl;

		//CGAL::make_conforming_Delaunay_2(cdt);
		//CGAL::make_conforming_Gabriel_2(cdt);
		
		// Runtime start
		auto t1 = std::chrono::high_resolution_clock::now();

		//steiner_mixed(&triangulation, &data, 5);

		for (int i = 0; i < data.get_optim_methods().size(); i++) {
			switch (data.get_optim_methods()[i]) {
				case om_my: {
					unsigned int depth = triangulation.size_obtuse() / 3;
					if (triangulation.size_obtuse() > 0 && depth < 5)
						depth = 5;

					triangulation.optim_mixed_recursive(depth);
					triangulation.set_progression_check(progression_less_equal);
					triangulation.optim_mixed_recursive(depth);
					break;
				}

				case om_ls:
					triangulation.optim_local_search();
					break;

				case om_sa:
					triangulation.optim_simulated_annealing();
					break;

				case om_ant:
					triangulation.optim_ant_colony();
					break;

				default:
					break;
			}
		}
		
		// Runtime end
		auto t2 = std::chrono::high_resolution_clock::now();
		auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		std::chrono::duration<double, std::milli> ms_double = t2 - t1;
		std::cout << ms_int.count() << "ms" << std::endl;
		std::cout << ms_double.count() << "ms" << std::endl;

		std::cout << "obtuse " << triangulation.size_obtuse() << std::endl;
		std::cout << "steiner " << triangulation.size_steiner() << std::endl;

		if (OUTPUT_TRIANGULATION) {
			data_out output(&data, &triangulation);

			CGAL::draw(triangulation.get_cdt());

			if (argc == 5) {
				std::ofstream outfile;
				outfile.open(argv[4]);
				pretty_print(outfile, output.get_jsonvalue());
				outfile.close();
			} else {
				pretty_print(std::cout, output.get_jsonvalue());
			}
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
