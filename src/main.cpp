#include <boost/json.hpp>
#include <boost/json/src.hpp>
#include <chrono>
#include <iostream>

#include "configuration.hpp"
#include "config_cgal.hpp"

#include "json.hpp"
#include "data_in.hpp"
#include "data.hpp"
#include "data_out.hpp"
#include "triangulation.hpp"

int main(int argc, char** argv) {
	std::srand(std::time(nullptr));
	char **inname = NULL;
	char **outname = NULL;
	if (argc == 3 && strcmp(argv[1], "-i") == 0) {
		inname = &(argv[2]);
	}

	if (argc == 5 && strcmp(argv[1], "-i") == 0 && strcmp(argv[3], "-o") == 0) {
		inname = &(argv[2]);
		outname = &(argv[4]);
	}
	if (argc == 5 && strcmp(argv[1], "-o") == 0 && strcmp(argv[3], "-i") == 0) {
		outname = &(argv[2]);
		inname = &(argv[4]);
	}

	if (inname == NULL) {
		std::cerr << "Usage: FourSix -i <input json>" << std::endl;
		std::cerr << "       FourSix -i <input json> -o <output json>" << std::endl;
		return EXIT_FAILURE;
	}
	if (outname != NULL && OUTPUT_TRIANGULATION == false) {
		std::cerr << "Usage: FourSix -i <input json>" << std::endl;
		std::cout << std::endl;
		std::cout << "Unavailable with current configuration.h:" << std::endl;
		std::cerr << "       FourSix -i <input json> -o <output json> | (needs OUTPUT_TRIANGULATION == true)" << std::endl;
		return EXIT_FAILURE;
	}

	try {
		// Parse the file as JSON
		auto const jv = parse_file(*inname);

		data_in input{jv};

		data_t data{input};
		triangulation_t triangulation(&data);
		std::cout << "Input Obtuse: " << triangulation.size_obtuse() << std::endl;
		
		// Runtime start
		auto t1 = std::chrono::high_resolution_clock::now();

		auto parameters = data.get_alg();
		for (auto it =parameters.begin(); it != parameters.end(); it++) {
			switch (it->method) {
				case om_my: {
					it->L = triangulation.size_obtuse() / 3;
					if (triangulation.size_obtuse() > 0 && it->L < 5)
						it->L = 5;

					triangulation.optim_mixed_recursive(it->L);
					triangulation.set_progression_check(progression_less_equal);
					triangulation.optim_mixed_recursive(it->L);
					break;
				}

				case om_ls:
					triangulation.optim_local_search(*it);
					break;

				case om_sa:
					triangulation.optim_simulated_annealing(*it);
					break;

				case om_ant:
					triangulation.optim_ant_colony(*it);
					break;

				default:
					break;
			}
		}
		
		// Runtime end
		auto t2 = std::chrono::high_resolution_clock::now();
		auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		std::chrono::duration<double, std::milli> ms_double = t2 - t1;
		std::cout << "ms " << ms_int.count() << std::endl;

		std::cout << "obtuse " << triangulation.size_obtuse() << std::endl;
		std::cout << "steiner " << triangulation.size_steiner() << std::endl;

		if (OUTPUT_TRIANGULATION) {
			data_out output(&data, &triangulation);

			CGAL::draw(triangulation.get_cdt());

			if (outname != NULL) {
				std::ofstream outfile;
				outfile.open(*outname);
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
