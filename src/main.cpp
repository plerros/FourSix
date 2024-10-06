#include <boost/json.hpp>
#include <boost/json/src.hpp>

#include <iostream>
//#include <fstream>

#include "json.hpp"

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cerr <<
			"Usage: pretty <filename>"
			<< std::endl;
		return EXIT_FAILURE;
	}

	try {
		// Parse the file as JSON
		auto const jv = parse_file( argv[1] );

		// Now pretty-print the value
		pretty_print(std::cout, jv);
	}
	catch(std::exception const& e) {
		std::cerr <<
			"Caught exception: "
			<< e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
