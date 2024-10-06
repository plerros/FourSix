#pragma once
#include <boost/json.hpp>

boost::json::value parse_file(char const* filename);

void pretty_print (
	std::ostream& os,
	boost::json::value const& jv,
	std::string* indent = nullptr);
