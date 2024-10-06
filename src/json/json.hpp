#ifndef UOA_ALGORITHM_PROJECT_JSON
#define UOA_ALGORITHM_PROJECT_JSON

#include <boost/json.hpp>

boost::json::value parse_file(char const* filename);

void pretty_print (
	std::ostream& os,
	boost::json::value const& jv,
	std::string* indent = nullptr);

#endif /* UOA_ALGORITHM_PROJECT_JSON */