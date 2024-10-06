#pragma once

#include <boost/json.hpp>

class indata_t
{
	private:
		std::string instance_uid;
		std::vector<std::vector<std::int64_t>> points;
		std::vector<size_t> region_boundary;
		std::vector<std::vector<std::size_t>> constraints;

	public:
		indata_t(boost::json::value const& jv);
		void print();
};