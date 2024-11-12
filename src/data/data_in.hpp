#pragma once

#include <boost/json.hpp>

class data_in
{
	private:
		std::string instance_uid;
		std::vector<std::pair<std::int64_t, std::int64_t>> points;
		std::vector<size_t> region_boundary;
		std::vector<std::pair<size_t, size_t>> constraints;

		/*
		 * delauney:
		 * True  == Start from delauney
		 * False == Start from my own triangulation
		 */
		bool delauney;

	public:
		data_in(boost::json::value const& jv);
		void print();
		std::string get_instance_uid();
		std::vector<std::pair<std::int64_t, std::int64_t>> get_points();
		std::vector<size_t>  get_boundary();
		std::vector<std::pair<size_t, size_t>> get_constraints();
};