#include <boost/json.hpp>
#include <iostream>

#include "data_in.hpp"

data_in::data_in(boost::json::value const& jv)
{
	this->instance_uid = jv.at("instance_uid").as_string();
	assert(jv.at("num_points").as_int64() >= 0);

	{
		auto tmp = jv.at("points_x").as_array();
		assert(tmp.size() == jv.at("num_points").as_int64());

		if (tmp.empty())
			goto skip_points_x;

		for (auto it = tmp.begin(); it != tmp.end(); it++)
			this->points.push_back({it->as_int64(), 0});
	}
skip_points_x:

	{
		auto tmp = jv.at("points_y").as_array();
		assert(tmp.size() == jv.at("num_points").as_int64());

		if (tmp.empty())
			goto skip_points_y;

		size_t i = 0;
		for (auto it = tmp.begin(); it != tmp.end() && i < this->points.size(); it++, i++)
			this->points[i].second = it->as_int64();
	}
skip_points_y:

	{
		auto tmp = jv.at("region_boundary").as_array();
		assert(tmp.size() <= jv.at("num_points").as_int64());

		if (tmp.empty())
			goto skip_region_boundary;

		for (auto it = tmp.begin(); it != tmp.end(); it++)
			this->region_boundary.push_back(it->as_int64());
	}
skip_region_boundary:

	assert(jv.at("num_constraints").as_int64() >= 0);

	{
		auto tmp = jv.at("additional_constraints").as_array();
		assert(tmp.size() == jv.at("num_constraints").as_int64());

		if (tmp.empty())
			goto skip_additional_constraints;

		for (auto it = tmp.begin(); it != tmp.end(); it++) {
			auto arr = it->as_array();
			assert(arr[0].as_int64() >= 0);
			assert(arr[1].as_int64() >= 0);

			size_t arr0 = arr[0].as_int64();
			size_t arr1 = arr[1].as_int64();

			this->constraints.push_back({arr0, arr1});
		}
	}
skip_additional_constraints:

	if (1) {}; // Needed for previous label
}

void data_in::print()
{
	std::cout << this->instance_uid << "\n";

	for (const auto& point : this->points)
		std::cout << "{" << point.first << ", " << point.second << "}\n";
	for (const auto& point : this->region_boundary)
		std::cout << point << "\n";
}

std::vector<std::pair<std::int64_t, std::int64_t>> data_in::get_points()
{
	return this->points;
}

std::vector<size_t>  data_in::get_boundary()
{
	return this->region_boundary;
}

std::vector<std::pair<size_t, size_t>> data_in::get_constraints()
{
	return this->constraints;
}
