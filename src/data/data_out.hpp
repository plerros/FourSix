#pragma once

#include "data.hpp"
#include "triangulation.hpp"

class data_out
{
	private:
		std::string content_type;
		std::string instance_uid;
		std::vector<std::string> steiner_points_x;
		std::vector<std::string> steiner_points_y;
		std::vector<std::pair<size_t, size_t>> edges;

	public:
		data_out(data_t *data, triangulation_t *triangulation);
		boost::json::value get_jsonvalue();
};