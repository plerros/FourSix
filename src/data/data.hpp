#pragma once

#include "configuration.hpp"
#include "data_in.hpp"

class data_t
{
	private:
		std::vector<CDT::Point> points;
		std::vector<CDT::Point> boundary_pts;
		std::vector<std::pair<CDT::Point, CDT::Point>> boundary;
		std::vector<std::pair<CDT::Point, CDT::Point>> constraints;

	public:
		data_t(data_in d);
		void print();
		std::vector<CDT::Point> get_points();
		std::vector<std::pair<CDT::Point, CDT::Point>> get_boundary();
		std::vector<std::pair<CDT::Point, CDT::Point>> get_constraints();
		bool inside(CDT::Point pt);
};