#pragma once

#include "configuration.hpp"
#include "data_in.hpp"

class data_t
{
	private:
		// Initial Input
		std::string instance_uid;
		std::vector<CDT::Point> points;
		std::vector<CDT::Point> boundary_pts;
		std::vector<std::pair<CDT::Point, CDT::Point>> boundary;
		std::vector<std::pair<CDT::Point, CDT::Point>> constraints;
		std::vector<CDT::Point> constraint_mid_pts;

	public:
		data_t(data_in d);
		void print();
		std::string get_instance_uid();
		std::vector<CDT::Point> get_points();
		std::vector<std::pair<CDT::Point, CDT::Point>> get_boundary();
		std::vector<std::pair<CDT::Point, CDT::Point>> get_constraints();
		std::vector<CDT::Point> get_constraint_mid_pts();
		bool inside(CDT::Point pt);

		void final_triangulation();
};