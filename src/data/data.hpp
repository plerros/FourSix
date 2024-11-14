#pragma once

#include "configuration.hpp"
#include "data_in.hpp"

enum optim_method_enum{om_local, om_sa, om_ant};

class data_t
{
	private:
		// Initial Input
		std::string instance_uid;
		std::vector<CDT::Point> points;
		std::vector<CDT::Point> boundary_pts;
		std::vector<std::pair<CDT::Point, CDT::Point>> boundary;
		std::vector<K::Segment_2> boundary_segments;
		std::vector<std::pair<CDT::Point, CDT::Point>> constraints;
		std::vector<K::Segment_2> constraint_segments;
		std::vector<CDT::Point> constraint_mid_pts;

		CGAL::Polygon_2<K> boundary_pgn;

		std::vector<int> optim_methods;

		double parameter_a;
		double parameter_b;
		double parameter_xi;
		double parameter_psi;
		double parameter_lambda;
		unsigned int parameter_kappa;
		unsigned int parameter_L;

	public:
		data_t(data_in d);
		void print();
		std::string get_instance_uid();
		std::vector<CDT::Point> get_points();
		std::vector<std::pair<CDT::Point, CDT::Point>> get_boundary();
		std::vector<K::Segment_2> get_boundary_segments();
		std::vector<std::pair<CDT::Point, CDT::Point>> get_constraints();
		std::vector<K::Segment_2> get_constraint_segments();
		std::vector<CDT::Point> get_constraint_mid_pts();
		bool inside(CDT::Point pt);
		bool on_boundary(CDT::Point pt);

		void final_triangulation();
};