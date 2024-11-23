#include <boost/json.hpp>
#include <iostream>

#include "config_cgal.hpp"

#include <CGAL/Quotient.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_integer.h>

#include "data.hpp"
#include "data_in.hpp"

data_t::data_t(data_in d)
{
	this->instance_uid = d.get_instance_uid();
	// Convert std::pair<u64 u64> to CDT::Point
	for (const auto& point : d.get_points()) {
		K::FT x(CGAL::Exact_integer(point.first));
		K::FT y(CGAL::Exact_integer(point.second));

		this->points.push_back(CDT::Point(x, y));
	}

	auto const boundary = d.get_boundary();
	K traits = K();
	this->boundary_pgn=CGAL::Polygon_2(traits);

	for (size_t i = 0; i < boundary.size(); i++) {
		std::pair<CDT::Point, CDT::Point> edge;
		edge.first = this->points[boundary[i]];
		edge.second = this->points[boundary[0]];

		if (i + 1 < boundary.size())
			edge.second = this->points[boundary[i+1]];

		this->boundary.push_back(edge);
		this->boundary_pts.push_back(edge.first);
		this->boundary_segments.push_back(K::Segment_2(edge.first, edge.second));
		this->boundary_pgn.push_back(edge.first);
	}

	for (const auto& edge_id : d.get_constraints()) {
		auto point0_id = edge_id.first;
		auto point1_id = edge_id.second;

		std::pair<CDT::Point, CDT::Point> edge;

		edge.first = this->points[point0_id];
		edge.second = this->points[point1_id];
		
		this->constraints.push_back(edge);
		this->constraint_segments.push_back(K::Segment_2(edge.first, edge.second));
	}

	// Remove points outside boundary
	for (auto it = this->points.begin(); it < this->points.end();) {
		CDT::Point point = *it;
		if (this->inside(point)) {
			it++;
			continue;
		}
		
		this->points.erase(it);
		it = this->points.begin();
	}

	for (auto it = this->constraints.begin(); it < this->constraints.end();) {
		std::pair<CDT::Point, CDT::Point> edge = *it;
		int inside = 0;
		if (std::find(this->points.begin(), this->points.end(), edge.first) != this->points.end())
			inside++;
		if (std::find(this->points.begin(), this->points.end(), edge.second) != this->points.end())
			inside++;

		if (inside == 2) {
			it++;
			continue;
		}

		this->constraints.erase(it);
		it = this->constraints.begin();
	}

	// Create constraint midpoints
	for (auto it = this->constraints.begin(); it < this->constraints.end(); it++)
		this->constraint_mid_pts.push_back(CGAL::midpoint(it->first, it->second));

	optim_alg_t default_val;
	default_val.method = om_none;
	default_val.a = 1.0;
	default_val.b = 1.0;
	default_val.psi = 0.0;
	default_val.lambda = 0.0;
	default_val.kappa = 0;
	default_val.L = 0;

	if (!d.get_delaunay()) {
		auto tmp = default_val;
		tmp.method = om_my;
		this->alg.push_back(tmp);
	}

	if (d.get_optim_method() == "ls") {
		auto tmp = default_val;
		tmp.method = om_ls;
		tmp.L      = d.get_parameter_L();
		this->alg.push_back(tmp);
	}
	if (d.get_optim_method() == "sa") {
		auto tmp = default_val;
		tmp.method = om_sa;
		tmp.a      = d.get_parameter_a();
		tmp.b      = d.get_parameter_b();
		tmp.L      = d.get_parameter_L();
		this->alg.push_back(tmp);
	}
	if (d.get_optim_method() == "ant") {
		auto tmp = default_val;
		tmp.method = om_ant;
		tmp.a      = d.get_parameter_a();
		tmp.b      = d.get_parameter_b();
		tmp.xi     = d.get_parameter_xi();
		tmp.psi    = d.get_parameter_psi();
		tmp.lambda = d.get_parameter_lambda();
		tmp.kappa  = d.get_parameter_kappa();
		tmp.L      = d.get_parameter_L();
		this->alg.push_back(tmp);
	}
}

data_t::data_t(K::Triangle_2 triangle, optim_alg_t alg)
{
	this->instance_uid = "dummy_uid";

	for (int i = 0; i < 3; i++) {
		this->boundary_pts.push_back(triangle.vertex(i));

		int j = (i+1) % 3;
		std::pair<CDT::Point, CDT::Point> edge;
		edge.first = triangle.vertex(i);
		edge.second = triangle.vertex(j);
		this->boundary.push_back(edge);

		this->boundary_segments.push_back(K::Segment_2(edge.first, edge.second));
		this->boundary_pgn.push_back(edge.first);
	}

	this->alg.push_back(alg);
}

void data_t::print()
{
}

std::string data_t::get_instance_uid()
{
	return this->instance_uid;
}

std::vector<CDT::Point> data_t::get_points()
{
	return this->points;
}

std::vector<std::pair<CDT::Point, CDT::Point>> data_t::get_boundary()
{
	return this->boundary;
}

std::vector<K::Segment_2> data_t::get_boundary_segments()
{
	return this->boundary_segments;
}

std::vector<std::pair<CDT::Point, CDT::Point>> data_t::get_constraints()
{
	return this->constraints;
}

std::vector<K::Segment_2> data_t::get_constraint_segments()
{
	return this->constraint_segments;
}

std::vector<CDT::Point> data_t::get_constraint_mid_pts()
{
	return this->constraint_mid_pts;
}

bool data_t::inside(CDT::Point pt)
{
	return (! this->boundary_pgn.has_on_unbounded_side(pt));
}

bool data_t::on_boundary(CDT::Point pt)
{
	return (this->boundary_pgn.has_on_boundary(pt));
}

std::vector<optim_alg_t> data_t::get_alg()
{
	return this->alg;
}
