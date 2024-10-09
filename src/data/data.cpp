#include <boost/json.hpp>
#include <iostream>

#include "configuration.hpp"

#include "data.hpp"
#include "data_in.hpp"

data_t::data_t(data_in d)
{
	// Convert std::pair<u64 u64> to CDT::Point
	for (const auto& point : d.get_points())
		this->points.push_back(Point(point.first, point.second));

	auto const boundary = d.get_boundary();

	for (size_t i = 0; i < boundary.size(); i++) {
		std::pair<CDT::Point, CDT::Point> edge;
		edge.first = this->points[boundary[i]];
		edge.second = this->points[boundary[0]];

		if (i + 1 < boundary.size())
			edge.second = this->points[boundary[i+1]];

		this->boundary.push_back(edge);
		this->boundary_pts.push_back(edge.first);
	}

	for (const auto& edge_id : d.get_constraints()) {
		auto point0_id = edge_id.first;
		auto point1_id = edge_id.second;

		std::pair<CDT::Point, CDT::Point> edge;

		edge.first = this->points[point0_id];
		edge.second = this->points[point1_id];
		
		this->constraints.push_back(edge);
	}
}

void data_t::print()
{
}

std::vector<CDT::Point> data_t::get_points()
{
	return this->points;
}

std::vector<std::pair<CDT::Point, CDT::Point>> data_t::get_boundary()
{
	return this->boundary;
}

std::vector<std::pair<CDT::Point, CDT::Point>> data_t::get_constraints()
{
	return this->constraints;
}

bool data_t::inside(CDT::Point pt)
{
	K traits = K();
	CDT::Point *front = &(this->boundary_pts.front());
	CDT::Point *back = &(this->boundary_pts.back());

	switch(CGAL::bounded_side_2(front, back, pt, traits)) {
		case CGAL::ON_BOUNDED_SIDE :
			return true;
		case CGAL::ON_BOUNDARY:
			return true;
		case CGAL::ON_UNBOUNDED_SIDE:
			return false;
	}
	return false;
}
