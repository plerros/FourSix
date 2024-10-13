#include <boost/json.hpp>
#include <iostream>

#include "configuration.hpp"

#include "data.hpp"
#include "data_in.hpp"

#include <CGAL/Polygon_2.h>

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
		if (std::find(this->points.begin(), this->points.end(), edge.first) != points.end())
			inside++;
		if (std::find(this->points.begin(), this->points.end(), edge.second) != points.end())
			inside++;

		if (inside == 2) {
			it++;
			continue;
		}

		this->constraints.erase(it);
		it = this->constraints.begin();
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

	CGAL::Polygon_2 pgn(traits);
	for (auto it = this->boundary_pts.begin(); it < this->boundary_pts.end(); it++)
		pgn.push_back(*it);

	return (! pgn.has_on_unbounded_side(pt));
}
