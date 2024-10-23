#include "configuration.hpp"

#include <CGAL/centroid.h>
#include <CGAL/draw_polygon_2.h> 
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Quotient.h>
#include <CGAL/Polygon_mesh_processing/refine.h>

#include "obtuse_polygon.hpp"

static bool is_in_triangle(K::Triangle_2 triangle, CDT::Point pt)
{
	for (size_t i = 0; i < 3; i++)
		if (pt == triangle.vertex(i))
			return true;
	return false;
}

static bool is_obtuse(K::Triangle_2 triangle)
{
	std::array <CGAL::Vector_2<K>, 3> edge;

	for (size_t i = 0; i < 3; i++) {
		size_t j = i + 1;
		if (j == 3)
			j = 0;
		
		edge[i] = CGAL::Vector_2<K>(triangle.vertex(i), triangle.vertex(j));
	}

	for (size_t i = 0; i < 3; i++) {
		size_t j = i + 1;
		if (j == 3)
			j = 0;

		switch (CGAL::angle(-edge[i], edge[j])) {
			case CGAL::OBTUSE:
				return true;
			default:
				break;
		}
	}
	return false;
}

obtuse_polygon_t::obtuse_polygon_t(data_t *data)
{
	this->data = data;
	this->boundary_valid = false;

	//K traits = K();
	//this->polygon(traits);
	this->polygon_valid = false;
};

void obtuse_polygon_t::recompute_boundary()
{
	this->boundary.clear();


	if (this->sequence.size() == 1) {
		K::Triangle_2 triangle = this->sequence[0];
		for (size_t i = 0; i < 3; i++)
			this->boundary.push_back(triangle.vertex(i));
		this->boundary_valid = true;
		return;
	}

	std::vector<CDT::Point> CounterClockWise;
	std::vector<CDT::Point> ClockWise;

	K::Triangle_2 *prev = NULL;
	K::Triangle_2 prev2;
	while (this->sequence.size() != 0) {
		K::Triangle_2 triangle = this->sequence.back();
		this->sequence.pop_back();

		if (prev == NULL) {
			// first loop
			std::pair<size_t, size_t> to_next(0, 1);

			for (;to_next.first < 3; to_next.first++, to_next.second = (to_next.first + 1) % 3)  {
				if (!is_in_triangle(this->sequence.back(), triangle.vertex(to_next.first)))
					continue;
				if (!is_in_triangle(this->sequence.back(), triangle.vertex(to_next.second)))
					continue;
				break;
			}

			size_t not_in_next = (to_next.second + 1) % 3;

			CounterClockWise.push_back(triangle.vertex(not_in_next));
			CounterClockWise.push_back(triangle.vertex(to_next.first));
			ClockWise.push_back(triangle.vertex(to_next.second));
			prev = &prev2;

		} else if (this->sequence.size() > 0) {
			// center of the sequence
			std::pair<size_t, size_t> to_next(0,1);
			std::pair<size_t, size_t> to_prev(0,1);

			for (;to_next.first < 3; to_next.first++, to_next.second = (to_next.first + 1) % 3)  {
				if (!is_in_triangle(this->sequence.back(), triangle.vertex(to_next.first)))
					continue;
				if (!is_in_triangle(this->sequence.back(), triangle.vertex(to_next.second)))
					continue;
				break;
			}

			for (;to_prev.first < 3; to_prev.first++, to_prev.second = (to_prev.first + 1) % 3)  {
				if (!is_in_triangle(*prev, triangle.vertex(to_prev.first)))
					continue;
				if (!is_in_triangle(*prev, triangle.vertex(to_prev.second)))
					continue;
				break;
			}

			size_t not_in_next = (to_next.second + 1) % 3;
			size_t not_in_prev = (to_prev.second + 1) % 3;

			if (to_prev.second == to_next.first) {
				/*
				 *    Next
				 *  *-------*
				 *   \     /
				 *    \   / Prev
				 *     \ /
				 *      * 
				 */
				ClockWise.push_back(triangle.vertex(not_in_prev));
			} else {
				/*
				 *         Next
				 *      *-------*
				 *       \     /
				 *   Prev \   /
				 *         \ /
				 *          * 
				 */
				CounterClockWise.push_back(triangle.vertex(not_in_prev));
			}
		} else {
			// end of the sequence
			std::pair<size_t, size_t> to_prev(0,1);

			for (;to_prev.first < 3; to_prev.first++, to_prev.second = (to_prev.first + 1) % 3)  {
				if (!is_in_triangle(*prev, triangle.vertex(to_prev.first)))
					continue;
				if (!is_in_triangle(*prev, triangle.vertex(to_prev.second)))
					continue;
				break;
			}
			size_t not_in_prev = (to_prev.second + 1) % 3;

			// Doesn't matter. Could go to ClockWise too
			CounterClockWise.push_back(triangle.vertex(not_in_prev));
		}
		*prev = triangle;
	}
	std::reverse(ClockWise.begin(), ClockWise.end());
	this->boundary.reserve(CounterClockWise.size() + ClockWise.size());
	this->boundary.insert(this->boundary.end(), CounterClockWise.begin(), CounterClockWise.end());
	this->boundary.insert(this->boundary.end(), ClockWise.begin(), ClockWise.end());
	this->boundary_valid = true;
}

void obtuse_polygon_t::recompute_polygon()
{
	K traits = K();
	this->polygon = CGAL::Polygon_2(traits);
	for (auto current = this->boundary.begin(); current < this->boundary.end(); current++) {
		auto next = current;
		next++;
		if (next == this->boundary.end())
			break;
		
		assert(*current != *next);
	}
	for (auto it = this->boundary.begin(); it < this->boundary.end(); it++)
		this->polygon.push_back(*it);
	
	assert(this->polygon.is_simple());
	this->polygon_valid = true;
}

void obtuse_polygon_t::draw()
{
	if (! this->boundary_valid)
		this->recompute_boundary();

	if (! this->polygon_valid)
		this->recompute_polygon();
	
	CGAL::draw(this->polygon);
}

CDT::Point obtuse_polygon_t::get_steiner()
{
	if (! this->boundary_valid)
		this->recompute_boundary();
	
	CDT::Point pt = CGAL::centroid(this->boundary.begin(), this->boundary.end(), CGAL::Dimension_tag<0>());
	return pt;
}

bool obtuse_polygon_t::internal(CDT::Point pt)
{
	if (! this->boundary_valid)
		this->recompute_boundary();

	if (! this->polygon_valid)
		this->recompute_polygon();

	if (this->polygon.has_on_unbounded_side(pt))
		return false;

	if (this->polygon.has_on_boundary(pt))
		return false;
	
	return true;
}

bool obtuse_polygon_t::has_constraint()
{
	std::vector<std::pair<CDT::Point, CDT::Point>> constraints = this->data->get_constraints();

	for (auto it = constraints.begin(); it < constraints.end(); it++) {
		if (this->internal(it->first))
			return true;
		if (this->internal(it->second))
			return true;
	}
	return false;
}

bool obtuse_polygon_t::is_convex()
{
	if (! this->boundary_valid)
		this->recompute_boundary();

	if (! this->polygon_valid)
		this->recompute_polygon();
	
	return this->polygon.is_convex();
}

size_t obtuse_polygon_t::size()
{
	return this->sequence.size();
}

bool obtuse_polygon_t::try_insert(K::Triangle_2 triangle)
{
	// Check if it's inside the boundary
	if (!this->data->inside(CGAL::centroid(triangle)))
		return false;

	if (!is_obtuse(triangle))
		return false;

	// Check if it's connected to previous
	if (this->sequence.size() != 0) {
		int common_edges = 0;
		for (size_t i = 0; i < 3; i++) {
			if (is_in_triangle(this->sequence.back(), triangle.vertex(i)))
				common_edges++;
		}
		if (common_edges != 2)
			return false;	
	}

	{
		struct obtuse_polygon_t tmp = *this;
		tmp.sequence.push_back(triangle);
		tmp.boundary_valid = false;
		tmp.polygon_valid = false;
		if(!tmp.is_convex())
			return false;

		if(tmp.has_constraint())
			return false;
	}

	this->sequence.push_back(triangle);
	this->boundary_valid = false;
	this->polygon_valid = false;

	return true;
}