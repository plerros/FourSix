#include <boost/json.hpp>
#include <iostream>

#include "config_cgal.hpp"

#include <CGAL/Quotient.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_integer.h>

#include <CGAL/draw_polygon_2.h> 

#include "data.hpp"
#include "data_in.hpp"

enum graph_constraint_enum {
	gconst_unconstrained,
	gconst_circle,
	gconst_other,
	gconst_end
};

static bool boundary_repeats (std::vector<CDT::Point> boundary_pts)
{
	for (auto current = boundary_pts.begin(); current < boundary_pts.end(); current++) {
		auto next = current;
		next++;
		if (next == boundary_pts.end())
			next = boundary_pts.begin();
		
		if(*current == *next)
			return true;
	}
	return false;
}

static bool boundary_convex (std::vector<CDT::Point> boundary_pts)
{
	assert(boundary_repeats(boundary_pts) == false);
	CGAL::Polygon_2<K> polygon;
	for (auto it = boundary_pts.begin(); it < boundary_pts.end(); it++)
		polygon.push_back(*it);
	
	assert(polygon.is_simple());
	return (polygon.is_convex());
}

static bool boundary_orthogonal (
	std::vector<CDT::Point> boundary_pts,
	std::vector<std::pair<CDT::Point, CDT::Point>> boundary)
{
	assert(boundary_repeats(boundary_pts) == false);
	for (auto it = boundary.begin(); it < boundary.end(); it++) {
		if (it->first.x() == it->second.x())
			continue;
		if (it->first.y() == it->second.y())
			continue;
		return false;
	}
	return true;
}

static bool graph_cycle (std::vector<std::pair<CDT::Point, CDT::Point>> edges)
{
	std::set<CDT::Point> visited;
	std::vector<CDT::Point> current_boundary;
	auto not_visited = edges;

	while (not_visited.size() > 0) {
		if (current_boundary.size() == 0) {
			current_boundary.push_back(not_visited.back().first);
			visited.insert(current_boundary[0]);
		}
		std::vector<CDT::Point> next;

		for (auto it = current_boundary.begin(); it < current_boundary.end(); it++) {
			std::vector<std::pair<CDT::Point, CDT::Point>> tmp;
			for (auto it2 = not_visited.begin(); it2 < not_visited.end(); it2++) {				
				if (it2->first == *it) {
					next.push_back(it2->second);
				}
				else if (it2->second == *it) {
					next.push_back(it2->first);
				}
				else {
					tmp.push_back(*it2);
					continue;
				}

				if (visited.find(next.back()) != visited.end())
					return true;
				visited.insert(next.back());
			}
			not_visited = tmp;
		}
		current_boundary = next;
	}
	return false;
}

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

	// Categorize input geometry
	bool convex = boundary_convex(this->boundary_pts);
	bool orthogonal = boundary_orthogonal(this->boundary_pts, this->boundary);
	int constraints_category = gconst_other;

	if (this->constraints.size() == 0) {
		constraints_category = gconst_unconstrained;
	} else {
		std::vector<std::pair<CDT::Point, CDT::Point>> tmp = this->constraints;
		assert(this->boundary.size() > 2);
		//tmp.insert(tmp.end(), this->boundary.begin(), this->boundary.end() - 1);
		if (graph_cycle(tmp))
			constraints_category = gconst_circle;
	}

	this->category = gc_E_other;
	if (convex) {
		switch (constraints_category) {
			case gconst_unconstrained:
				this->category = gc_A_convex_simple;
				break;
			case gconst_other:
				this->category = gc_B_convex_line;
				break;
			case gconst_circle:
				this->category = gc_C_convex_circle;
				break;
		}
	}
	else if (orthogonal) {
		this->category = gc_D_ortho_simple;
	}

	int category_om = om_none;
	switch (this->category) {
		case gc_A_convex_simple:
			category_om = om_ant;
			break;
		case gc_B_convex_line:
			category_om = om_ant;
			break;
		case gc_C_convex_circle:
			category_om = om_ant;
			break;
		case gc_D_ortho_simple:
			category_om = om_ant;
			break;
		case gc_E_other:
			category_om = om_ant;
			break;
	}

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

	{
		auto tmp = default_val;
		if (d.get_optim_method() == "ls")
			tmp.method = om_ls;
		else if (d.get_optim_method() == "sa")
			tmp.method = om_sa;
		else if (d.get_optim_method() == "ant")
			tmp.method = om_ant;
		else if (d.get_optim_method() == "auto")
			tmp.method = category_om;

		switch (tmp.method) {
			case om_ant:
				tmp.xi     = d.get_parameter_xi();
				tmp.psi    = d.get_parameter_psi();
				tmp.lambda = d.get_parameter_lambda();
				tmp.kappa  = d.get_parameter_kappa();
				// Allow fall through
			case om_sa:
				tmp.a      = d.get_parameter_a();
				tmp.b      = d.get_parameter_b();
				// Allow fall through
			case om_ls:
				tmp.L      = d.get_parameter_L();
				this->alg.push_back(tmp);
			case om_none:
				break;
		}
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

int data_t::get_category()
{
	return this->category;
}
