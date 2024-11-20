#include <iterator>
#include <random>
#include <limits>

#include "configuration.hpp"
#include "config_cgal.hpp"

#include <CGAL/intersections.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Quotient.h>
#include <CGAL/Polygon_mesh_processing/refine.h>

#include "triangulation.hpp"
#include "obtuse_polygon.hpp"
#include "helper.hpp"

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>   Mesher;


std::array<std::set<std::tuple<CDT::Point, CDT::Point, CDT::Point>>, st_end> global_tried;

// Helper functions

static size_t count_obtuse(CDT *cdt, data_t *data)
{
	size_t ret = 0;
	for (auto it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {
		K::Triangle_2 triangle = cdt->triangle(it);
		if (!(data->inside(CGAL::centroid(triangle))))
			continue;
		if (is_obtuse(triangle))
			ret++;
	}
	return ret;
}

static size_t count_all_obtuse(CDT *cdt, data_t *data)
{
	size_t ret = 0;
	for (auto it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {
		K::Triangle_2 triangle = cdt->triangle(it);
		if (is_obtuse(triangle))
			ret++;
	}
	return ret;
}

static void print_st_method(int method)
{
	#if PRINT_METHODS
	std::string name[st_end];
	std::fill_n(name, st_end, "");
	name[st_centroid] = "centroid";
	name[st_circumcenter] = "circumcenter";
	name[st_polygon_centroid] = "polygon";
	name[st_projection_outward] = "projection_outward";
	name[st_projection_all] = "projection_all";
	name[st_neighbor_random] = "neighbor";
	name[st_constraint_random] = "constraint";

	if (method < st_end)
		std::cout << name[method] << std::endl;
	#endif
}

// Class functions

triangulation_t::triangulation_t(data_t *data)
{
	for (const auto& point : data->get_points())
		this->cdt.insert(point);
	for (const auto& edge : data->get_boundary())
		this->cdt.insert_constraint(edge.first, edge.second);
	for (const auto& constraint : data->get_constraints())
		this->cdt.insert_constraint(constraint.first, constraint.second);

	this->data = data;
	this->start_obtuse = count_obtuse(&(this->cdt), this->data);
	this->obtuse = this->start_obtuse;
	this->update_outside_obtuse();
	this->steiner = 0;
	this->progression_check = progression_less;

	this->tried = &global_tried;
}

void triangulation_t::update_outside_obtuse()
{
	this->outside_obtuse = 0;
	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		K::Triangle_2 triangle = this->cdt.triangle(it);
		if (
			!(data->inside(CGAL::centroid(triangle)))
			&& is_obtuse(triangle)
		)
			this->outside_obtuse++;	
	}
}

static inline bool never_on_boundary(int method)
{
	switch (method) {
		case st_centroid:
			return true;
		case st_constraint_random:
			return true;
		case st_neighbor_random:
			return true;
		case st_polygon_centroid:
			return true;
		default:
			return false;
	}
}

void triangulation_t::insert(CDT::Point steiner, int method)
{
	// If we are sure the point won't be on the boundary, use a faster method
	if (never_on_boundary(method)) {
		if (!this->data->inside(steiner))
			return;

		this->cdt.insert(steiner);
		this->obtuse = count_all_obtuse(&(this->cdt), this->data) - this->outside_obtuse;
	} else {
		if (!this->data->inside(steiner)
		&& !this->data->on_boundary(steiner))
			return;

		this->cdt.insert(steiner);
		if (this->data->on_boundary(steiner)) {
			//this->obtuse = count_obtuse(&(this->cdt), this->data);
			this->update_outside_obtuse();
			this->obtuse = count_all_obtuse(&(this->cdt), this->data) - this->outside_obtuse;

		} else {
			this->obtuse = count_all_obtuse(&(this->cdt), this->data) - this->outside_obtuse;
		}
	}

	this->steiner++;
	this->history.push_back(method);
}

static inline std::tuple<CDT::Point, CDT::Point, CDT::Point> triangle_to_tuple(K::Triangle_2 triangle)
{
	std::tuple<CDT::Point, CDT::Point, CDT::Point> ret(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));
	return ret;
}

void triangulation_t::set_progression_check(int check)
{
	this->progression_check = check;
}

size_t triangulation_t::size_obtuse()	
{
	return this->obtuse;
}

size_t triangulation_t::size_steiner()	
{
	return this->steiner;
}

CDT triangulation_t::get_cdt()
{
	return this->cdt;
}

/*
 * At least this or one of its neighbors are obtuse
 */

static inline bool detect_obtuse(
	CDT *cdt,
	CDT::Face_iterator *it,
	K::Triangle_2 triangle)
{
	if(is_obtuse(triangle))
		return true;

	for (size_t i = 0; i < 3; i++) {
		auto neighbor = (*it)->neighbor(i);
		if (cdt->is_infinite(neighbor))
			continue;

		K::Triangle_2 tmp = cdt->triangle(neighbor);
		if(is_obtuse(tmp))
			return true;
	}
	return false;
}

/*
 * At least this or one of its neighbors haven't been tried
 */

static inline bool detect_nottried(
	CDT *cdt,
	std::set<std::tuple<CDT::Point, CDT::Point, CDT::Point>> *tried,
	CDT::Face_iterator *it,
	K::Triangle_2 triangle)
{
	//return true;
	if (tried->find(triangle_to_tuple(triangle)) == tried->end())
		return true;

	for (size_t i = 0; i < 3; i++) {
		auto neighbor = (*it)->neighbor(i);
		if (cdt->is_infinite(neighbor))
			continue;

		K::Triangle_2 tmp = cdt->triangle(neighbor);
		if (tried->find(triangle_to_tuple(tmp)) == tried->end())
			return true;
	}
	return false;
}

// Steiner methods

void triangulation_t::steiner_centroid(
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions)
{
	const int method = st_centroid;

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);
		
		if (!detect_obtuse(&(this->cdt), &it, triangle))
			continue;
		if (!detect_nottried(&(this->cdt), &((*(this->tried))[method]), &it, triangle))
			continue;

		std::pair<K::Triangle_2, std::vector<CDT::Point>> solution;
		solution.first = triangle;
		solution.second.push_back(CGAL::centroid(triangle));
		solutions->push_back(solution);
	}
}

void triangulation_t::steiner_circumcenter(
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions)
{
	const int method = st_circumcenter;

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);
		if (!detect_obtuse(&(this->cdt), &it, triangle))
			continue;
		if (!detect_nottried(&(this->cdt), &((*(this->tried))[method]), &it, triangle))
			continue;

		std::pair<K::Triangle_2, std::vector<CDT::Point>> solution;
		solution.first = triangle;
		solution.second.push_back(CGAL::circumcenter(triangle));

		solutions->push_back(solution);
	}
}

void triangulation_t::steiner_midpoint(std::vector<CDT::Point> *steiner_pts)
{
	std::vector<std::pair<K::Segment_2, double>> longest;
	for (auto it = this->cdt.finite_edges_begin(); it != this->cdt.finite_edges_end(); it++) {
		CDT::Edge e = *it;

		auto v1 = e.first->vertex( (e.second+1)%3 );
		auto v2 = e.first->vertex( (e.second+2)%3 );
	
		CDT::Point p1 = v1->point();
		CDT::Point p2 = v2->point();

		K::Segment_2 seg(p1, p2);
		double length = CGAL::to_double(seg.squared_length());

		if (longest.size() == 0 || length > longest.back().second) {
			std::pair<K::Segment_2, double> tmp;
			tmp.first = seg;
			tmp.second = length;
			longest.push_back(tmp);
		}
	}
	if (longest.size() > 0)
		steiner_pts->push_back(CGAL::midpoint(longest.back().first));
}

void triangulation_t::steiner_polygon_centroid(std::vector<CDT::Point> *steiner_pts)
{
	std::set<std::tuple<CDT::Point, CDT::Point, CDT::Point>> visited;

	while (1) {
		std::vector<CDT::Point> current_steiner_pts;
		obtuse_polygon_t obtuse_polygon(this->data);
		// Initialize first element
		for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
			K::Triangle_2 triangle = this->cdt.triangle(it);
			if (std::find(visited.begin(), visited.end(), triangle_to_tuple(triangle)) != visited.end())
				continue;

			if (!is_obtuse(triangle))
				continue;

			std::vector<K::Triangle_2> obtuse_neighbors;
			for (size_t i = 0; i < 3; i++) {
				auto neighbor = it->neighbor(i);
				if (this->cdt.is_infinite(neighbor))
					continue;

				K::Triangle_2 tmp = this->cdt.triangle(neighbor);
				if(!is_obtuse(tmp))
					continue;

				if (!this->data->inside(CGAL::centroid(tmp)))
					continue;

				if (std::find(visited.begin(), visited.end(), triangle_to_tuple(tmp)) != visited.end())
					continue;

				obtuse_neighbors.push_back(tmp);
			}
			if (obtuse_neighbors.size() == 1) {
				obtuse_polygon.try_insert(triangle);
				visited.insert(triangle_to_tuple(triangle));
				break;
			}
		}
		if (obtuse_polygon.size() == 0)
			break;

		// Continue the sequence:
		bool exit_flag = false;
		while (!exit_flag) {
			auto visited_before = visited;
			for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
				assert(obtuse_polygon.size() > 0);
				K::Triangle_2 triangle = this->cdt.triangle(it);
				if (std::find(visited.begin(), visited.end(), triangle_to_tuple(triangle)) != visited.end())
					continue;

				if(!obtuse_polygon.try_insert(triangle))
					continue;

				CDT::Point steiner_pt = obtuse_polygon.get_steiner();
				current_steiner_pts.push_back(steiner_pt);
				visited.insert(triangle_to_tuple(triangle));

				std::vector<K::Triangle_2> obtuse_neighbors;
				for (size_t i = 0; i < 3; i++) {
					auto neighbor = it->neighbor(i);
					if (this->cdt.is_infinite(neighbor))
						continue;

					K::Triangle_2 tmp = this->cdt.triangle(neighbor);
					if(!is_obtuse(tmp))
						continue;

					if (!this->data->inside(CGAL::centroid(tmp)))
						continue;

					obtuse_neighbors.push_back(tmp);
				}

				if (obtuse_neighbors.size() == 1) {
					exit_flag = true;
				}
				if (obtuse_neighbors.size() == 2) {
				}
				if (obtuse_neighbors.size() == 3){
					exit_flag = true;
				}
			}
			if (visited_before == visited)
				break;
		}

		CDT::Point steiner_pt = obtuse_polygon.get_steiner();
		current_steiner_pts.push_back(steiner_pt);

		std::reverse(current_steiner_pts.begin(), current_steiner_pts.end());
		steiner_pts->insert(steiner_pts->end(), current_steiner_pts.begin(), current_steiner_pts.end());
	}
}

void triangulation_t::steiner_projection_internal(
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *inward,
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *outward,
	const int method)
{
	assert (
		method == st_projection
		|| method == st_projection_all
		|| method == st_projection_outward
		|| method == st_projection_outward_all
		|| method == st_unused
	);

	std::vector<K::Segment_2> boundary_segments = this->data->get_boundary_segments();

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);
		if (!detect_obtuse(&(this->cdt), &it, triangle))
			continue;
		if (!detect_nottried(&(this->cdt), &((*(this->tried))[method]), &it, triangle))
			continue;

		std::array <CGAL::Vector_2<K>, 3> edge;

		for (size_t i = 0; i < 3; i++) {
			size_t j = i + 1;
			if (j == 3)
				j = 0;
			
			edge[i] = CGAL::Vector_2<K>(triangle.vertex(i), triangle.vertex(j));
		}

		std::array <int, 3> angle_type;

		for (size_t i = 0; i < 3; i++) {
			size_t j = i + 1;
			size_t k = j + 1;
			if (j > 2)
				j -= 3;
			if (k > 2)
				k -= 3;
			
			angle_type[i] = CGAL::angle(-edge[i], edge[j]);

			if (angle_type[i] != CGAL::OBTUSE)
				continue;

			CGAL::Line_2<K> l(triangle.vertex(i), triangle.vertex(k));
			CDT::Point steiner = l.projection(triangle.vertex(j));
			
			// Sort to inward / outward
			K::Segment_2 closest_segment_src = boundary_segments[0];
			auto distance_src = CGAL::squared_distance(triangle.vertex(j), closest_segment_src);

			for (auto it = boundary_segments.begin(); it < boundary_segments.end(); it++) {
				auto tmp = CGAL::squared_distance(triangle.vertex(j), *it);
				if (tmp < distance_src) {
					closest_segment_src = *it;
					distance_src = tmp;
				}
			}
			
			K::Segment_2 closest_segment_dst = boundary_segments[0];
			auto distance_dst = CGAL::squared_distance(steiner, closest_segment_dst);
			for (auto it = boundary_segments.begin(); it < boundary_segments.end(); it++) {
				auto tmp = CGAL::squared_distance(steiner, *it);
				if (tmp < distance_dst) {
					closest_segment_dst = *it;
					distance_dst = tmp;
				}
			}

			if (distance_src > distance_dst || distance_dst == 0) {
				if (outward != NULL) {
					std::pair<K::Triangle_2, std::vector<CDT::Point>> solution;
					solution.first = triangle;
					solution.second.push_back(steiner);
					outward->push_back(solution);
				}
			}
			else {
				if (inward != NULL) {
					std::pair<K::Triangle_2, std::vector<CDT::Point>> solution;
					solution.first = triangle;
					solution.second.push_back(steiner);
					inward->push_back(solution);
				}
			}
			break;
		}
	}
}

void triangulation_t::steiner_projection_inward(
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions)
{
	const int method = st_unused;

	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> local;
	steiner_projection_internal(&local, NULL, method);

	if (solutions != NULL) {
		solutions->insert(solutions->end(), local.begin(), local.end());
		local.clear();
	}

	for (size_t i = 0; i < local.size(); i++) {
		for (size_t j = 0; j < local[i].second.size(); j++) {
			this->insert(local[i].second[j], method);
			(*(this->tried))[method].insert(triangle_to_tuple(local[i].first));

			print_st_method(method);
			if (this->obtuse == 0)
				return;
		}
	}
}

void triangulation_t::steiner_projection_outward(
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions)
{
	int method = st_projection_outward;
	if (solutions == NULL)
		method = st_projection_outward_all;

	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> local;
	steiner_projection_internal(NULL, &local, method);

	if (solutions != NULL) {
		solutions->insert(solutions->end(), local.begin(), local.end());
		local.clear();
	}

	for (size_t i = 0; i < local.size(); i++) {
		for (size_t j = 0; j < local[i].second.size(); j++) {
			this->insert(local[i].second[j], method);
			(*(this->tried))[method].insert(triangle_to_tuple(local[i].first));

			print_st_method(method);
			if (this->obtuse == 0)
				return;
		}
	}
}

void triangulation_t::steiner_projection(
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions)
{
	int method = st_projection;
	if (solutions == NULL)
		method = st_projection_all;

	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> local;

	steiner_projection_internal(&local, &local, method);

	if (solutions != NULL) {
		solutions->insert(solutions->end(), local.begin(), local.end());
		local.clear();
	}

	for (size_t i = 0; i < local.size(); i++) {
		for (size_t j = 0; j < local[i].second.size(); j++) {
			this->insert(local[i].second[j], method);
			(*(this->tried))[method].insert(triangle_to_tuple(local[i].first));

			print_st_method(method);
			if (this->obtuse == 0)
				return;
		}
	}
}

void triangulation_t::steiner_constraint_random(std::vector<CDT::Point> *steiner_pts)
{
	std::vector<std::pair<CDT::Point, CDT::Point>> constraints = this->data->get_constraints();

	for (auto it = constraints.begin(); it < constraints.end(); it++) {
		CGAL::Random_points_on_segment_2<CDT::Point> generator(it->first, it->second);

		K::Segment_2 seg(it->first, it->second);
		double squared_length = CGAL::to_double(seg.squared_length());
		size_t length = std::round(std::sqrt(squared_length));
		size_t amount = 100;

		if (length / 10 > 100)
			amount = length / 10;

		//std::cout << amount << std::endl;
		std::copy_n(generator, 100, std::back_inserter(*steiner_pts));
	}
}

void triangulation_t::steiner_neighbor_random(
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions)
{
	const int method = st_neighbor_random;
	std::vector<K::Triangle_2> triangles;

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);

		if (!detect_obtuse(&(this->cdt), &it, triangle))
			continue;
		if (!detect_nottried(&(this->cdt), &((*(this->tried))[method]), &it, triangle))
			continue;

		size_t area = std::round(CGAL::to_double(triangle.area()));
		if (area < 100)
			continue;

		for (size_t i = 0; i < 3; i++) {
			auto neighbor = it->neighbor(i);	
			if (this->cdt.is_infinite(neighbor))
				continue;

			auto triangle = this->cdt.triangle(neighbor);
			triangles.push_back(triangle);
		}
	}

	std::vector<K::Segment_2> constraint_segments = this->data->get_constraint_segments();

	for (size_t i = 0; i < triangles.size(); i++){
		std::vector<CDT::Point> random_pts;
		CDT cdt2;
		K::Triangle_2 triangle = triangles[i];

		bool touch_constraint = false;
		for (auto it = constraint_segments.begin(); it < constraint_segments.end(); it++) {
			const auto result = CGAL::intersection(*it, triangle);
			if (result) {
				touch_constraint = true;
				break;
			}
		}
		if (touch_constraint) {
			continue;
		}

		K::FT area_exact = triangle.area();

		size_t area = std::round(CGAL::to_double(area_exact));
		size_t amount = 100;
		if (area / 10 < amount)
			amount = area / 10;

		CGAL::Random_points_in_triangle_2<CDT::Point> generator(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));

		std::vector<CDT::Point> steiner_pts;
		std::copy_n(generator, amount, std::back_inserter(steiner_pts));

		std::pair<K::Triangle_2, std::vector<CDT::Point>> solution;
		solution.first = triangle;
		solution.second = steiner_pts;

		solutions->push_back(solution);
	}
}

void triangulation_t::steiner_add(const int method, size_t max)
{
	size_t amount = 0;
	std::random_device rd;
	std::mt19937 g(rd());
	std::vector<CDT::Point> steiner_pts;
	std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> solutions;
	switch (method){
		case st_centroid:
			this->steiner_centroid(&solutions);
			break;
		case st_circumcenter:
			this->steiner_circumcenter(&solutions);
			break;
		case st_constraint_random:
			this->steiner_constraint_random(&steiner_pts);
			break;
		case st_midpoint:
			this->steiner_midpoint(&steiner_pts);
			break;
		case st_neighbor_random:
			this->steiner_neighbor_random(&solutions);
			break;
		case st_polygon_centroid:
			this->steiner_polygon_centroid(&steiner_pts);
			break;
		case st_projection_outward:
			this->steiner_projection_outward(&solutions);
			break;
		case st_projection_outward_all:
			if (max == SIZE_MAX) {
				this->steiner_projection_outward(NULL);
			}
			else {
				this->steiner_projection_outward(&solutions);
			}
			break;
		case st_projection:
			this->steiner_projection(&solutions);
			break;
		case st_projection_all:
			if (max == SIZE_MAX) {
				this->steiner_projection(NULL);
			}
			else {
				this->steiner_projection(&solutions);
			}
			break;
		default:
			// maybe error?
			break;
	}
	//std::shuffle(steiner_pts.begin(), steiner_pts.end(), g);

	for (size_t i = 0; i < steiner_pts.size(); i++) {
		triangulation_t current = *this;
		current.insert(steiner_pts[i], method);

		if (this->progression_check == progression_less
			&& current.obtuse < this->obtuse) {
			*this = current;
			print_st_method(method);
			amount++;

			if (amount >= max)
				return;
		}

		if (this->progression_check == progression_less_equal
			&& current.obtuse <= this->obtuse
			&& method != st_constraint_random
			&& method != st_neighbor_random) {
			*this = current;
			print_st_method(method);
			amount++;

			if (amount >= max)
				return;
		}

		if (this->obtuse == 0)
			return;
	}


	for (size_t i = 0; i < solutions.size(); i++) {
		(*(this->tried))[method].insert(triangle_to_tuple(solutions[i].first));

		for (size_t j = 0; j < solutions[i].second.size(); j++) {
			triangulation_t current = *this;
			current.insert(solutions[i].second[j], method);


			if (this->progression_check == progression_less
				&& current.obtuse < this->obtuse) {
				*this = current;
				print_st_method(method);
				amount++;

				if (amount >= max)
					return;
			}

			if (this->progression_check == progression_less_equal
				&& current.obtuse <= this->obtuse
				&& method != st_constraint_random
				&& method != st_neighbor_random){
				*this = current;
				print_st_method(method);
				amount++;

				if (amount >= max)
					return;
			}

			if (this->obtuse == 0)
				return;
		}
	}
}
void triangulation_t::steiner_add(const int method)
{
	this->steiner_add(method, SIZE_MAX);
}

bool triangulation_t::exit_early()
{
	if (this->obtuse == 0)
		return true;
	
	if (this->steiner > this->start_obtuse * 3)
		return true;
	if (this->start_obtuse * 3 - this->steiner < this->obtuse)
		return true;
	
	return false;
}

std::vector<std::pair<std::string, std::string>> triangulation_t::get_steiner_str()
{
	std::vector<std::pair<std::string, std::string>> ret;
	std::vector<CDT::Point> input_points = this->data->get_points();
	std::vector<std::pair<CDT::Point, CDT::Point>> input_edges = this->data->get_constraints();

	for (auto it = this->cdt.finite_vertices_begin(); it != this->cdt.finite_vertices_end(); it++) {
		CDT::Point pt = it->point();
		if (std::find(input_points.begin(), input_points.end(), pt) != input_points.end())
			continue;
		
		std::stringstream steiner_x;
		std::stringstream steiner_y;

		const auto exact_x = CGAL::exact(pt.x());
		steiner_x << exact_x.get_num() << "/" << exact_x.get_den();

		const auto exact_y = CGAL::exact(pt.y());
		steiner_y << exact_y.get_num() << "/" << exact_y.get_den();

		std::pair<std::string, std::string> element;
		element.first = steiner_x.str();
		element.second = steiner_y.str();

		ret.push_back(element);
	}
	return ret;
}
	
std::vector<std::pair<size_t, size_t>> triangulation_t::get_edges()
{
	std::vector<std::pair<size_t, size_t>> ret;
	auto input_boundary = this->data->get_boundary();
	auto input_constraints = this->data->get_constraints();

	for (auto it = this->cdt.finite_edges_begin(); it != this->cdt.finite_edges_end(); it++) {
		CDT::Edge e = *it;

		auto v1 = e.first->vertex( (e.second+1)%3 );
		auto v2 = e.first->vertex( (e.second+2)%3 );
	
		CDT::Point p1 = v1->point();
		CDT::Point p2 = v2->point();
		
		CDT::Point center = CGAL::midpoint(p1, p2);
		if (!this->data->inside(center))
			continue;
		
		std::pair<size_t, size_t> id;
		size_t i = 0;
		for (auto it2 = cdt.finite_vertices_begin(); it2 != cdt.finite_vertices_end(); it2++, i++) {
			if (it2->point() == p1)
				id.first = i;
			if (it2->point() == p2)
				id.second = i;
		}
		ret.push_back(id);
	}
	return ret;
}
