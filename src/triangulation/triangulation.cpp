#include "configuration.hpp"

#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Quotient.h>
#include <CGAL/Polygon_mesh_processing/refine.h>

#include "triangulation.hpp"
#include <CGAL/draw_polygon_2.h> 
#include <CGAL/centroid.h>

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>   Mesher;

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

triangulation_t::triangulation_t(data_t *data)
{
	for (const auto& point : data->get_points())
		this->cdt.insert(point);
	for (const auto& edge : data->get_boundary())
		this->cdt.insert_constraint(edge.first, edge.second);
	for (const auto& constraint : data->get_constraints())
		this->cdt.insert_constraint(constraint.first, constraint.second);

	this->data = data;
	this->obtuse = count_obtuse(&(this->cdt), this->data);
	this->progression_check = progression_less;
}

void triangulation_t::set_progression_check(int check)
{
	this->progression_check = check;
}

size_t triangulation_t::get_obtuse()	
{
	return this->obtuse;
}

CDT triangulation_t::get_cdt()
{
	return this->cdt;
}

void triangulation_t::steiner_centroid()
{
	std::vector<CDT::Point> steiner_pts;

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);
		CDT::Point steiner = CGAL::centroid(triangle);
		if (!this->data->inside(steiner))
			continue;
		
		steiner_pts.push_back(steiner);
	}

	for (size_t i = 0; i < steiner_pts.size(); i++) {
		struct triangulation_t current = *this;
		current.cdt.insert(steiner_pts[i]);
		current.obtuse = count_obtuse(&(current.cdt), this->data);
		
		if (this->obtuse > current.obtuse)
			*this = current;
	}
}

void triangulation_t::steiner_circumcenter()
{
	std::vector<CDT::Point> steiner_pts;

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);
		CDT::Point steiner = CGAL::circumcenter(triangle);
		if (!this->data->inside(steiner))
			continue;
		
		steiner_pts.push_back(steiner);
	}

	for (size_t i = 0; i < steiner_pts.size(); i++) {
		struct triangulation_t current = *this;
		current.cdt.insert(steiner_pts[i]);
		current.obtuse = count_obtuse(&(current.cdt), this->data);

		if (this->progression_check == progression_less
			&& current.obtuse < this->obtuse)
			*this = current;

		if (this->progression_check == progression_less_equal
			&& current.obtuse <= this->obtuse)
			*this = current;

		if (current.obtuse == 0)
			return;
	}
}

static bool is_in_triangle(K::Triangle_2 triangle, CDT::Point pt)
{
	for (size_t i = 0; i < 3; i++)
		if (pt == triangle.vertex(i))
			return true;
	return false;
}

struct boundary_edge_t {
	CDT::Point a;
	CDT::Point b;
	bool has_neighbor;
	K::Triangle_2 neighbor;
};

void triangulation_t::steiner_polygon_centroid()
{
	std::vector<K::Triangle_2> visited;
	std::vector<CDT::Point> steiner_pts;

	std::vector<K::Triangle_2> obtuse_sequence;
	while (1) {
		obtuse_sequence.clear();
		// Initialize first element
		for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
			K::Triangle_2 triangle = this->cdt.triangle(it);
			if (std::find(visited.begin(), visited.end(), triangle) != visited.end())
				continue;

			if (!this->data->inside(CGAL::centroid(triangle)))
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

				if (std::find(visited.begin(), visited.end(), tmp) != visited.end())
					continue;

				obtuse_neighbors.push_back(tmp);
			}
			if (obtuse_neighbors.size() == 1) {
				obtuse_sequence.push_back(triangle);
				visited.push_back(triangle);
				break;
			}
		}
		if (obtuse_sequence.size() == 0)
			break;

		// Continue the sequence:
		bool exit_flag = false;
		while (!exit_flag) {
			std::vector<K::Triangle_2> visited_before = visited;
			for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
				assert(obtuse_sequence.size() > 0);
				K::Triangle_2 triangle = this->cdt.triangle(it);
				if (std::find(visited.begin(), visited.end(), triangle) != visited.end())
					continue;

				if (!this->data->inside(CGAL::centroid(triangle)))
					continue;

				if (!is_obtuse(triangle))
					continue;

				std::vector<K::Triangle_2> obtuse_neighbors;
				bool connected = false;
				for (size_t i = 0; i < 3; i++) {
					auto neighbor = it->neighbor(i);
					if (this->cdt.is_infinite(neighbor))
						continue;

					K::Triangle_2 tmp = this->cdt.triangle(neighbor);
					if(!is_obtuse(tmp))
						continue;

					if (!this->data->inside(CGAL::centroid(tmp)))
						continue;

					if (tmp == obtuse_sequence.back())
						connected = true;
					obtuse_neighbors.push_back(tmp);
				}
				if (!connected)
					continue;

				assert(obtuse_neighbors.size() != 0);
				if (obtuse_neighbors.size() == 1) {
					obtuse_sequence.push_back(triangle);
					visited.push_back(triangle);
					exit_flag = true;
					break;
				}
				if (obtuse_neighbors.size() == 2) {
					obtuse_sequence.push_back(triangle);
					visited.push_back(triangle);
					break;
				}
				if (obtuse_neighbors.size() == 3){
					obtuse_sequence.push_back(triangle);
					visited.push_back(triangle);
					exit_flag = true;
					break;
				}
			}
			assert(visited_before != visited);
		}

		assert(obtuse_sequence.size() > 1);
		// Generate the boundary
		std::vector<CDT::Point> pointsA;
		std::vector<CDT::Point> pointsB;

		K::Triangle_2 *prev = NULL;
		K::Triangle_2 prev2;
		while (obtuse_sequence.size() != 0) {
			K::Triangle_2 triangle = obtuse_sequence.back();
			obtuse_sequence.pop_back();
			size_t i = 0;
			size_t j = 1;
			size_t k = 2;

			assert(obtuse_sequence.size() != 0 || prev != NULL);
			for (; i < 3; i++) {
				j = i + 1;
				if (j == 3)
					j = 0;
				k = j + 1;
				if (k == 3)
					k = 0;
				if (
					(obtuse_sequence.size() > 0)
					&& is_in_triangle(obtuse_sequence.back(), triangle.vertex(i))
					&& is_in_triangle(obtuse_sequence.back(), triangle.vertex(j))
				)
					break;
				
				if (
					(obtuse_sequence.size() == 0)
					&& is_in_triangle(*prev, triangle.vertex(i))
					&& is_in_triangle(*prev, triangle.vertex(j))
				)
					break;
			}
			if (prev == NULL) {
				prev = &prev2;
				pointsA.push_back(triangle.vertex(k));
				pointsA.push_back(triangle.vertex(i));
				pointsB.push_back(triangle.vertex(j));
			} else {
				if (is_in_triangle(*prev, triangle.vertex(i)))
					pointsB.push_back(triangle.vertex(k));
				else
					pointsA.push_back(triangle.vertex(k));

				pointsA.push_back(triangle.vertex(i));
				pointsB.push_back(triangle.vertex(j));
			}
			*prev = triangle;

		}
		std::reverse(pointsB.begin(), pointsB.end());
		std::vector<CDT::Point> boundary;
		boundary.reserve(pointsA.size() + pointsB.size());
		boundary.insert(boundary.end(), pointsA.begin(), pointsA.end());
		boundary.insert(boundary.end(), pointsB.begin(), pointsB.end());

		/*
		K traits = K();
		CGAL::Polygon_2 pgn(traits);
		for (auto it = boundary.begin(); it < boundary.end(); it++)
			pgn.push_back(*it);
		
		CGAL::draw(pgn);
		*/

		CDT::Point steiner_pt = CGAL::centroid(boundary.begin(), boundary.end(), CGAL::Dimension_tag<0>());
		steiner_pts.push_back(steiner_pt);
	}
	for (size_t i = 0; i < steiner_pts.size(); i++) {
		struct triangulation_t current = *this;
		current.cdt.insert(steiner_pts[i]);
		current.obtuse = count_obtuse(&(current.cdt), this->data);
		
		if (this->obtuse > current.obtuse) {
			std::cout << "polycent" << std::endl;
			*this = current;
		}
	}	
}

void triangulation_t::steiner_projection()
{
	std::vector<CDT::Point> steiner_pts;

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);
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

			switch (angle_type[i]) {
				case CGAL::OBTUSE: {
					CGAL::Line_2<K> l(triangle.vertex(i), triangle.vertex(k));

					CDT::Point steiner = l.projection(triangle.vertex(j));
					
					if (this->data->inside(steiner))
						steiner_pts.push_back(steiner);
					break;
				}
				default:
					break;
			}
		}
	}
	for (size_t i = 0; i < steiner_pts.size(); i++)
		this->cdt.insert(steiner_pts[i]);

	this->obtuse = count_obtuse(&(this->cdt), this->data);
}

void triangulation_t::steiner_random()
{
	std::vector<CDT::Point> random_pts;
	{
		CDT cdt2 = this->cdt;
		Mesher mesher(cdt2);
		mesher.set_criteria(Criteria(0, 0));
		mesher.refine_mesh();
		CGAL::Random_points_in_triangle_mesh_2<CDT::Point, CDT> generator(cdt2);
		std::copy_n(generator, 100, std::back_inserter(random_pts));
		assert(random_pts.size() == 100);

	}

	for (size_t i = 0; i < random_pts.size(); i++) {
		struct triangulation_t current = *this;
		current.cdt.insert(random_pts[i]);
		current.obtuse = count_obtuse(&(current.cdt), this->data);
		if (current.obtuse < this->obtuse) {
			std::cout << "randpoint" << std::endl;
			*this = current;
		}
	}
}

void triangulation_t::steiner_neighbor_random()
{
	std::vector<std::array<CDT::Point, 3>> triangles;

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);
		std::array<CDT::Point, 3> pts;
		for (size_t i = 0; i < 3; i++)
			pts[i] = triangle.vertex(i);
	
		bool obtuse = false;
		for (size_t i = 0; i < 3; i++) {
			size_t j = i + 1;
			if (j > 2)
				j -= 3;
			size_t k = j + 1;
			if (k > 2)
				k -= 3;
			
			if (CGAL::angle(pts[i], pts[j], pts[k]) == CGAL::OBTUSE) {
				obtuse = true;
				break;
			}
		}

		if (obtuse)
			continue;
		
		//triangles.push_back(pts);

		for (size_t i = 0; i < 3; i++) {
			auto neighbor = it->neighbor(i);	
			if (this->cdt.is_infinite(neighbor))
				continue;

			auto triangle = this->cdt.triangle(neighbor);

			for (size_t i = 0; i < 3; i++)
				pts[i] = triangle.vertex(i);
			triangles.push_back(pts);
		}
	}

	for (size_t i = 0; i < triangles.size(); i++){
		std::vector<CDT::Point> random_pts;
		CDT cdt2;
		std::array<CDT::Point, 3> triangle = triangles[i];

		CGAL::Random_points_in_triangle_2<CDT::Point> generator(triangle[0], triangle[1], triangle[2]);
		std::copy_n(generator, 100, std::back_inserter(random_pts));


		for (size_t i = 0; i < random_pts.size(); i++) {
			struct triangulation_t current = *this;
			current.cdt.insert(random_pts[i]);
			current.obtuse = count_obtuse(&(current.cdt), this->data);
			if (current.obtuse < this->obtuse) {
				std::cout << "randpoint" << std::endl;
				*this = current;
			}
		}
	}
}

void triangulation_t::steiner_mixed(unsigned int retries)
{
	struct triangulation_t original = *this;
	for (unsigned int i = 0; i < retries; i++) {
		struct triangulation_t current = original;
		current.steiner_random();
		for (int j = 0; j < 20; j++) {
			struct triangulation_t current2 = current;
			current2.steiner_circumcenter();
			current2.steiner_projection();
			if (current.obtuse > current2.obtuse)
				current = current2;
		}
		if (this->obtuse > current.obtuse)
			*this = current;
	}
}

enum steiner_method_enum {st_circumcenter, st_projection, st_polygon_centroid, st_centroid, st_random};

void cout_space(size_t n)
{
	for (unsigned int i = 0; i < n; i++)
		std::cout << " ";
}

void triangulation_t::steiner_mixed_recursive(unsigned int depth)
{
	if (PRINT_RECURSION_TREE) {
		cout_space(depth);
		std::cout << depth << std::endl;
	}
	if (depth == 0)
		return;
	if (this->obtuse == 0)
		return;

	struct triangulation_t best = *this;
	struct triangulation_t current = *this;


	for (int method = st_circumcenter; method <= st_random; method ++) {
		current = *this;
		switch(method) {
			case st_circumcenter:
				current.steiner_circumcenter();
				break;
			case st_projection:
				current.steiner_projection();
				break;
			case st_polygon_centroid:
				if (this->progression_check == progression_less
					&& (best.obtuse != this->obtuse))
					break;
				current.steiner_polygon_centroid();
				break;
			case st_centroid:
				if (this->progression_check == progression_less
					&& (best.obtuse != this->obtuse))
					break;
				current.steiner_centroid();
				break;
			case st_random:
				if (this->progression_check == progression_less
					&& (best.obtuse != this->obtuse))
					break;
				current.steiner_neighbor_random();
				break;
		}

		if (this->progression_check == progression_less
			&& (current.obtuse < best.obtuse))
			best = current;
		if (this->progression_check == progression_less_equal
			&& (current.obtuse <= best.obtuse))
			best = current;

		if (best.obtuse == 0)
			break;

		if (current.obtuse != this->obtuse) {
			current.steiner_mixed_recursive(depth - 1);
			if (this->progression_check == progression_less
				&& (current.obtuse < best.obtuse))
				best = current;
			if (this->progression_check == progression_less_equal
				&& (current.obtuse <= best.obtuse))
				best = current;
			if (best.obtuse == 0)
				break;
		}
	}
	*this = best;
}

void triangulation_t::steiner_benchmark(unsigned int depth)
{
	for (unsigned int i = 0; i < depth; i++)
		this->steiner_centroid();
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
		std::pair<CDT::Point, CDT::Point> points1 (p1, p2);
		std::pair<CDT::Point, CDT::Point> points2 (p2, p1);

		if (std::find(input_boundary.begin(), input_boundary.end(), points1) != input_boundary.end())
			continue;
		if (std::find(input_boundary.begin(), input_boundary.end(), points2) != input_boundary.end())
			continue;

		if (std::find(input_constraints.begin(), input_constraints.end(), points1) != input_constraints.end())
			continue;
		if (std::find(input_constraints.begin(), input_constraints.end(), points2) != input_constraints.end())
			continue;

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