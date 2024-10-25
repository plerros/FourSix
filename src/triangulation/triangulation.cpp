#include "configuration.hpp"

#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Quotient.h>
#include <CGAL/Polygon_mesh_processing/refine.h>

#include "triangulation.hpp"
#include "obtuse_polygon.hpp"
#include "helper.hpp"

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>   Mesher;

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
		
		if (this->progression_check == progression_less
			&& current.obtuse < this->obtuse)
			*this = current;

		if (this->progression_check == progression_less_equal
			&& current.obtuse <= this->obtuse)
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

	while (1) {
		obtuse_polygon_t obtuse_polygon(this->data);
		// Initialize first element
		for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
			K::Triangle_2 triangle = this->cdt.triangle(it);
			if (std::find(visited.begin(), visited.end(), triangle) != visited.end())
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
				obtuse_polygon.try_insert(triangle);
				visited.push_back(triangle);
				break;
			}
		}
		if (obtuse_polygon.size() == 0)
			break;

		// Continue the sequence:
		bool exit_flag = false;
		while (!exit_flag) {
			std::vector<K::Triangle_2> visited_before = visited;
			for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
				assert(obtuse_polygon.size() > 0);
				K::Triangle_2 triangle = this->cdt.triangle(it);
				if (std::find(visited.begin(), visited.end(), triangle) != visited.end())
					continue;

				if(!obtuse_polygon.try_insert(triangle))
					continue;

				visited.push_back(triangle);

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
		steiner_pts.push_back(steiner_pt);
	}
	for (size_t i = 0; i < steiner_pts.size(); i++) {
		struct triangulation_t current = *this;
		current.cdt.insert(steiner_pts[i]);
		current.obtuse = count_obtuse(&(current.cdt), this->data);
		

		if (this->progression_check == progression_less
			&& current.obtuse < this->obtuse) {
			std::cout << "polycent" << std::endl;
			*this = current;
		}

		if (this->progression_check == progression_less_equal
			&& current.obtuse <= this->obtuse){
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
	std::vector<K::Triangle_2> triangles;

	for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
		auto triangle = this->cdt.triangle(it);

		size_t area = std::round(CGAL::to_double(triangle.area()));
		if (area < 100)
			continue;

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
			triangles.push_back(triangle);
		}
	}

	for (size_t i = 0; i < triangles.size(); i++){
		std::vector<CDT::Point> random_pts;
		CDT cdt2;
		K::Triangle_2 triangle = triangles[i];

		K::FT area_exact = triangle.area();

		size_t area = std::round(CGAL::to_double(area_exact));
		size_t amount = 100;
		if (area / 10 < amount)
			amount = area / 10;

		CGAL::Random_points_in_triangle_2<CDT::Point> generator(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));
		std::copy_n(generator, amount, std::back_inserter(random_pts));

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
		if (this->progression_check == progression_less_equal)
			std::cout << ".";
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

		if (current.cdt != this->cdt) {
			if (current.obtuse < this->obtuse)
				current.steiner_mixed_recursive(depth - 1);
			if (current.obtuse == this->obtuse)
				current.steiner_mixed_recursive(depth);
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
		this->steiner_polygon_centroid();
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