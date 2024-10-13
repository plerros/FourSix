#include <boost/json.hpp>
#include <boost/json/src.hpp>
#include <chrono>
#include <iostream>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Line_2.h>

#include "configuration.hpp"

#include "json.hpp"
#include "data_in.hpp"
#include "data.hpp"

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>   Mesher;

size_t count_obtuse(CDT *cdt)
{
	size_t ret = 0;
	for (auto it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {
		auto triangle = cdt->triangle(it);
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
			if (j == 3)
				j = 0;
			
			angle_type[i] = CGAL::angle(-edge[i], edge[j]);

			switch (angle_type[i]) {
				case CGAL::OBTUSE:
					ret++;
					break;
				default:
					break;
			}
		}
	}
	return ret;
}

struct triangulation_t
{
	CDT cdt;
	size_t obtuse;
};

void steiner_random(struct triangulation_t *ptr)
{
	std::vector<Point> random_pts;
	{
		CDT cdt2 = ptr->cdt;
		Mesher mesher(cdt2);
		mesher.refine_mesh();
		CGAL::Random_points_in_triangle_mesh_2<Point, CDT> generator(cdt2);
		std::copy_n(generator, 100, std::back_inserter(random_pts));
		assert(random_pts.size() == 100);
	}

	for (size_t i = 0; i < random_pts.size(); i++) {
		struct triangulation_t current = *ptr;
		current.cdt.insert(random_pts[i]);
		current.obtuse = count_obtuse(&(current.cdt));
		if (current.obtuse < ptr->obtuse)
			*ptr = current;
	}
}

void steiner_projection(struct triangulation_t *ptr)
{
	std::vector<Point> steiner_pts;

	for (auto it = ptr->cdt.finite_faces_begin(); it != ptr->cdt.finite_faces_end(); it++) {
		auto triangle = ptr->cdt.triangle(it);
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

					Point steiner = l.projection(triangle.vertex(j));
					steiner_pts.push_back(steiner);
					break;
				}
				default:
					break;
			}
		}
	}
	for (size_t i = 0; i < steiner_pts.size(); i++)
		ptr->cdt.insert(steiner_pts[i]);
	ptr->obtuse = count_obtuse(&(ptr->cdt));
}

void steiner_circumcenter(struct triangulation_t *ptr, data_t *data)
{
	std::vector<Point> steiner_pts;

	for (auto it = ptr->cdt.finite_faces_begin(); it != ptr->cdt.finite_faces_end(); it++) {
		auto triangle = ptr->cdt.triangle(it);
		std::array <CGAL::Vector_2<K>, 3> edge;

		Point steiner = CGAL::circumcenter(triangle);
		if (!data->inside(steiner))
			continue;
		
		steiner_pts.push_back(steiner);
	}

	for (size_t i = 0; i < steiner_pts.size(); i++) {
		struct triangulation_t current = *ptr;
		current.cdt.insert(steiner_pts[i]);
		current.obtuse = count_obtuse(&(current.cdt));
		if (ptr->obtuse > current.obtuse)
			*ptr = current;
	}
}

void steiner_mixed(struct triangulation_t *ptr, data_t *data, unsigned int retry)
{
	struct triangulation_t original = *ptr;
	for (unsigned int i = 0; i < retry; i++) {
		struct triangulation_t current = original;
		steiner_random(&current);
		for (int j = 0; j < 20; j++) {
			struct triangulation_t current2 = current;
			steiner_circumcenter(&current2, data);
			steiner_projection(&(current2));
			if (current.obtuse > current2.obtuse)
				current = current2;
		}
		if (ptr->obtuse > current.obtuse)
			*ptr = current;
	}
}

void steiner_mixed_recursive(struct triangulation_t *ptr, data_t *data, unsigned int depth)
{
	if (depth == 0)
		return;

	struct triangulation_t best = *ptr;
	struct triangulation_t current = *ptr;
	steiner_circumcenter(&current, data);
	steiner_mixed_recursive(&current, data, depth - 1);
	if (current.obtuse == 0) {
		*ptr = current;
		return;
	}
	if (best.obtuse > current.obtuse)
		best = current;

	current = *ptr;
	steiner_projection(&current);
	steiner_mixed_recursive(&current, data, depth - 1);
	if (current.obtuse == 0) {
		*ptr = current;
		return;
	}
	if (best.obtuse > current.obtuse)
		best = current;
/*
	current = *ptr;
	steiner_random(&current);
	steiner_mixed_recursive(&current, data, depth - 1);
	if (best.obtuse > current.obtuse)
		best = current;
*/
	*ptr = best;
}

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cerr <<
			"Usage: pretty <filename>"
			<< std::endl;
		return EXIT_FAILURE;
	}

	try {
		// Parse the file as JSON
		auto const jv = parse_file(argv[1]);

		// Now pretty-print the value
		//pretty_print(std::cout, jv);

		data_in input{jv};
		//input.print();

		CDT cdt;
		data_t data{input};
		for (const auto& point : data.get_points())
			cdt.insert(point);
		for (const auto& edge : data.get_boundary())
			cdt.insert_constraint(edge.first, edge.second);
		for (const auto& constraint : data.get_constraints())
			cdt.insert_constraint(constraint.first, constraint.second);

		struct triangulation_t triangulation;
		triangulation.cdt = cdt;
		triangulation.obtuse = count_obtuse(&cdt); 

		//CGAL::make_conforming_Delaunay_2(cdt);
		//CGAL::make_conforming_Gabriel_2(cdt);
		//steiner_mixed(&cdt, &data, 5);
		using std::chrono::high_resolution_clock;
		using std::chrono::duration_cast;
		using std::chrono::duration;
		using std::chrono::milliseconds;
		auto t1 = high_resolution_clock::now();

		//steiner_mixed(&triangulation, &data, 5);
		steiner_mixed_recursive(&triangulation, &data, triangulation.obtuse);

		auto t2 = high_resolution_clock::now();
		auto ms_int = duration_cast<milliseconds>(t2 - t1);
		duration<double, std::milli> ms_double = t2 - t1;
		std::cout << ms_int.count() << "ms\n";
		std::cout << ms_double.count() << "ms\n";

		std::cout << count_obtuse(&(triangulation.cdt)) << std::endl;
		CGAL::draw(triangulation.cdt);	
	}
	catch(std::exception const& e) {
		std::cerr <<
			"Caught exception: "
			<< e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
