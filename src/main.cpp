#include <boost/json.hpp>
#include <boost/json/src.hpp>
#include <iostream>

#include "configuration.hpp"

#include "json.hpp"
#include "data_in.hpp"
#include "data.hpp"

#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Line_2.h>

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

void steiner_random(CDT *cdt)
{
	std::vector<Point> random_pts;
	{
		CDT cdt2 = *cdt;
		Mesher mesher(cdt2);
		mesher.refine_mesh();
		CGAL::Random_points_in_triangle_mesh_2<Point, CDT> generator(cdt2);
		std::copy_n(generator, 1000, std::back_inserter(random_pts));
		assert(random_pts.size() == 1000);
	}

	for (size_t i = 0; i < random_pts.size(); i++) {
		CDT cdt2 = *cdt;
		cdt2.insert(random_pts[i]);
		if (count_obtuse(&cdt2) < count_obtuse(cdt))
			cdt->insert(random_pts[i]);
	}
}

void steiner_projection(CDT *cdt)
{
	std::vector<Point> steiner_pts;

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
		cdt->insert(steiner_pts[i]);
}

void steiner_circumcenter(CDT *cdt, data_t *data)
{
	std::vector<Point> steiner_pts;

	for (auto it = cdt->finite_faces_begin(); it != cdt->finite_faces_end(); it++) {
		auto triangle = cdt->triangle(it);
		std::array <CGAL::Vector_2<K>, 3> edge;

		Point steiner = CGAL::circumcenter(triangle);
		if (!data->inside(steiner))
			continue;
		
		steiner_pts.push_back(steiner);
	}
	for (size_t i = 0; i < steiner_pts.size(); i++) {
		CDT cdt2 = *cdt;
		cdt2.insert(steiner_pts[i]);
		if (count_obtuse(cdt) > count_obtuse(&cdt2))
			*cdt = cdt2;
	}
}

void steiner_random_projection(CDT *ptr, data_t *data, unsigned int retry)
{
	CDT original = *ptr;
	for (unsigned int i = 0; i < retry; i++) {
		CDT cdt = original;
		steiner_random(&cdt);
		for (int j = 0; j < 20; j++) {
			CDT cdt2 = cdt;
			steiner_circumcenter(&cdt, data);
			steiner_projection(&cdt2);
			if (count_obtuse(&cdt) > count_obtuse(&cdt2))
				cdt = cdt2;
			//else if (count_obtuse(&cdt) < count_obtuse(&cdt2))
			//	break;
		}
		if (count_obtuse(ptr) > count_obtuse(&cdt))
			*ptr = cdt;
	}
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

		//CGAL::make_conforming_Delaunay_2(cdt);
		//CGAL::make_conforming_Gabriel_2(cdt);
		steiner_random_projection(&cdt, &data, 5);

		std::cout << count_obtuse(&cdt) << std::endl;
		CGAL::draw(cdt);	
	}
	catch(std::exception const& e) {
		std::cerr <<
			"Caught exception: "
			<< e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
