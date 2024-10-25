#include "configuration.hpp"

#include "helper.hpp"

bool is_in_triangle(K::Triangle_2 triangle, CDT::Point pt)
{
	for (size_t i = 0; i < 3; i++)
		if (pt == triangle.vertex(i))
			return true;
	return false;
}

bool is_obtuse(K::Triangle_2 triangle)
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