#pragma once

#include "configuration.hpp"
#include "data.hpp"


class obtuse_polygon_t
{
	private:
		data_t *data;
		std::vector<K::Triangle_2> sequence;
		std::vector<CDT::Point> boundary;
		bool boundary_valid;

		CGAL::Polygon_2<K> polygon;
		bool polygon_valid;

		void recompute_boundary();
		void recompute_polygon();

		bool internal(CDT::Point pt);
	public:
		obtuse_polygon_t(data_t *data);

		void draw();
		CDT::Point get_steiner();
		bool has_constraint();
		bool is_convex();
		size_t size();
		bool try_insert(K::Triangle_2 triangle);
};