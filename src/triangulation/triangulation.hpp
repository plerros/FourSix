#pragma once

#include "configuration.hpp"
#include "data.hpp"

enum progression_check_enum{progression_less, progression_less_equal};
enum steiner_method_enum {st_start, st_circumcenter, st_projection, st_polygon_centroid, st_centroid, st_constraint_random, st_neighbor_random, st_end};

class triangulation_t
{
	private:
		CDT cdt;
		size_t start_obtuse;
		size_t obtuse;
		size_t steiner;
		data_t *data;
		int progression_check;
		int method_performance[st_end];

	public:
		triangulation_t(data_t *data);
		void set_progression_check(int check);

		size_t size_obtuse();
		size_t size_steiner();
		CDT get_cdt();

		// Steiner operations
		void steiner_centroid(std::vector<CDT::Point> *steiner_pts);
		void steiner_circumcenter(std::vector<CDT::Point> *steiner_pts);
		void steiner_polygon_centroid(std::vector<CDT::Point> *steiner_pts);
		void steiner_projection(std::vector<CDT::Point> *steiner_pts);
		void steiner_random(std::vector<CDT::Point> *steiner_pts);
		void steiner_constraint_random(std::vector<CDT::Point> *steiner_pts);
		void steiner_neighbor_random(std::vector<CDT::Point> *steiner_pts);
		void steiner_add(int method);
		void steiner_mixed(unsigned int retries);
		void steiner_mixed_recursive(unsigned int depth);

		void steiner_benchmark(unsigned int depth);

		std::vector<std::pair<std::string, std::string>> get_steiner_str();
		std::vector<std::pair<size_t, size_t>> get_edges();
};