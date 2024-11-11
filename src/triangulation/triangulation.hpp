#pragma once

#include "configuration.hpp"
#include "data.hpp"

enum progression_check_enum{progression_less, progression_less_equal};

enum steiner_method_enum {
	st_start,
	st_circumcenter,
	st_projection,
	st_projection_all,
	st_polygon_centroid,
	st_centroid,
	st_constraint_random,
	st_neighbor_random,
	st_end,
	st_unused
};

class triangulation_t
{
	private:
		CDT cdt;
		size_t start_obtuse;
		size_t obtuse;
		size_t steiner;
		data_t *data;
		int progression_check;
		std::vector<int> history;
		int method_performance[st_end];

		void insert(CDT::Point steiner, int method);
		bool exit_early();
		void steiner_projection_internal(std::vector<CDT::Point> *inward, std::vector<CDT::Point> *outward);

	public:
		triangulation_t(data_t *data);
		void set_progression_check(int check);

		size_t size_obtuse();
		size_t size_steiner();
		CDT get_cdt();
		std::vector<std::pair<std::string, std::string>> get_steiner_str();
		std::vector<std::pair<size_t, size_t>> get_edges();

		// Steiner operations
		void steiner_centroid(std::vector<CDT::Point> *steiner_pts);
		void steiner_circumcenter(std::vector<CDT::Point> *steiner_pts);
		void steiner_polygon_centroid(std::vector<CDT::Point> *steiner_pts);

		void steiner_projection(std::vector<CDT::Point> *steiner_pts);
		void steiner_projection_inward(std::vector<CDT::Point> *steiner_pts);
		void steiner_projection_outward(std::vector<CDT::Point> *steiner_pts);

		void steiner_random(std::vector<CDT::Point> *steiner_pts);
		void steiner_constraint_random(std::vector<CDT::Point> *steiner_pts);
		void steiner_neighbor_random(std::vector<CDT::Point> *steiner_pts);

		void steiner_add(int method);

		void steiner_mixed(unsigned int retries);
		void steiner_mixed_recursive(unsigned int depth);

};