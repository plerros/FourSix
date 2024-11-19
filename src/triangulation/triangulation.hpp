#pragma once

#include "configuration.hpp"
#include "data.hpp"

enum progression_check_enum{progression_less, progression_less_equal};

/*
 * To iterate, start from st_start+1 and end at st_end
 *
 * Values outside that range shouldn't be iterated over.
 */

enum steiner_method_enum {
	st_start,
	st_circumcenter,
	st_projection_outward,
	st_projection_all,
	st_midpoint,
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
		size_t outside_obtuse;
		size_t steiner;
		data_t *data;
		int progression_check;
		std::vector<int> history;
		std::array<std::set<std::tuple<CDT::Point, CDT::Point, CDT::Point>>, st_end> *tried;

		void update_outside_obtuse();
		void insert(CDT::Point steiner, int method);
		bool exit_early();
		void steiner_projection_internal(
			std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *inward,
			std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *outward,
			const int method);

	public:
		triangulation_t(data_t *data);
		void set_progression_check(int check);

		size_t size_obtuse();
		size_t size_steiner();
		CDT get_cdt();
		std::vector<std::pair<std::string, std::string>> get_steiner_str();
		std::vector<std::pair<size_t, size_t>> get_edges();

		// Steiner operations
		void steiner_centroid(std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *steiner_pts);
		void steiner_circumcenter(std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions);
		void steiner_midpoint(std::vector<CDT::Point> *steiner_pts);
		void steiner_polygon_centroid(std::vector<CDT::Point> *steiner_pts);

		void steiner_projection(std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions);
		void steiner_projection_inward(std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions);
		void steiner_projection_outward(std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions);

		void steiner_random(std::vector<CDT::Point> *steiner_pts);
		void steiner_constraint_random(std::vector<CDT::Point> *steiner_pts);
		void steiner_neighbor_random(std::vector<std::pair<K::Triangle_2, std::vector<CDT::Point>>> *solutions);

		void steiner_add(int method);
		void steiner_add(int method, size_t max); // Add up to [max] steiner points using [method]

		// Optimization Methods
		void optim_mixed_recursive(const unsigned int depth); // My initial attempt
		void optim_local_search();
		void optim_simulated_annealing();
		void optim_ant_colony();
};
