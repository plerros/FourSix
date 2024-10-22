#pragma once

#include "configuration.hpp"
#include "data.hpp"

enum progression_check_enum{progression_less, progression_less_equal};

class triangulation_t
{
	private:
		CDT cdt;
		size_t obtuse;
		data_t *data;
		int progression_check;

	public:
		triangulation_t(data_t *data);
		void set_progression_check(int check);

		size_t get_obtuse();
		CDT get_cdt();

		// Steiner operations
		void steiner_centroid();
		void steiner_circumcenter();
		void steiner_projection();
		void steiner_random();
		void steiner_neighbor_random();
		void steiner_mixed(unsigned int retries);
		void steiner_mixed_recursive(unsigned int depth);

		std::vector<std::pair<std::string, std::string>> get_steiner_str();
		std::vector<std::pair<size_t, size_t>> get_edges();
};