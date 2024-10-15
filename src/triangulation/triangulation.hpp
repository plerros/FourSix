#pragma once

#include "configuration.hpp"
#include "data.hpp"

class triangulation_t
{
	private:
		CDT cdt;
		size_t obtuse;
		data_t *data;

	public:
		triangulation_t(data_t *data);
		size_t get_obtuse();
		CDT get_cdt();

		// Steiner operations
		void steiner_centroid();
		void steiner_circumcenter(bool first_run);
		void steiner_projection();
		void steiner_random();
		void steiner_mixed(unsigned int retries);
		void steiner_mixed_recursive(unsigned int depth, bool first_run);

		std::vector<std::pair<std::string, std::string>> get_steiner_str();
		std::vector<std::pair<size_t, size_t>> get_edges();
};