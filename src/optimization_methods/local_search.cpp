#include "configuration.hpp"

#include "helper.hpp"
#include "triangulation.hpp"

const std::array st_lsearch = {
	st_circumcenter,
	st_projection,
	st_midpoint,
	st_centroid,
	st_polygon_centroid
};

void triangulation_t::optim_local_search(optim_alg_t parameters)
{
	const unsigned int depth = parameters.L;

	for (unsigned int i = 0; i < depth && this->obtuse > 0; i++) {
		if (PRINT_PROGRESS) {
			std::cerr << this->obtuse << " \t| ";
			std::cerr << this->steiner << " \t| ";
			if (this->progression_check == progression_less_equal)
				std::cerr << ".";
			std::cerr << i << std::endl;
		}
		triangulation_t best = *this;
		for (int j = 0; j < st_lsearch.size(); j++) {
			int method = st_lsearch[j];
			// prevent projection_all from running after projection.
			if (method == st_projection_all
				&& this->history.size() > 0
				&& this->history.back() == st_projection_outward)
				continue;
			
			triangulation_t current = *this;
			// Add 1 steiner point on each step
			current.steiner_add(method, false, 1);
			if (this->progression_check == progression_less
				&& (current.obtuse < best.obtuse))
				best = current;
			if (this->progression_check == progression_less_equal
				&& (current.obtuse <= best.obtuse))
				best = current;
		}
		if (this->obtuse > best.obtuse) {
			*this = best;
		} else {
			this->steiner_add(st_neighbor_random);
		}
	}
}