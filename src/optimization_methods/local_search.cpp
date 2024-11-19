#include "helper.hpp"
#include "triangulation.hpp"

void triangulation_t::optim_local_search()
{
	const unsigned int depth = this->data->get_parameter_L();

	for (unsigned int i = 0; i < depth && this->obtuse > 0; i++) {
		if (PRINT_PROGRESS) {
			std::cout << this->obtuse << " \t| ";
			std::cout << this->steiner << " \t| ";
			if (this->progression_check == progression_less_equal)
				std::cout << ".";
			std::cout << i << std::endl;
		}
		triangulation_t best = *this;
		for (int method = st_start + 1; method < st_end; method++) {
			// prevent projection_all from running after projection.
			if (method == st_projection_all
				&& this->history.size() > 0
				&& this->history.back() == st_projection_outward)
				continue;
			
			triangulation_t current = *this;
			// Add 1 steiner point on each step
			current.steiner_add(method, 1);
			if (this->progression_check == progression_less
				&& (current.obtuse < best.obtuse))
				best = current;
			if (this->progression_check == progression_less_equal
				&& (current.obtuse <= best.obtuse))
				best = current;
			}
		*this = best;
	}
}