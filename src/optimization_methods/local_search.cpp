#include "helper.hpp"
#include "triangulation.hpp"

/*
 * Local Search:
 *
 * Input
 * 	L == depth
 */

void triangulation_t::optim_local_search(unsigned int depth)
{
	for (; depth > 0; depth--) {
		if (PRINT_RECURSION_TREE) {
			std::cout << this->obtuse << " \t| ";
			std::cout << this->steiner << " \t|";
			cout_space(depth);
			if (this->progression_check == progression_less_equal)
				std::cout << ".";
			std::cout << depth << std::endl;
		}
		triangulation_t best = *this;
		for (int method = st_start + 1; method < st_end; method++) {
			// prevent projection_all from running after projection.
			if (method == st_projection_all
				&& this->history.size() > 0
				&& this->history.back() == st_projection)
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