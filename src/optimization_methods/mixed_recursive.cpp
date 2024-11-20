#include "configuration.hpp"

#include "helper.hpp"
#include "triangulation.hpp"

void triangulation_t::optim_mixed_recursive(const unsigned int depth)
{
	if (depth == 0)
		return;
	if (this->exit_early())
		return;

	if (PRINT_PROGRESS) {
		std::cout << this->obtuse << " \t| ";
		std::cout << this->steiner << " \t|";
		cout_space(depth);
		if (this->progression_check == progression_less_equal)
			std::cout << ".";
		std::cout << depth << std::endl;
	}
	
	triangulation_t best = *this;
	triangulation_t current = *this;

	for (int method = st_start + 1; method < st_end; method++) {
		//if (this->history.size() > 0 && this->history.back() == method)
		//	continue;

		current = *this;

		bool skip_flag = false;

		switch(method) {
			// allow statements to fall through
			case st_midpoint:
			case st_polygon_centroid:
			case st_centroid:
			case st_constraint_random:
			case st_neighbor_random:
				if (this->progression_check == progression_less
					&& (best.obtuse != this->obtuse))
					skip_flag = true;
				break;
			case st_projection_all:
				// prevent projection_all from running after projection.
				if (this->history.size() > 0  && this->history.back() == st_projection_outward)
					skip_flag = true;
				break;
			default:
				break;
		}
		if (!skip_flag)
			current.steiner_add(method);

		//this->tried = current.tried;

		if (this->progression_check == progression_less
			&& (current.obtuse < best.obtuse))
			best = current;
		if (this->progression_check == progression_less_equal
			&& (current.obtuse <= best.obtuse))
			best = current;

		if (best.exit_early())
			break;

		if (current.cdt != this->cdt) {
			if (current.obtuse < this->obtuse)
				current.optim_mixed_recursive(depth - 1);
			if (current.obtuse == this->obtuse)
				current.optim_mixed_recursive(depth);
			if (this->progression_check == progression_less
				&& (current.obtuse < best.obtuse))
				best = current;
			if (this->progression_check == progression_less_equal
				&& (current.obtuse <= best.obtuse))
				best = current;
			if (best.exit_early())
				break;
		}
	}
	//best.tried = this->tried;

	*this = best;
}