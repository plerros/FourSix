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
		std::cerr << this->obtuse << " \t| ";
		std::cerr << this->steiner << " \t|";
		cout_space(depth);
		if (this->progression_check == progression_less_equal)
			std::cerr << ".";
		std::cerr << depth << std::endl;
	}
	
	triangulation_t best = *this;
	triangulation_t current = *this;

	for (int i = 0; i < st_active.size(); i++) {
		int method = st_active[i];
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
	*this = best;
}