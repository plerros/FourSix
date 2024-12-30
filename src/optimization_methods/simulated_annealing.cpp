#include "configuration.hpp"

#include "helper.hpp"
#include "triangulation.hpp"

const std::array st_sannealing = {
	st_circumcenter,
	st_projection,
	st_midpoint,
	st_centroid,
	st_polygon_centroid
};

#define SA_RANDOM 2

void triangulation_t::optim_simulated_annealing(optim_alg_t parameters)
{
	const double obtuse_constant  = parameters.a;
	const double steiner_constant = parameters.b;
	const unsigned int depth      = parameters.L;

	double temperature = 1.0;

	while (temperature >= 0.0 && this->obtuse > 0) {
		int method = st_sannealing[std::rand() % st_sannealing.size()];

		triangulation_t current = *this;
		/*
		 * To randomly select a triangle and insert 1 steiner point:
		 * current.steiner_add(method, true, 1);
		 */
		current.steiner_add(method);

		double delta_energy = current.get_energy(parameters) - this->get_energy(parameters);

		// Metropolis Criterion
		double random = ((double) std::rand()) / ((double) RAND_MAX);
		double probability = std::exp(- delta_energy / temperature);

#ifndef SA_RANDOM
		if (delta_energy < 0 || random < probability) {
			*this = current;
		}

#elif SA_RANDOM == 1	// attempt 1
		if (delta_energy < 0 || random < probability) {
			*this = current;
		} else {
			this->steiner_add(st_constraint_random);
			this->steiner_add(st_neighbor_random, true, 1);
		}

#elif SA_RANDOM == 2	// attempt 2
		if (delta_energy < 0) {
			*this = current;
		}
		else if (random < probability) {
			//this->steiner_add(st_constraint_random);
			this->steiner_add(st_neighbor_random);
		}
#endif
		temperature = temperature - (1.0 / ((double) depth));
	}
}