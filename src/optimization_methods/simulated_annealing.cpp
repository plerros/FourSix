#include "configuration.hpp"

#include "helper.hpp"
#include "triangulation.hpp"

const std::array st_sannealing = {
	st_circumcenter,
	st_projection,
	st_midpoint,
	st_centroid
};

void triangulation_t::optim_simulated_annealing(optim_alg_t parameters)
{
	const double obtuse_constant  = parameters.a;
	const double steiner_constant = parameters.b;
	const unsigned int depth      = parameters.L;

	double temperature = 1.0;

	while (temperature >= 0.0 && this->obtuse > 0) {
		int method = st_sannealing[std::rand() % st_sannealing.size()];

		triangulation_t current = *this;
		current.steiner_add(method);

		double delta_energy = current.get_energy(parameters) - this->get_energy(parameters);

		// Metropolis Criterion
		double random = 1.0 / ((double) std::rand());
		double probability = std::exp(- delta_energy / temperature);

		if ( delta_energy < 0 || random < probability) {
			*this = current;
		}

		temperature = temperature - (1.0 / ((double) depth));
	}
}