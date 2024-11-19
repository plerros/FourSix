#include "helper.hpp"
#include "triangulation.hpp"

void triangulation_t::optim_ant_colony(optim_alg_t parameters)
{
	const double obtuse_constant     = parameters.a;
	const double steiner_constant    = parameters.b;
	const double heuristic_constant  = parameters.xi;
	const double pherormone_constant = parameters.psi;
	const double evaporation_rate    = parameters.lambda;
	const double ants_amount         = parameters.kappa;
	const unsigned int depth         = parameters.L;
}