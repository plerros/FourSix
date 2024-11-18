#include "helper.hpp"
#include "triangulation.hpp"

/*
 * Ant Colony:
 *
 * Input
 * 	a          == obtuse_constant
 * 	b          == steiner_constant
 * 	xi / χ     == heuristic_constant
 * 	psi / ψ    == pherormone_constant
 * 	lambda (λ) == evaporation_rate
 * 	kappa      == ants_amount
 * 	L          == depth
 */

void triangulation_t::optim_ant_colony (
	double obtuse_constant,
	double steiner_constant,
	double heuristic_constant,
	double pherormone_constant,
	double evaporation_rate,
	double ants_amount,
	unsigned int depth)
{

}