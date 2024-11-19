#include "helper.hpp"
#include "triangulation.hpp"

/*
 * Simulated annealing:
 *
 * Input
 * 	a == obtuse_constant
 * 	b == steiner_constant
 * 	L == depth
 */

void triangulation_t::optim_simulated_annealing()
{
	const double obtuse_constant = this->data->get_parameter_a();
	const double steiner_constant = this->data->get_parameter_b();
	const unsigned int depth = this->data->get_parameter_L();


}