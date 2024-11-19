#include "helper.hpp"
#include "triangulation.hpp"

void triangulation_t::optim_ant_colony()
{
	const double obtuse_constant     = this->data->get_parameter_a();
	const double steiner_constant    = this->data->get_parameter_b();
	const double heuristic_constant  = this->data->get_parameter_xi();
	const double pherormone_constant = this->data->get_parameter_psi();
	const double evaporation_rate    = this->data->get_parameter_lambda();
	const double ants_amount         = this->data->get_parameter_kappa();
	const unsigned int depth         = this->data->get_parameter_L();
}