#pragma once

#include <boost/json.hpp>

class data_in
{
	private:
		std::string instance_uid;
		std::vector<std::pair<std::int64_t, std::int64_t>> points;
		std::vector<size_t> region_boundary;
		std::vector<std::pair<size_t, size_t>> constraints;

		int optim_method;

		/*
		 * delauney:
		 * True  == Start from delauney
		 * False == Start from my own triangulation
		 */
		bool delauney;

		double parameter_a;
		double parameter_b;
		double parameter_xi;
		double parameter_psi;
		double parameter_lambda;
		unsigned int parameter_kappa;
		unsigned int parameter_L;

	public:
		data_in(boost::json::value const& jv);
		void print();
		std::string get_instance_uid();
		std::vector<std::pair<std::int64_t, std::int64_t>> get_points();
		std::vector<size_t>  get_boundary();
		std::vector<std::pair<size_t, size_t>> get_constraints();

		int get_optim_method();
		bool get_delauney();
		double get_parameter_a();
		double get_parameter_b();
		double get_parameter_xi();
		double get_parameter_psi();
		double get_parameter_lambda();
		unsigned int get_parameter_kappa();
		unsigned int get_parameter_L();
};