#include <boost/json.hpp>
#include <iostream>

#include "configuration.hpp"
#include "data_in.hpp"

data_in::data_in(boost::json::value const& jv)
{
	boost::json::object const& root{jv.as_object()};

	this->instance_uid = jv.at("instance_uid").as_string();
	assert(jv.at("num_points").as_int64() >= 0);

	{
		auto tmp = jv.at("points_x").as_array();
		assert(tmp.size() == jv.at("num_points").as_int64());

		if (tmp.empty())
			goto skip_points_x;

		for (auto it = tmp.begin(); it != tmp.end(); it++)
			this->points.push_back({it->as_int64(), 0});
	}
skip_points_x:

	{
		auto tmp = jv.at("points_y").as_array();
		assert(tmp.size() == jv.at("num_points").as_int64());

		if (tmp.empty())
			goto skip_points_y;

		size_t i = 0;
		for (auto it = tmp.begin(); it != tmp.end() && i < this->points.size(); it++, i++)
			this->points[i].second = it->as_int64();
	}
skip_points_y:

	{
		auto tmp = jv.at("region_boundary").as_array();
		assert(tmp.size() <= jv.at("num_points").as_int64());

		if (tmp.empty())
			goto skip_region_boundary;

		for (auto it = tmp.begin(); it != tmp.end(); it++)
			this->region_boundary.push_back(it->as_int64());
	}
skip_region_boundary:

	assert(jv.at("num_constraints").as_int64() >= 0);

	{
		auto tmp = jv.at("additional_constraints").as_array();
		assert(tmp.size() == jv.at("num_constraints").as_int64());

		if (tmp.empty())
			goto skip_additional_constraints;

		for (auto it = tmp.begin(); it != tmp.end(); it++) {
			auto arr = it->as_array();
			assert(arr[0].as_int64() >= 0);
			assert(arr[1].as_int64() >= 0);

			size_t arr0 = arr[0].as_int64();
			size_t arr1 = arr[1].as_int64();

			this->constraints.push_back({arr0, arr1});
		}
	}
skip_additional_constraints:

	this->optim_method = optim_none;
	if (root.if_contains("method")) {
		auto tmp = jv.at("method").as_string();
		if (tmp == "ls")
			this->optim_method = optim_local_search;
		if (tmp == "sa")
			this->optim_method = optim_simulated_annealing;
		if (tmp == "ant")
			this->optim_method = optim_ant_colony;
	}

	this->delauney = false;
	if (root.if_contains("delauney")) {
		auto tmp = jv.at("delauney").as_bool();
		this->delauney=tmp;
	}

	this->parameter_a      = 1.0;
	this->parameter_b      = 1.0;
	this->parameter_xi     = 0.0;
	this->parameter_psi    = 0.0;
	this->parameter_lambda = 0.0;
	this->parameter_kappa  = 0;
	this->parameter_L      = 0;

	if (root.if_contains("parameters")) {
		auto parameters = jv.at("parameters").as_object();
		auto pjv = jv.at("parameters");

		if (parameters.if_contains("alpha")) {
			auto tmp = pjv.at("alpha").as_double();
			this->parameter_a = tmp;
		}
		
		if (parameters.if_contains("beta")) {
			auto tmp = pjv.at("beta").as_double();
			this->parameter_b = tmp;
		}

		if (parameters.if_contains("xi")) {
			auto tmp = pjv.at("xi").as_double();
			this->parameter_xi = tmp;
		}

		if (parameters.if_contains("psi")) {
			auto tmp = pjv.at("psi").as_double();
			this->parameter_psi = tmp;
		}

		if (parameters.if_contains("lambda")) {
			auto tmp = pjv.at("lambda").as_double();
			this->parameter_lambda = tmp;
		}

		if (parameters.if_contains("kappa")) {
			assert(pjv.at("kappa").as_int64() >= 0);
			auto tmp = pjv.at("kappa").as_int64();
			this->parameter_kappa = tmp;
		}

		if (parameters.if_contains("L")) {
			assert(pjv.at("L").as_int64() >= 0);
			auto tmp = pjv.at("L").as_int64();
			this->parameter_L = tmp;
		}
	}

	if (1) {}; // Needed for previous label
}

void data_in::print()
{
	std::cout << this->instance_uid << "\n";

	for (const auto& point : this->points)
		std::cout << "{" << point.first << ", " << point.second << "}\n";
	for (const auto& point : this->region_boundary)
		std::cout << point << "\n";
}

std::string data_in::get_instance_uid()
{
	return this->instance_uid;
}

std::vector<std::pair<std::int64_t, std::int64_t>> data_in::get_points()
{
	return this->points;
}

std::vector<size_t>  data_in::get_boundary()
{
	return this->region_boundary;
}

std::vector<std::pair<size_t, size_t>> data_in::get_constraints()
{
	return this->constraints;
}

int data_in::get_optim_method()
{
	return this->optim_method;
}

bool data_in::get_delauney()
{
	return this->delauney;
}

double data_in::get_parameter_a()
{
	return this->parameter_a;
}

double data_in::get_parameter_b()
{
	return this->parameter_b;
}

double data_in::get_parameter_xi()
{
	return this->parameter_xi;
}

double data_in::get_parameter_psi()
{
	return this->parameter_psi;
}

double data_in::get_parameter_lambda()
{
	return this->parameter_lambda;
}

unsigned int data_in::get_parameter_kappa()
{
	return this->parameter_kappa;
}

unsigned int data_in::get_parameter_L()
{
	return this->parameter_L;
}
