#include <boost/json.hpp>
#include <iterator>
#include <ostream>
#include <string>

#include "data_out.hpp"

namespace json = boost::json;

data_out::data_out(data_t *data, triangulation_t *triangulation, optim_alg_t parameters)
{
	this->content_type = "CG_SHOP_2025_Solution";
	this->instance_uid = data->get_instance_uid();

	std::vector<std::pair<std::string, std::string>> steiner_str = triangulation->get_steiner_str();
	for (auto it = steiner_str.begin(); it < steiner_str.end(); it++) {
		this->steiner_points_x.push_back(it->first);
		this->steiner_points_y.push_back(it->second);
	}
	this->edges = triangulation->get_edges();
	this->obtuse_count = triangulation->size_obtuse();
	this->parameters = parameters;
}

boost::json::value data_out::get_jsonvalue()
{
	boost::json::value ret;

	ret = {{"content_type", this->content_type}};
	ret.as_object().emplace("instance_uid", this->instance_uid);

	{
		boost::json::array arr;
		for (auto it = this->steiner_points_x.begin(); it < this->steiner_points_x.end(); it++)
			arr.emplace_back(*it);
		ret.as_object().emplace("steiner_points_x", arr);
	}
	{
		boost::json::array arr;
		for (auto it = this->steiner_points_y.begin(); it < this->steiner_points_y.end(); it++)
			arr.emplace_back(*it);
		ret.as_object().emplace("steiner_points_y", arr);
	}
	{
		boost::json::array arr;
		for (auto it = this->edges.begin(); it < this->edges.end(); it++) {
			boost::json::array arr2;
			arr2.emplace_back(it->first);
			arr2.emplace_back(it->second);
			arr.emplace_back(arr2);
		}
		ret.as_object().emplace("edges", arr);
	}

	ret.as_object().emplace("obtuse_count", this->obtuse_count);

	if (this->parameters.method == om_ls
		|| this->parameters.method == om_sa
		|| this->parameters.method == om_ant)
	{
		boost::json::object parameters;
		switch (this->parameters.method) {
			case om_ls:
				ret.as_object().emplace("method", "ls");
				parameters.emplace("L", this->parameters.L);
				break;
			case om_sa:
				ret.as_object().emplace("method", "sa");
				parameters.emplace("L", this->parameters.L);
				parameters.emplace("a", this->parameters.a);
				parameters.emplace("b", this->parameters.b);
				break;
			case om_ant:
				ret.as_object().emplace("method", "ant");
				parameters.emplace("L", this->parameters.L);
				parameters.emplace("a", this->parameters.a);
				parameters.emplace("b", this->parameters.b);
				parameters.emplace("xi", this->parameters.xi);
				parameters.emplace("psi", this->parameters.psi);
				parameters.emplace("lambda", this->parameters.lambda);
				parameters.emplace("kappa", this->parameters.kappa);
				break;
		}
		ret.as_object().emplace("parameters", parameters);

	}

	return ret;
}