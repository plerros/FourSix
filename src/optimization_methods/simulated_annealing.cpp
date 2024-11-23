#include "configuration.hpp"

#include "helper.hpp"
#include "triangulation.hpp"

const std::array st_sannealing = {
	st_circumcenter,
	st_projection,
	st_midpoint,
	st_centroid
};


static inline std::tuple<CDT::Point, CDT::Point, CDT::Point> triangle_to_tuple(K::Triangle_2 triangle)
{
	std::tuple<CDT::Point, CDT::Point, CDT::Point> ret(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));
	return ret;
}

static inline bool detect_nottried(
	CDT *cdt,
	std::set<std::tuple<CDT::Point, CDT::Point, CDT::Point>> *tried,
	CDT::Face_iterator *it,
	K::Triangle_2 triangle)
{
	//return true;
	if (tried->find(triangle_to_tuple(triangle)) == tried->end())
		return true;

	for (size_t i = 0; i < 3; i++) {
		auto neighbor = (*it)->neighbor(i);
		if (cdt->is_infinite(neighbor))
			continue;

		K::Triangle_2 tmp = cdt->triangle(neighbor);
		if (tried->find(triangle_to_tuple(tmp)) == tried->end())
			return true;
	}
	return false;
}

void triangulation_t::optim_simulated_annealing(optim_alg_t parameters)
{
	const double obtuse_constant  = parameters.a;
	const double steiner_constant = parameters.b;
	const unsigned int depth      = parameters.L;

	double temperature = 1.0;

	while (temperature >= 0.0 && this->obtuse > 0) {
/*
		for (auto it = this->cdt.finite_faces_begin(); it != this->cdt.finite_faces_end(); it++) {
			auto triangle = this->cdt.triangle(it);
			if (!this->data->inside(CGAL::centroid(triangle)))
				continue;
			if (!is_obtuse(triangle))
				continue;

			int method = st_sannealing[std::rand() % st_sannealing.size()];

			// Workaround to generate steiner points for just the current triangle
			data_t dummy_data(triangle, parameters);
			std::vector<CDT::Point> steiner_pts;

			triangulation_t tmp(&dummy_data);
			tmp.set_progression_check(progression_less_equal);
			tmp.steiner_add(method);

			// Find which point got added to the triangulation
			for (auto it = tmp.cdt.finite_vertices_begin(); it != tmp.cdt.finite_vertices_end(); it++) {
				CDT::Point pt(it->point());
				if (pt == triangle.vertex(0) || pt == triangle.vertex(1) || pt == triangle.vertex(2))
					continue;
				steiner_pts.push_back(pt);
			}

			if (steiner_pts.size() == 0)
				continue;

			triangulation_t current = (*this);
			for (size_t i = 0; i < steiner_pts.size(); i++) {
				current.insert(steiner_pts[i], method);
			}

			// Here we deviate from the alg
			double delta_energy = current.get_energy(parameters) - initial.get_energy(parameters);

			// Metropolis Criterion
			double random = 1.0 / ((double) std::rand());
			double probability = std::exp(- delta_energy / temperature);

			if ( delta_energy < 0 || random < probability) {
				*this = current;
				break;
			}

		}
*/

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