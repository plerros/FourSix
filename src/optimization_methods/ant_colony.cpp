#include <cmath>

#include "configuration.hpp"

#include "helper.hpp"
#include "triangulation.hpp"

static inline size_t count_obtuse_neighbors(
	CDT *cdt,
	CDT::Face_iterator *it,
	K::Triangle_2 triangle)
{
	size_t ret = 0;
	for (size_t i = 0; i < 3; i++) {
		auto neighbor = (*it)->neighbor(i);
		if (cdt->is_infinite(neighbor))
			continue;

		K::Triangle_2 tmp = cdt->triangle(neighbor);
		if(is_obtuse(tmp))
			ret++;
	}
	return ret;
}

const std::array st_ant = {
	st_circumcenter,
	st_projection,
	st_midpoint,
	st_centroid,
	st_polygon_centroid
};

static double get_radius_to_height(K::Triangle_2 triangle)
{
	const CDT::Point circumcenter = CGAL::circumcenter(triangle);
	double circumradius = 0.0;
	{
		K::Segment_2 initial_segment(triangle.vertex(0), triangle.vertex(1));
		auto distance =  CGAL::squared_distance(circumcenter, initial_segment);

		for (int i = 1; i < 3; i++) {
			int j = (i + 1) % 3;
			K::Segment_2 seg(triangle.vertex(i), triangle.vertex(j));
			auto current = CGAL::squared_distance(circumcenter, seg);
			if (current < distance)
				distance = current;
		}
		circumradius = std::sqrt(CGAL::to_double(distance));
	}

	double shortest_height = 1.0;
	{
		K::Segment_2 longest(triangle.vertex(0), triangle.vertex(1));
		CDT::Point opposite = triangle.vertex(2);
		for (int i = 1; i < 3; i++) {
			int j = (i + 1) % 3;
			int k = (j + 1) % 3;
			K::Segment_2 seg(triangle.vertex(i), triangle.vertex(j));
			CDT::Point pt(triangle.vertex(k));
			
			if (longest.squared_length() < seg.squared_length()) {
				longest = seg;
				opposite = pt;
			}
		}
		shortest_height = std::sqrt(CGAL::to_double(CGAL::squared_distance(opposite, longest)));
	}

	return (circumradius / shortest_height);
}

static void init_heuristic(double *ptr, double radius_to_height)
{
	for (size_t i = 0; i < st_end; i++)
		ptr[i] = 0.0;

	ptr[st_circumcenter] = radius_to_height / (2.0 + radius_to_height);
	{
		double tmp = (3.0 - 2.0 * radius_to_height) / 3.0;
		if (tmp > 0.0)
			ptr[st_midpoint] = tmp;
	}
	ptr[st_polygon_centroid] = 1.0;
	{
		double tmp = (1.0 - radius_to_height) / radius_to_height;
		if (tmp > 0.0)
			ptr[st_projection] = tmp;
	}
}

static int random_method(
	K::Triangle_2 triangle,
	double *pherormones,
	const double pherormone_constant,
	const double heuristic_constant,
	size_t obtuse_neighbors)
{
	double radius_to_height = get_radius_to_height(triangle);

	double heuristic[st_end];
	init_heuristic(heuristic, radius_to_height);

	if (obtuse_neighbors > 1)
		heuristic[st_polygon_centroid] = 1.0;

	double probability[st_end]; // NOT AN ACTUAL PROBABILTY, stretching the terms here
	for (size_t i = 0; i < st_end; i++)
		probability[i] = 0.0;
	
	double sum = 0.0;
	for (size_t i = 0; i < st_ant.size(); i++) {
		int method = st_ant[i];
		probability[method] = pow(pherormones[method], pherormone_constant)
		                    + pow(heuristic[method], heuristic_constant);
		sum += probability[method];
	}

	double random = ((double) std::rand()) / ((double) RAND_MAX);
	random *= sum; // Rescale


	int method = st_end;
	for (size_t i = 0; i < st_ant.size(); i++) {
		if (random < probability[i]) {
			method = i;
			break;
		}
		random -= probability[i];
	}
	return method;
}

class ant_data_t {
	public:
		triangulation_t triangulation;
		int method;
		double delta_energy;

		ant_data_t(triangulation_t tri) : triangulation(tri), method(st_unused), delta_energy(0.0) {}
};

void triangulation_t::optim_ant_colony(optim_alg_t parameters)
{
	const double obtuse_constant     = parameters.a;
	const double steiner_constant    = parameters.b;
	const double heuristic_constant  = parameters.xi;
	const double pherormone_constant = parameters.psi;
	const double evaporation_rate    = parameters.lambda;
	const unsigned int ants_amount   = parameters.kappa;
	const unsigned int depth         = parameters.L;

	double pherormones[st_end];
	for (size_t i = 0; i < st_end; i++)
		pherormones[i] = 0.0;
	
	for (size_t i = 0; i < st_ant.size(); i++)
		pherormones[st_ant[i]] = 1.0;

	for (unsigned int cycle = 0; cycle < depth; cycle++) {
		std::vector<ant_data_t> current;
		for (unsigned int ant = 0; ant < ants_amount; ant++) {
			if (this->obtuse == 0)
				break;
			
			ant_data_t local(*this);
			// Improve triangulation (ant)

			const size_t face_count = std::distance(this->cdt.finite_faces_begin(), this->cdt.finite_faces_end());

			// TODO rewrite to improve complexity
			K::Triangle_2 triangle;
			size_t obtuse_neighbors = 0;
			while (1) {
				const size_t random_id = std::rand() % face_count;

				auto it = this->cdt.finite_faces_begin();
				std::advance(it, random_id);

				triangle = this->cdt.triangle(it);
				if (!is_obtuse(triangle))
					continue;
				obtuse_neighbors = count_obtuse_neighbors(&(this->cdt), &it, triangle);
				if (this->data->inside(CGAL::centroid(triangle)))
					break;
			}
			const int method = random_method(triangle, pherormones, pherormone_constant, heuristic_constant, obtuse_neighbors);
			local.method = method;

			// Workaround to generate steiner points for just the current triangle
			data_t dummy_data(triangle, parameters);
			std::vector<CDT::Point> steiner_pts;

			triangulation_t tmp(&dummy_data);
			tmp.set_progression_check(progression_less_equal);
			tmp.tried = NULL;
			tmp.steiner_add(method);

			// Find which point got added to the triangulation
			for (auto it = tmp.cdt.finite_vertices_begin(); it != tmp.cdt.finite_vertices_end(); it++) {
				CDT::Point pt(it->point());
				if (pt == triangle.vertex(0) || pt == triangle.vertex(1) || pt == triangle.vertex(2))
					continue;
				steiner_pts.push_back(pt);
			}
			
			for (size_t i = 0; i < steiner_pts.size(); i++)
				local.triangulation.insert(steiner_pts[i], method);

			// Evaluate triangulation (ant)
			local.delta_energy = local.triangulation.get_energy(parameters) - this->get_energy(parameters);
			current.push_back(local);
		}
		if (current.size() == 0)
			continue;

		// Save best triangulation (cycle)
		ant_data_t best_ant = current[0];
		for (size_t i = 1; i < current.size(); i++) {
			if (current[i].delta_energy < best_ant.delta_energy) {
				best_ant = current[i];
			}
		}
		(*this) = best_ant.triangulation;

		// Update pherormones (cycle)
		pherormones[best_ant.method] *= 1.0 - evaporation_rate;
		pherormones[best_ant.method] += 1.0 / (1.0 + best_ant.triangulation.get_energy(parameters));
	}
}