#include "FluidSim.h"

using namespace std;
using Eigen::Vector3d;

Particle::Particle(Vector3d position) {

	// intrinsic properties
	m_position = position;
	m_delta_position = Vector3d(0, 0, 0);
	m_position_star = Vector3d(0, 0, 0);
	m_lambda = 0.d;
	m_velocity = Vector3d(0, 0, 0);
	m_acceleration = Vector3d(0, 0, 0);
	m_force = Vector3d(0, 0, 0);

	m_density = 0.d;
	m_constr_gradient = 0.d;

}
/*
 * Example simulation that changes the colors of a cube.
 */
FluidSim::FluidSim() : Simulation() {
	init();
}

void FluidSim::init() {
	// setup
	m_dt = TIME_STEP_SIZE;

	// create 3D grid of particles
	m_particleRadius = 0.25f;
	int gridSize = 5;
	double scale = 1.5;
	Vector3d offset = Vector3d(-1, 1, -1) * (scale * (gridSize - 1) / 2.0);

	for (int x = 0; x < gridSize; ++x) {
		for (int y = 0; y < gridSize; ++y) {
			for (int z = 0; z < gridSize; ++z) {
				m_particles.push_back(Particle(Vector3d(x, y, z) * scale + offset));
			}
		}
	}

	// initially calculate geometry
	updateRenderGeometry();
}

void FluidSim::resetMembers() {
	m_C.setZero();
	m_C.col(0).setOnes();
}

void FluidSim::updateRenderGeometry() {
	/// Rendering simple bounding box (cube) for each particle for now
	const int n = m_particles.size();

	// 8 vertices per particle s.t. particle position is in center of a bounding box with sides 2*radius
	int v = 8;
	m_V.resize(v * n, 3);
	for (int i = 0; i < n; ++i) {
		const Vector3d pos = m_particles[i].m_position;

		m_V.row(v*i)   = pos + Vector3d(-1, -1, -1) * m_particleRadius;
		m_V.row(v*i+1) = pos + Vector3d(-1, -1, 1) * m_particleRadius;
		m_V.row(v*i+2)   = pos + Vector3d(-1, 1, -1) * m_particleRadius;
		m_V.row(v*i+3) = pos + Vector3d(-1, 1, 1) * m_particleRadius;
		m_V.row(v*i+4) = pos + Vector3d(1, -1, -1) * m_particleRadius;
		m_V.row(v*i+5) = pos + Vector3d(1, -1, 1) * m_particleRadius;
		m_V.row(v*i+6) = pos + Vector3d(1, 1, -1) * m_particleRadius;
		m_V.row(v*i+7) = pos + Vector3d(1, 1, 1) * m_particleRadius;
	}

	// 2 triangles per face of particle bounding box
	int f = 2 * 6;
	m_F.resize(f * n, 3);
	for (int i = 0; i < n; ++i) {
		m_F.row(f*i)   = Eigen::Vector3i(v*i, v*i+1, v*i+3);
		m_F.row(f*i+1) = Eigen::Vector3i(v*i, v*i+3, v*i+2);

		m_F.row(f*i+2) = Eigen::Vector3i(v*i+1, v*i+5, v*i+7);
		m_F.row(f*i+3) = Eigen::Vector3i(v*i+1, v*i+7, v*i+3);

		m_F.row(f*i+4) = Eigen::Vector3i(v*i+5, v*i+4, v*i+6);
		m_F.row(f*i+5) = Eigen::Vector3i(v*i+5, v*i+6, v*i+7);

		m_F.row(f*i+6) = Eigen::Vector3i(v*i+4, v*i, v*i+2);
		m_F.row(f*i+7) = Eigen::Vector3i(v*i+4, v*i+2, v*i+6);

		m_F.row(f*i+8) = Eigen::Vector3i(v*i, v*i+5, v*i+1);
		m_F.row(f*i+9) = Eigen::Vector3i(v*i, v*i+4, v*i+5);

		m_F.row(f*i+10) = Eigen::Vector3i(v*i+2, v*i+3, v*i+7);
		m_F.row(f*i+11) = Eigen::Vector3i(v*i+2, v*i+7, v*i+6);
	}

	// face colors
	m_C.resize(f * n, 3);

	// update rendering
	m_renderV = m_V;
	m_renderF = m_F;
	m_renderC = m_C;
}

bool FluidSim::advance() {
	/*// do next step of some color animation
	int speed = 60;
	int decColor = (m_step / speed) % 3;
	int incColor = (decColor + 1) % 3;

	for (int i = 0; i < m_C.rows(); i++) {
		m_C(i, decColor) = (m_C(i, decColor) * speed - 1) / speed;
		m_C(i, incColor) = (m_C(i, incColor) * speed + 1) / speed;
	}

	// advance step
	m_step++;
	return false;*/

	/// Simulation step
	const int n = m_particles.size();
	// apply forces & predict position
	for (int i = 0; i < n; ++i) {
		Particle & p = m_particles[i];
		// apply gravity
		p.m_velocity += m_dt * Vector3d(0, -9.81, 0);
		// predict position
		p.m_position_star = p.m_position + m_dt * p.m_velocity;
	}
	// find neighbours

	// newton iterations
	int solverIterations = 3;
	for (int iter = 0; iter < solverIterations; ++iter) {
		// calculate lambda_i (equation 11)
		for (int i = 0; i < n; ++i) {
			Particle & pi = m_particles[i];
			// calculate density
			pi.m_density = 0;
			for (int j = 0; j < n; ++j) {
				Particle & pj = m_particles[j];
				pi.m_density += poly6Kernel(pi.m_position - pj.m_position, NEIGHBOURHOOD_RADIUS); // equation (2)

			}
			double C_i = pi.m_density / REST_DENSITY - 1.d;
			// calculate denominator, equation (11)
			double denominator = 0;

			for (int k = 0; k < n; ++k) {
				Particle & pk = m_particles[k];
				// calculate gradient d_pk_Ci, equation (8)
				double d_pk_Ci = 0;
				if (k == i) {
					double d_pk_W = 0;
					for (int j = 0; j < n; ++j) {
						Particle & pj = m_particles[j];
						d_pk_W += spikyKernel(pi.m_position - pj.m_position, NEIGHBOURHOOD_RADIUS);
					}
					d_pk_Ci = d_pk_W;
				} else {
					d_pk_Ci = - spikyKernel(pi.m_position - pk.m_position, NEIGHBOURHOOD_RADIUS);
				}
				d_pk_Ci *= 1 / REST_DENSITY;

				denominator += d_pk_Ci * d_pk_Ci;
			}
			denominator += EPSILON;

			pi.m_lambda = - C_i / denominator; // equation (11)

		}

		// calculate delta_p_i, equation (12, rsp. 14 for tensile stability)
		for (int i = 0; i < n; ++i) {
			Particle & pi = m_particles[i];
			Vector3d delta_pi = Vector3d(0, 0, 0);
			for (int j = 0; j < n; ++j) {
				Particle & pj = m_particles[j];
				//delta_pi += (pi.m_lambda + pj.m_lambda) * gradSpikyKernel(pi.m_position - pj.m_position, NEIGHBOURHOOD_RADIUS);
				const double k = 0.1;
				const Vector3d delta_q = Vector3d(0.1, 0.1, 0.1) / 3 * NEIGHBOURHOOD_RADIUS;
				double s_corr = -k * std::pow(poly6Kernel(pi.m_position - pj.m_position, NEIGHBOURHOOD_RADIUS) / poly6Kernel(delta_q, NEIGHBOURHOOD_RADIUS), 4.d);
				delta_pi += (pi.m_lambda + pj.m_lambda + s_corr) * gradSpikyKernel(pi.m_position - pj.m_position, NEIGHBOURHOOD_RADIUS);
			}
			delta_pi *= 1 / REST_DENSITY;
			pi.m_delta_position = delta_pi;
		}
		// collision detection & response
		// no collisions for now

		// update position
		for (int i = 0; i < n; ++i) {
			Particle & pi = m_particles[i];
			pi.m_position_star += pi.m_delta_position;
		}
	}

	for (int i = 0; i < n; ++i) {
		Particle & pi = m_particles[i];
		// update velocity
		pi.m_velocity = 1 / m_dt * (pi.m_position_star - pi.m_position);

		// OPTIONAL: apply vorticity confinement & XSPH viscosity
		// not for now

		// update position
		pi.m_position = pi.m_position_star;
	}

	// advance step
	m_step++;
	m_time += m_dt;
	return false;
}

void FluidSim::renderRenderGeometry(
	igl::opengl::glfw::Viewer &viewer) {
	viewer.data().set_mesh(m_renderV, m_renderF);
	viewer.data().set_colors(m_renderC);
}

// density estimation
double FluidSim::poly6Kernel(Vector3d r, double h) {
	double r_norm = r.norm();
	if (r_norm <= h) {
		return 315 / (64 * M_PI * std::pow(h, 9.d)) * std::pow(h*h - r_norm*r_norm, 3.d);
	} else {
		return 0;
	}
}

// pressure gradient
double FluidSim::spikyKernel(Vector3d r, double h) {
	double r_norm = r.norm();
	if (r_norm <= h) {
		return 15 / (M_PI * std::pow(h, 6.d)) * std::pow(h - r_norm, 3.d);
	} else {
		return 0;
	}
}

// spiky gradient
Vector3d FluidSim::gradSpikyKernel(Vector3d r, double h) {
	double r_norm = r.norm();
	if (r_norm <= h) {
		return -45 / (M_PI * std::pow(h, 6.d)) * std::pow(h - r_norm, 2.d) * r;
	} else {
		return Vector3d(0, 0, 0);
	}
}
