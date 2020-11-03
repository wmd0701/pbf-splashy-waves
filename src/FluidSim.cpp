#include "FluidSim.h"

using namespace std;
using Eigen::Vector3d;

Particle::Particle(Vector3d position) {

	// intrinsic properties
	m_position = position;
	m_velocity = Vector3d(0, 0, 0);
	m_acceleration = Vector3d(0, 0, 0);
	m_force = Vector3d(0, 0, 0);

	m_density = 0.f;

	// external forces
	m_gravitationForce = Vector3d(0, 0, 0);
}
/*
 * Example simulation that changes the colors of a cube.
 */
FluidSim::FluidSim() : Simulation() {
	init();
}

void FluidSim::init() {
	// create 3D grid of particles
	m_particleRadius = 0.25f;
	int gridSize = 5;
	for (int x = 0; x < gridSize; ++x) {
		for (int y = 0; y < gridSize; ++y) {
			for (int z = 0; z < gridSize; ++z) {
				m_particles.push_back(Particle(Vector3d(x, y, z)));
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
	// do next step of some color animation
	int speed = 60;
	int decColor = (m_step / speed) % 3;
	int incColor = (decColor + 1) % 3;

	for (int i = 0; i < m_C.rows(); i++) {
		m_C(i, decColor) = (m_C(i, decColor) * speed - 1) / speed;
		m_C(i, incColor) = (m_C(i, incColor) * speed + 1) / speed;
	}

	// advance step
	m_step++;
	return false;
}

void FluidSim::renderRenderGeometry(
	igl::opengl::glfw::Viewer &viewer) {
	viewer.data().set_mesh(m_renderV, m_renderF);
	viewer.data().set_colors(m_renderC);
}
