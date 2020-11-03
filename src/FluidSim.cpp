#include "FluidSim.h"

using namespace std;
using Eigen::Vector3f;

Particle::Particle(Vector3f position) {

	// intrinsic properties
	m_Position = position;
	m_Velocity = Vector3f(0, 0, 0);
	m_Acceleration = Vector3f(0, 0, 0);
	m_Force = Vector3f(0, 0, 0);

	m_Density = 0.f;

	// external forces
	m_GravitationForce = Vector3f(0, 0, 0);
}
/*
 * Example simulation that changes the colors of a cube.
 */
FluidSim::FluidSim() : Simulation() {
	init();
}

void FluidSim::init() {
	// create a cube on [-1,1]^3
	// vertices
	m_V.resize(8, 3);
	int v = 5;
	m_V << -v, -v, -v, v, -v, -v, -v, v, -v, v, v, -v, -v, -v, v, v, -v, v,
		-v, v, v, v, v, v;

	// faces
	m_F.resize(12, 3);
	m_F << 0, 2, 1, 2, 3, 1, 1, 3, 5, 3, 7, 5, 2, 6, 3, 6, 7, 3, 5, 7, 4, 7,
		6, 4, 4, 6, 0, 6, 2, 0, 0, 4, 1, 4, 5, 1;

	// face colors
	m_C.resize(12, 3);

	reset();



	// create 3D grid of particles
	int gridSize = 5;
	for (int x = 0; x < gridSize; ++x) {
		for (int y = 0; y < gridSize; ++y) {
			for (int z = 0; z < gridSize; ++z) {
				m_Particles.push_back(Particle(Eigen::Vector3f(x, y, z)));
			}
		}
	}


}

void FluidSim::resetMembers() {
	m_C.setZero();
	m_C.col(0).setOnes();
}

void FluidSim::updateRenderGeometry() {
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
