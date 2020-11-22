#include "FluidSim.h"

using namespace std;

FluidSim::FluidSim() : Simulation() {
	init();
}

void FluidSim::init() {
	// just draw a ground plane to satisfy igl... with no geometry it seems to crash!
	// vertices
	float floorHalfSize = PARTICLES_PER_CUBE_SIDE * PARTICLE_RADIUS * 2.0f;
	m_renderV.resize(4, 3);
	m_renderV << -floorHalfSize, 0.0f, -floorHalfSize,
				  floorHalfSize, 0.0f, -floorHalfSize,
		          floorHalfSize, 0.0f, floorHalfSize,
			     -floorHalfSize, 0.0f, floorHalfSize;
	// faces
	m_renderF.resize(2, 3);
	m_renderF << 2, 1, 0, 3, 2, 0;
	// face colors
	m_renderC.resize(2, 3);
	m_renderC << 1.0f, 1.0f, 1.0f, 
		1.0f, 1.0f, 1.0f;
	// uv coordinates for grid texture
	m_renderUV.resize(4, 2);
	m_renderUV << 0, 0, 
		2.0f * floorHalfSize, 0, 
		2.0f * floorHalfSize, 2.0f * floorHalfSize, 
		0, 2.0f * floorHalfSize;

	// setup
	m_dt = TIME_STEP_SIZE;

	// populates the particle cube and all particles initial data.
	resetMembers();

	// The "fast" instanced viewer needs to be initialized once opengl is initialized.
	// There is no other point for unobtrusive initialization that I found but in renderRenderGeometry.
	initializedInstancedViewer = false;
	p_iviewer = new InstancedViewer(positions, renderColors);
}

void FluidSim::resetMembers() {

	// repopulate the particle cube, reset colors and all other particle data.
	int n = PARTICLES_PER_CUBE_SIDE;
	positions1.resize(n * n * n, 3);
	positions2.resize(n * n * n, 3);
	positionsStar = &positions1;
	positions = &positions2;

	colors1.resize(n * n * n);
	colors2.resize(n * n * n);
	updateColors = &colors1;
	renderColors = &colors2;

	float particleDiameter = 2.0f * PARTICLE_RADIUS;
	float positionOffset = (n * particleDiameter) / 2.0f + PARTICLE_RADIUS;
	float x = -positionOffset;
	for (int i = 0; i < n; i++)
	{
		x += particleDiameter;
		float y = positionOffset;
		for (int j = 0; j < n; j++)
		{
			y += particleDiameter;
			float z = -positionOffset;
			for (int k = 0; k < n; k++)
			{
				z += particleDiameter;
				positions1.row(i*n*n+j*n+k) << x, y, z;
				colors1[i * n * n + j * n + k] = (float)(i * n * n + j * n) / (float)(n * n * n);
			}
		}
	}

	colors2 = colors1;
	positions2 = positions1;

	velocities = Eigen::Matrix<float, -1, -1, Eigen::RowMajor>::Zero(n * n * n, 3);
	densities = Eigen::VectorXf::Zero(n * n * n);
	lambdas = Eigen::VectorXf::Zero(n * n * n);
}

void FluidSim::updateRenderGeometry() {
	// just swap the pointers, no copying or even rewriting of all the data (it's done on the gpu)
	auto temp = positions;
	auto temp2 = renderColors;
	positions = positionsStar;
	positionsStar = temp;
	renderColors = updateColors;
	updateColors = temp2;
}

bool FluidSim::advance() {
	/// Simulation step
	const int n = positions->rows();
	// delta q from equation (13)
	const Eigen::Vector3f delta_q = Eigen::Vector3f(0.1, 0.1, 0.1) / 3 * NEIGHBOURHOOD_RADIUS;
	const float poly6kernel_delta_q = poly6Kernel(delta_q, NEIGHBOURHOOD_RADIUS);
	// apply forces & predict position
	for (int i = 0; i < n; ++i) {
		// apply gravity
		velocities.row(i) += m_dt * Eigen::Vector3f(0.0f, -9.81f, 0.0f);
		// predict position
		positionsStar->row(i) = positions->row(i) + m_dt * velocities.row(i);
	}
	// TODO: find neighbours
	
	// newton iterations
	int solverIterations = 3;
	for (int iter = 0; iter < solverIterations; ++iter) {
		// calculate lambda_i (equation 11)
		for (int i = 0; i < n; ++i) {
			// calculate density
			densities[i] = 0.0f;
			const Eigen::RowVector3f& position_i = positionsStar->row(i);
			for (int j = 0; j < n; ++j) {
				densities[i] += poly6Kernel(position_i - positionsStar->row(j), NEIGHBOURHOOD_RADIUS); // equation (2)
			}
			float C_i = densities[i] / REST_DENSITY - 1.f;
			// calculate denominator, equation (11)
			float denominator = 0.0f;

			for (int k = 0; k < n; ++k) {
				// calculate gradient d_pk_Ci, equation (8)
				Eigen::Vector3f d_pk_Ci = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
				if (k == i) {
					for (int j = 0; j < n; ++j) {
						d_pk_Ci += gradSpikyKernel(position_i - positionsStar->row(j), NEIGHBOURHOOD_RADIUS);
					}
				} else {
					d_pk_Ci = -1.0f * gradSpikyKernel(position_i - positionsStar->row(k), NEIGHBOURHOOD_RADIUS);
				}
				d_pk_Ci *= 1.0f / REST_DENSITY;

				denominator += d_pk_Ci.squaredNorm();
			}
			denominator += EPSILON;
			lambdas[i] = - C_i / denominator; // equation (11)
		}
		// calculate delta_p_i, equation (12, rsp. 14 for tensile stability)
		for (int i = 0; i < n; ++i) {
			Eigen::Vector3f delta_pi = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
			for (int j = 0; j < n; ++j) {
				//delta_pi += (lambdas[i] + lambdas[j]) * gradSpikyKernel(renderPositions->row(i) - renderPositions->row(j), NEIGHBOURHOOD_RADIUS); // equation (12)
				const float k = 0.1f;
				float s_corr = -k * std::pow(poly6Kernel(positionsStar->row(i) - positionsStar->row(j), NEIGHBOURHOOD_RADIUS) / poly6kernel_delta_q, 4); // equation (13)
				delta_pi += (lambdas[i] + lambdas[j] + s_corr) * gradSpikyKernel(positionsStar->row(i) - positionsStar->row(j), NEIGHBOURHOOD_RADIUS); // equation (14)
			}
			delta_pi *= 1 / REST_DENSITY;
			positionsStar->row(i) += delta_pi;
		}
		// collision detection & response
		// no collisions for now
	}

	for (int i = 0; i < n; ++i) {
		// update velocity
		velocities.row(i) = 1 / m_dt * (positionsStar->row(i) - positions->row(i));

		// OPTIONAL: apply vorticity confinement & XSPH viscosity
		// not for now

		// update position (already happens automatically, due to double buffering)
	}

	// advance step
	m_step++;
	m_time += m_dt;
	return false;
}

void FluidSim::renderRenderGeometry(
	igl::opengl::glfw::Viewer &viewer) {
	//just keeping it here at every frame because its necessary, or it crashes the igl viewer.
	viewer.data().set_mesh(m_renderV, m_renderF);

	if (!initializedInstancedViewer) { 
		/* do initialization things once
		 * It is required to be initialized here, as otherwise opengl is not initialized.
		 * To do it elsewhere, we need to rewrite the gui/simulator/simulation to adjust for that.
		 * initializing here is therefore ugly, but requires no changes to other (library) code.
		 */

		//viewer.core.background_color = Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
		viewer.data().grid_texture();
		viewer.data().set_colors(m_renderC);
		viewer.data().set_uv(m_renderUV);
		viewer.data().show_texture = true;
		viewer.data().shininess = 0.5;
		viewer.core.animation_max_fps = 60; // still seems to ceil at 30 no matter what.

		p_iviewer->init(&viewer);
		p_iviewer->setParticleSize(PARTICLE_RADIUS);
		//p_iviewer->setPerInstanceColor(false);
		initializedInstancedViewer = true;
	}

	p_iviewer->updatePositions(positions);
	//p_iviewer->updateColors(renderColors);
	p_iviewer->drawInstanced();
}

// density estimation
float FluidSim::poly6Kernel(Eigen::Vector3f r, float h) {
	float r_norm = r.norm();
	if (r_norm <= h) {
		return 315.f / (64.f * M_PI * std::pow(h, 9)) * std::pow(h*h - r_norm*r_norm, 3);
	} else {
		return 0;
	}
}

// pressure gradient
float FluidSim::spikyKernel(Eigen::Vector3f r, float h) {
	float r_norm = r.norm();
	if (r_norm <= h) {
		return 15.f / (M_PI * std::pow(h, 6)) * std::pow(h - r_norm, 3);
	} else {
		return 0;
	}
}

// spiky gradient
Eigen::Vector3f FluidSim::gradSpikyKernel(Eigen::Vector3f r, float h) {
	float r_norm = r.norm();
	if (1e-12 <= r_norm && r_norm <= h) {
		return -45.f / (M_PI * std::pow(h, 6)) * std::pow(h - r_norm, 2) * r / r_norm;
	} else {
		return Eigen::Vector3f(0.0f, 0.0f, 0.0f);
	}
}
