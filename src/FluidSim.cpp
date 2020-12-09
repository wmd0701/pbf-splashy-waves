#include "FluidSim.h"

using namespace std;

// TODO: do not have them as global variables, and properly terminate the thread function if possible.
bool finished = false;
Barrier barrier(ThreadCount);
Barrier start_barrier(ThreadCount + 1);
Barrier end_barrier(ThreadCount + 1);
Barrier preprefix_barrier(ThreadCount + 1);
Barrier prefix_barrier(ThreadCount + 1);

// The actual Advance method when computing with the threads
void perThreadAdvance(int start_index, int end_index, FluidSim* fs)
{
	while (!finished) // TODO: add proper termination (now its infinite -> reset might not work)
	{
		// reset members for neighborhood search
		start_barrier.wait();

		// delta q from equation (13)
		const Eigen::Vector3f delta_q = Eigen::Vector3f(0.1, 0.1, 0.1) / 3 * NEIGHBOURHOOD_RADIUS;
		const float poly6kernel_delta_q = fs->poly6Kernel(delta_q);
		// apply forces & predict position

		for (int i = start_index; i < end_index; i++) {
			// apply gravity
			fs->velocities->row(i) += fs->m_dt * Eigen::Vector3f(0.0f, -9.81f, 0.0f);
			// predict position
			fs->positionsStar->row(i) = fs->positions->row(i) + fs->m_dt * fs->velocities->row(i);
		}

		for (int i = start_index; i < end_index; i++)
		{
			const Eigen::Vector3f& pos = fs->positionsStar->row(i);
			int gridPosX = std::floor(fs->relPos * (pos.x() + fs->posOffsetxz));
			int gridPosY = std::floor(fs->relPos * (pos.y() + fs->posOffsety));
			int gridPosZ = std::floor(fs->relPos * (pos.z() + fs->posOffsetxz));
			int entryPos = gridPosX + gridPosY * fs->gridWidth + gridPosZ * fs->gridWidth * fs->gridWidth;
			fs->bin_index[i] = entryPos;
			fs->bin_sub_index[i] = fs->bin_count[entryPos].fetch_add(1); // out of range value!!!!
		}

		preprefix_barrier.wait();
		prefix_barrier.wait();

		//reinsert particles according to prefix sum and sub index.
		for (int i = start_index; i < end_index; i++)
		{
			int starting_location = fs->bin_prefix_sum[fs->bin_index[i]];
			int offset = fs->bin_sub_index[i];
			fs->neighbor_bin_index[starting_location + offset] = i;
		}
		// neighbor search is done by now.

		barrier.wait();

		// newton iterations
		int solverIterations = 1;
		for (int iter = 0; iter < solverIterations; ++iter) {

			// calculate lambda_i (equation 11)
			for (int i = start_index; i < end_index; ++i) {
				// calculate density
				fs->densities[i] = 0.0f;
				const Eigen::RowVector3f& pos_i = fs->positionsStar->row(i);
				int p_bin = fs->bin_index[i];

				// calculate denominator, equation (11)
				float denominator = 0.0f;

				// all 27 neighboring cells...
				Eigen::Vector3f d_pk_Ci = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
				Eigen::Vector3f d_pk_Ci_same = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
				for (int z = -1; z < 2; z++)
				{
					for (int y = -1; y < 2; y++)
					{
						int actual_bin = p_bin - 1 + y * fs->gridWidth + z * fs->gridWidth * fs->gridWidth;
						int bin_start = fs->bin_prefix_sum[actual_bin];
						int bin_end = fs->bin_prefix_sum[actual_bin + 1 + 2];
						for (int n = bin_start; n < bin_end; n++) // for each neighbor
						{
							int n_index = fs->neighbor_bin_index[n];
							const Eigen::RowVector3f& pos_j = fs->positionsStar->row(n_index);

							if (n_index != i) // omit same particle
							{
								fs->densities[i] += fs->poly6Kernel(pos_i - pos_j); // equation (2)
								d_pk_Ci += fs->gradSpikyKernel(pos_i - pos_j);
								d_pk_Ci_same = -1.0f * fs->gradSpikyKernel(pos_i - pos_j);
								d_pk_Ci_same *= 1.0f / REST_DENSITY;
								denominator += d_pk_Ci_same.squaredNorm();

							}
						}
					}
				}

				d_pk_Ci *= 1.0f / REST_DENSITY;
				denominator += d_pk_Ci.squaredNorm();
				float C_i = fs->densities[i] / REST_DENSITY - 1.f;

				denominator += EPSILON;
				fs->lambdas[i] = -C_i / denominator; // equation (11)


				// calculate delta_p_i, equation (12, rsp. 14 for tensile stability)
				Eigen::Vector3f delta_pi = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
				for (int z = -1; z < 2; z++)
				{
					for (int y = -1; y < 2; y++)
					{
						int actual_bin = p_bin - 1 + y * fs->gridWidth + z * fs->gridWidth * fs->gridWidth;
						int bin_start = fs->bin_prefix_sum[actual_bin];
						int bin_end = fs->bin_prefix_sum[actual_bin + 1 + 2];
						for (int n = bin_start; n < bin_end; n++) // for each neighbor
						{
							int n_index = fs->neighbor_bin_index[n];
							const Eigen::RowVector3f& pos_j = fs->positionsStar->row(n_index);

							if (n_index != i) // omit same particle
							{
								delta_pi += (fs->lambdas[i] + fs->lambdas[n_index]) * fs->gradSpikyKernel(pos_i - pos_j); // equation (12)
							}
						}
					}
				}
				delta_pi *= 1 / REST_DENSITY;

				// boundary
				Eigen::RowVector3f contactPoint;
				Eigen::Vector3f normal; // unit surface normal
				// TODO: DANGEROUS -> too high velocities will push particles out of border too far -> there are no grids for neighbors in those areas and program will crash
				// a very basic solution could be to use a smaller time step
				if (fs->collision(pos_i, contactPoint, normal)) {
					//std::cout << "collision" << std::endl;
					const Eigen::Vector3f velocity = fs->velocities->row(i);
					float penetrationDepth = (pos_i - contactPoint).norm();
					fs->velocities->row(i) = velocity - normal * (1 + 0.5f * penetrationDepth / (fs->m_dt * velocity.norm())) * velocity.dot(normal);
					fs->positionsStar->row(i) = contactPoint;
				}

				// update position
				fs->positionsStar->row(i) += delta_pi;

				// An idea for even faster collision detection... directly project to valid position... but requires more solver iterations in order to not break stuff
				/*
				const Eigen::RowVector3f& poss_i = fs->positionsStar->row(i);
				const float bouncing = -0.8f;
				if (poss_i.x() < fs->simBoundary[0]) { (*fs->positionsStar)(i, 0) = fs->simBoundary[0]; (*fs->velocities)(i, 0) *= bouncing; }
				if (poss_i.x() > fs->simBoundary[1]) { (*fs->positionsStar)(i, 0) = fs->simBoundary[1]; (*fs->velocities)(i, 0) *= bouncing; }
				if (poss_i.y() < fs->simBoundary[2]) { (*fs->positionsStar)(i, 1) = fs->simBoundary[2]; (*fs->velocities)(i, 1) *= bouncing; }
				if (poss_i.y() > fs->simBoundary[3]) { (*fs->positionsStar)(i, 1) = fs->simBoundary[3]; (*fs->velocities)(i, 1) *= bouncing; }
				if (poss_i.z() < fs->simBoundary[4]) { (*fs->positionsStar)(i, 2) = fs->simBoundary[4]; (*fs->velocities)(i, 2) *= bouncing; }
				if (poss_i.z() > fs->simBoundary[5]) { (*fs->positionsStar)(i, 2) = fs->simBoundary[5]; (*fs->velocities)(i, 2) *= bouncing; }
				*/
			}
		}

		for (int i = start_index; i < end_index; ++i) {
			// update velocity
			fs->velocities->row(i) = 1 / fs->m_dt * (fs->positionsStar->row(i) - fs->positions->row(i));
			
			// OPTIONAL: apply vorticity confinement
			// not for now

			// OPTIONAL: apply XSPH viscosity, equation (17)
			const float c = 0.01f;
			fs->velocitiesStar->row(i) = fs->velocities->row(i);
			int p_bin = fs->bin_index[i];
			for (int z = -1; z < 2; z++)
			{
				for (int y = -1; y < 2; y++)
				{
					int actual_bin = p_bin - 1 + y * fs->gridWidth + z * fs->gridWidth * fs->gridWidth;
					int bin_start = fs->bin_prefix_sum[actual_bin];
					int bin_end = fs->bin_prefix_sum[actual_bin + 1 + 2];
					for (int n = bin_start; n < bin_end; n++) // for each neighbor
					{
						int n_index = fs->neighbor_bin_index[n];
						const Eigen::Vector3f& v_ij= fs->velocities->row(n_index) - fs->velocities->row(i);
						fs->velocitiesStar->row(i) += c * v_ij * fs->poly6Kernel(fs->positionsStar->row(i) - fs->positionsStar->row(n_index));
					}
				}
			}

		}

		end_barrier.wait();
	}
}

FluidSim::FluidSim() : Simulation() {
	init();
}

void FluidSim::init() {
	
	simBoundary.resize(6);
	simBoundary << -halfBoundarySize, halfBoundarySize,
					0.0f, 2.0f * halfBoundarySize,
					-halfBoundarySize, halfBoundarySize;

	// just draw a ground plane to satisfy igl... with no geometry it seems to crash!
	// vertices
	float floorHalfSize = halfBoundarySize;
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
	p_iviewer = new InstancedViewer(positions, renderColors, NUM_FLUID_PARTICLES);



	// for threaded advance, sets range for each thread and starts the "pool"
	indices = std::vector<int>(ThreadCount + 1, NUM_FLUID_PARTICLES);
	indices[0] = 0;
	int amount = NUM_FLUID_PARTICLES / ThreadCount; // work per thread
	int upto = amount;
	for (int i = 1; i < indices.size() - 1; i++)
	{
		indices[i] = upto;
		upto += amount;
	}
	indices[indices.size() - 1] = NUM_FLUID_PARTICLES;

	for (int i = 0; i < ThreadCount; i++)
	{
		threads.push_back(std::thread(perThreadAdvance, indices[i], indices[i + 1], this));
	}
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
				positions1.row(i * n * n + j * n + k) << x, y, z;
				colors1[i * n * n + j * n + k] = (float)(i * n * n + j * n) / (float)(n * n * n);
			}
		}
	}

	colors2 = colors1;
	positions2 = positions1;

	velocities1 = Eigen::Matrix<float, -1, -1, Eigen::RowMajor>::Zero(n * n * n, 3);
	velocities2 = Eigen::Matrix<float, -1, -1, Eigen::RowMajor>::Zero(n * n * n, 3);
	velocitiesStar = &velocities1;
	velocities = &velocities2;
	densities = Eigen::VectorXf::Zero(n * n * n);
	lambdas = Eigen::VectorXf::Zero(n * n * n);

	// for fast neighborhood search
	int totalGridCells = gridWidth * gridWidth * gridWidth;
	bin_index = std::vector<unsigned int>(NUM_FLUID_PARTICLES);
	bin_sub_index = std::vector<unsigned int>(NUM_FLUID_PARTICLES, 0);
	//bin_count = std::vector<unsigned int>(totalGridCells, 0);
	bin_count = std::vector<std::atomic<unsigned int>>(totalGridCells);
	for (int i = 0; i < totalGridCells; i++) bin_count[i].store(0);
	bin_prefix_sum = std::vector<unsigned int>(totalGridCells, 0);
	// actually stores the neighbor indices in sorted grid
	neighbor_bin_index = std::vector<unsigned int>(NUM_FLUID_PARTICLES);
}

void FluidSim::updateRenderGeometry() {
	// just swap the pointers, no copying or even rewriting of all the data (it's done on the gpu)
	auto temp = positions;
	auto temp2 = renderColors;
	auto temp3 = velocities;
	positions = positionsStar;
	positionsStar = temp;
	renderColors = updateColors;
	updateColors = temp2;
	velocities = velocitiesStar;
	velocitiesStar = temp3;
}

// Multithreaded advance method...
// TODO: possibly parallelize prefix sum and resetting of vectors if it helps
bool FluidSim::advance()
{
	std::fill(bin_sub_index.begin(), bin_sub_index.end(), 0);
	//std::fill(bin_count.begin(), bin_count.end(), 0);
	for (int i = 0; i < gridWidth * gridWidth * gridWidth; i++) bin_count[i].store(0);
	std::fill(bin_prefix_sum.begin(), bin_prefix_sum.end(), 0);

	// do stuff in parallel
	start_barrier.wait(); // should get everything to start...

	preprefix_barrier.wait(); // prefix sum is only thing non-parallel
	int sum = 0;
	for (int i = 0; i < gridWidth * gridWidth * gridWidth; i++)
	{
		bin_prefix_sum[i] = sum;
		sum += bin_count[i];
	}
	prefix_barrier.wait(); // gets the others running again...

	end_barrier.wait(); // waits until the others are finished

	// advance step
	m_step++;
	m_time += m_dt;
	return false;
}
// Single Threaded Advance method
// ATTENTION: requires the bin_count array to be made to a non-atomic variant in FluidSim.h and in the reset_members function.
/*
bool FluidSim::advance() {

	/// Simulation step
	const int n = positions->rows();
	// delta q from equation (13)
	const Eigen::Vector3f delta_q = Eigen::Vector3f(0.1, 0.1, 0.1) / 3 * NEIGHBOURHOOD_RADIUS;
	const float poly6kernel_delta_q = poly6Kernel(delta_q);
	// apply forces & predict position

	for (int i = 0; i < n; ++i) {
		// apply gravity
		velocities.row(i) += m_dt * Eigen::Vector3f(0.0f, -9.81f, 0.0f);
		// predict position
		positionsStar->row(i) = positions->row(i) + m_dt * velocities.row(i);
	}


	// reset members for neighborhood search
	std::fill(bin_sub_index.begin(), bin_sub_index.end(), 0);
	std::fill(bin_count.begin(), bin_count.end(), 0);
	std::fill(bin_prefix_sum.begin(), bin_prefix_sum.end(), 0);

	for (int i = 0; i < PARTICLES_PER_CUBE_SIDE*PARTICLES_PER_CUBE_SIDE*PARTICLES_PER_CUBE_SIDE; i++)
	{
		const Eigen::Vector3f& pos = positions->row(i);
		int gridPosX = std::floor(relPos * (pos.x() + posOffsetxz));
		int gridPosY = std::floor(relPos * (pos.y() + posOffsety));
		int gridPosZ = std::floor(relPos * (pos.z() + posOffsetxz));
		int entryPos = gridPosX + gridPosY * gridWidth + gridPosZ * gridWidth * gridWidth;
		bin_index[i] = entryPos;
		bin_sub_index[i] = bin_count[entryPos]++;
	}

	// compute prefix sum
	int sum = 0;
	for (int i = 0; i < bin_prefix_sum.size(); i++)
	{
		bin_prefix_sum[i] = sum;
		sum += bin_count[i];
	}

	//reinsert particles according to prefix sum and sub index.
	for (int i = 0; i < PARTICLES_PER_CUBE_SIDE*PARTICLES_PER_CUBE_SIDE*PARTICLES_PER_CUBE_SIDE; i++)
	{
		int starting_location = bin_prefix_sum[bin_index[i]];
		int offset = bin_sub_index[i];
		neighbor_bin_index[starting_location + offset] = i;
	}
	// neighbor search is done by now.



	// newton iterations
	int solverIterations = 1;
	for (int iter = 0; iter < solverIterations; ++iter) {

		// calculate lambda_i (equation 11)
		for (int i = 0; i < n; ++i) {
			// calculate density
			densities[i] = 0.0f;
			const Eigen::RowVector3f& pos_i = positionsStar->row(i);
			int p_bin = bin_index[i];

			// calculate denominator, equation (11)
			float denominator = 0.0f;

			// all 27 neighboring cells...
			Eigen::Vector3f d_pk_Ci = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
			Eigen::Vector3f d_pk_Ci_same = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
			for (int z = -1; z < 2; z++)
			{
				for (int y = -1; y < 2; y++)
				{
					int actual_bin = p_bin - 1 + y * gridWidth + z * gridWidth * gridWidth;
					int bin_start = bin_prefix_sum[actual_bin];
					int bin_end = bin_prefix_sum[actual_bin + 1 + 2];
					for (int n = bin_start; n < bin_end; n++) // for each neighbor
					{
						int n_index = neighbor_bin_index[n];
						const Eigen::RowVector3f& pos_j = positions->row(n_index);

						if (n_index != i) // omit same particle
						{
							densities[i] += poly6Kernel(pos_i - pos_j); // equation (2)
							d_pk_Ci += gradSpikyKernel(pos_i - pos_j);
							d_pk_Ci_same = -1.0f * gradSpikyKernel(pos_i - pos_j);
							d_pk_Ci_same *= 1.0f / REST_DENSITY;
							denominator += d_pk_Ci_same.squaredNorm();

						}
					}
				}
			}

			d_pk_Ci *= 1.0f / REST_DENSITY;
			denominator += d_pk_Ci.squaredNorm();
			float C_i = densities[i] / REST_DENSITY - 1.f;

			denominator += EPSILON;
			lambdas[i] = - C_i / denominator; // equation (11)


			// calculate delta_p_i, equation (12, rsp. 14 for tensile stability)
			Eigen::Vector3f delta_pi = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
			for (int z = -1; z < 2; z++)
			{
				for (int y = -1; y < 2; y++)
				{
					int actual_bin = p_bin - 1 + y * gridWidth + z * gridWidth * gridWidth;
					int bin_start = bin_prefix_sum[actual_bin];
					int bin_end = bin_prefix_sum[actual_bin + 1 + 2];
					for (int n = bin_start; n < bin_end; n++) // for each neighbor
					{
						int n_index = neighbor_bin_index[n];
						const Eigen::RowVector3f& pos_j = positions->row(n_index);

						if (n_index != i) // omit same particle
						{
							delta_pi += (lambdas[i] + lambdas[n_index]) * gradSpikyKernel(pos_i - pos_j); // equation (12)
						}
					}
				}
			}
			delta_pi *= 1 / REST_DENSITY;


			Eigen::RowVector3f contactPoint;
			Eigen::Vector3f normal; // unit surface normal

			if (collision(pos_i, contactPoint, normal)) {
				//std::cout << "collision" << std::endl;
				const Eigen::Vector3f velocity = velocities.row(i);
				float penetrationDepth = (pos_i - contactPoint).norm();
				velocities.row(i) = velocity - normal * (1 + 0.5f * penetrationDepth / (m_dt * velocity.norm())) * velocity.dot(normal);
				positionsStar->row(i) = contactPoint;
			}

			// update position
			positionsStar->row(i) += delta_pi;
		}
	}

	for (int i = 0; i < n; ++i) {
		// update velocity
		velocities.row(i) = 1 / m_dt * (positionsStar->row(i) - positions->row(i));

		// OPTIONAL: apply vorticity confinement & XSPH viscosity
		// not for now
	}

	// advance step
	m_step++;
	m_time += m_dt;
	return false;
}
*/

void FluidSim::renderRenderGeometry(
	igl::opengl::glfw::Viewer& viewer) {
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

	// only render the first NUM_FLUID_PARTICLES rows from positions. The rows left are boundary particles.
	//p_iviewer->updatePositions(&Eigen::Matrix<float, -1, -1, Eigen::RowMajor>(positions->topRows(NUM_FLUID_PARTICLES)));
	p_iviewer->updatePositions(positions);
	//p_iviewer->updateColors(renderColors);
	p_iviewer->drawInstanced();
}

// density estimation
float FluidSim::poly6Kernel(Eigen::Vector3f r) {
	float r_norm = r.norm();
	if (r_norm <= NEIGHBOURHOOD_RADIUS) {
		return pow_h_9 * std::pow(NEIGHBOURHOOD_RADIUS * NEIGHBOURHOOD_RADIUS - r_norm * r_norm, 3);
	}
	else {
		return 0;
	}
}

// pressure gradient
float FluidSim::spikyKernel(Eigen::Vector3f r) {
	float r_norm = r.norm();
	if (r_norm <= NEIGHBOURHOOD_RADIUS) {
		return pow_h_6_1 * std::pow(NEIGHBOURHOOD_RADIUS - r_norm, 3);
	}
	else {
		return 0;
	}
}

// spiky gradient
Eigen::Vector3f FluidSim::gradSpikyKernel(Eigen::Vector3f r) {
	float r_norm = r.norm();
	if (0.0001 <= r_norm && r_norm <= NEIGHBOURHOOD_RADIUS) {
		return pow_h_6_2 * std::pow(NEIGHBOURHOOD_RADIUS - r_norm, 2) * r / r_norm;
	}
	else {
		return Eigen::Vector3f(0.0f, 0.0f, 0.0f);
	}
}

// collision detection
bool FluidSim::collision(const Eigen::Vector3f pos, Eigen::RowVector3f& contactPoint, Eigen::Vector3f& normal) {
	bool collided = false;
	normal = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
	contactPoint = Eigen::Vector3f(pos.x(), pos.y(), pos.z());
	if (pos.x() < simBoundary[0])
	{
		normal.x() += 1.0f;
		contactPoint.x() = simBoundary[0];
		collided = true;
	}
	if (pos.x() > simBoundary[1])
	{
		normal.x() += -1.0f;
		contactPoint.x() = simBoundary[1];
		collided = true;
	}
	if (pos.y() < simBoundary[2])
	{
		normal.y() += 1.0f;
		contactPoint.y() = simBoundary[2];
		collided = true;
	}
	if (pos.y() > simBoundary[3])
	{
		normal.y() += -1.0f;
		contactPoint.y() = simBoundary[3];
		collided = true;
	}
	if (pos.z() < simBoundary[4])
	{
		normal.z() += 1.0f;
		contactPoint.z() = simBoundary[4];
		collided = true;
	}
	if (pos.z() > simBoundary[5])
	{
		normal.z() += -1.0f;
		contactPoint.z() = simBoundary[5];
		collided = true;
	}
	return collided;
}
