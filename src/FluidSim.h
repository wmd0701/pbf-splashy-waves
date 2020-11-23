#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "Simulation.h"
#include "InstancedViewer.h"
#include <vector>
#include <atomic>
#include <condition_variable>

#define TIME_STEP_SIZE 0.05f
#define REST_DENSITY 1.f
#define NEIGHBOURHOOD_RADIUS 2.0f
// user specified relaxation parameter (equation 11):
// Bigger -> slower contraction of particles: influences the negative pressure inverse proportionally
#define EPSILON 1.0f

// neighbourhood_radius and particle_radius must be the same.
// else the visualization and grid generation generates particles too close
#define PARTICLES_PER_CUBE_SIDE 12
#define PARTICLE_RADIUS 0.5f

#define halfBoundarySize 2.0f*PARTICLES_PER_CUBE_SIDE*PARTICLE_RADIUS

#define ThreadCount 1

// Reusable Barrier class to synchronize all threads. C++20 has its own std::barrier but it doesn't work for me.
class Barrier
{

public:
	Barrier(int count)
		: thread_count(count)
		, counter(0)
		, waiting(0)
	{}

	void wait()
	{
		//fence mechanism
		std::unique_lock<std::mutex> lk(m);
		++counter;
		++waiting;
		cv.wait(lk, [&] {return counter >= thread_count; });
		cv.notify_all();
		--waiting;
		if (waiting == 0)
		{
			//reset barrier
			counter = 0;
		}
		lk.unlock();
	}

private:
	std::mutex m;
	std::condition_variable cv;
	int counter;
	int waiting;
	int thread_count;
};

class FluidSim : public Simulation {
public:
	FluidSim();

	virtual void init() override;
	virtual void resetMembers() override;
	virtual void updateRenderGeometry() override;
	virtual bool advance() override;
	virtual void renderRenderGeometry(igl::opengl::glfw::Viewer& viewer) override;

	// simulation boundary AABB [minx,miny,minz, maxx,maxy,maxz] (y is up)
	Eigen::VectorXf simBoundary;
	InstancedViewer* p_iviewer;

	//NeighborhoodSearch member variables. A loot and very badly named
	const int envWidth = 2.0f * halfBoundarySize + NEIGHBOURHOOD_RADIUS * 6.0f; // TODO: border of 3! cells that are "useless" reduce to 1 when properly computing collision!
	const int gridWidth = std::ceil(envWidth / NEIGHBOURHOOD_RADIUS);
	const float relPos = gridWidth / (float)envWidth;
	const float posOffsety = NEIGHBOURHOOD_RADIUS * 3.0f;
	const float posOffsetxz = envWidth / 2.0f;
	std::vector<unsigned int> bin_index;
	std::vector<unsigned int> bin_sub_index;
	//std::vector<unsigned int> bin_count;
	std::vector<std::atomic<unsigned int>> bin_count;
	std::vector<unsigned int> bin_prefix_sum;
	std::vector<unsigned int> neighbor_bin_index;


	std::vector<std::thread> threads;
	std::vector<int> indices;

	// Precomputed Powers
	const float pow_h_9 = 315.f / (64.f * M_PI * std::pow(NEIGHBOURHOOD_RADIUS, 9));
	const float pow_h_6_1 = 15.f / (M_PI * std::pow(NEIGHBOURHOOD_RADIUS, 6));
	const float pow_h_6_2 = -45.f / (M_PI * std::pow(NEIGHBOURHOOD_RADIUS, 6));

	float m_dt;

	/* PARTICLE DATA
	 * positions (and colors) are double buffered,
	 * as one (positions) needs to be available for the gpu.
	 * Therefore only write to positionsStar
	 * Reading from positions is safe and stores current position.
	 * Velocities, densities and lambdas can be modified at will.
	 * Each particle entry is stored in a single row of each matrix.
	 */
	bool initializedInstancedViewer = false;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positionsStar;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positions;

	Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* velocitiesStar;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* velocities;
	Eigen::VectorXf densities;
	Eigen::VectorXf lambdas;

	/* vectors to store individual colors for each particle.
	 * requires setPerinstanceColor(true) in instancedViewer.
	 * only write to updateColors, renderColors are read by gpu.
	 */
	Eigen::VectorXf* updateColors;
	Eigen::VectorXf* renderColors;

	// Actual storage Vectors for positions and colors. Do not modify.
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> positions1;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> positions2;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> velocities1;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> velocities2;
	Eigen::VectorXf colors1;
	Eigen::VectorXf colors2;


	// slow igl renderer data. Used for boundaries/obstacles and floor.
	Eigen::MatrixXd m_renderV;
	Eigen::MatrixXi m_renderF;
	Eigen::MatrixXd m_renderC;
	Eigen::MatrixXd m_renderUV;

	float poly6Kernel(Eigen::Vector3f distance);
	float spikyKernel(Eigen::Vector3f distance);
	Eigen::Vector3f gradSpikyKernel(Eigen::Vector3f r);
	bool collision(const Eigen::Vector3f pos, Eigen::RowVector3f& contactPoint, Eigen::Vector3f& normal);
};
#endif