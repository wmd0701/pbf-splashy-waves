#include "Simulation.h"
#include "InstancedViewer.h"
#include <vector>

#define TIME_STEP_SIZE 0.1f
#define REST_DENSITY 1.f
#define NEIGHBOURHOOD_RADIUS 5.f
// user specified relaxation parameter (equation 11):
// Bigger -> slower contraction of particles: influences the negative pressure inverse proportionally
#define EPSILON 1.0f

// neighbourhood_radius and particle_radius must be the same.
// else the visualization and grid generation generates particles too close
#define PARTICLES_PER_CUBE_SIDE 7
#define PARTICLE_RADIUS 0.5f


class FluidSim : public Simulation {
public:
	FluidSim();

	virtual void init() override;
	virtual void resetMembers() override;
	virtual void updateRenderGeometry() override;
	virtual bool advance() override;
	virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) override;

	InstancedViewer* p_iviewer;
private:
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
	
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> velocities;
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
	Eigen::VectorXf colors1;
	Eigen::VectorXf colors2;
	

	// slow igl renderer data. Used for boundaries/obstacles and floor.
	Eigen::MatrixXd m_renderV;
	Eigen::MatrixXi m_renderF;
	Eigen::MatrixXd m_renderC;
	Eigen::MatrixXd m_renderUV;

	float poly6Kernel(Eigen::Vector3f distance, float h);
	float spikyKernel(Eigen::Vector3f distance, float h);
	Eigen::Vector3f gradSpikyKernel(Eigen::Vector3f r, float h);
	bool collision(const Eigen::Vector3f pos, Eigen::RowVector3f &contactPoint, Eigen::Vector3f &normal);
};
