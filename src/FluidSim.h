#include "Simulation.h"
#include <vector>

#define TIME_STEP_SIZE 0.1d
#define REST_DENSITY 0.1d
#define EPSILON 0.0000000001d
#define NEIGHBOURHOOD_RADIUS 2.0d

using Eigen::Vector3d;

class Particle {
public:
	Particle(Vector3d);

	// intrinsic properties
	Vector3d m_position;
	Vector3d m_delta_position;
	Vector3d m_position_star;
	double m_lambda;
	Vector3d m_velocity;
	Vector3d m_acceleration;
	Vector3d m_force; // eg. due to collision

	double m_density;
	double m_constr_gradient;
};


/*
 * Example simulation that changes the colors of a cube.
 */
class FluidSim : public Simulation {
public:
	FluidSim();

	virtual void init() override;
	virtual void resetMembers() override;
	virtual void updateRenderGeometry() override;
	virtual bool advance() override;
	virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) override;

private:
	Eigen::MatrixXd m_V;  // vertex positions
	Eigen::MatrixXi m_F;  // face indices
	Eigen::MatrixXd m_C;  // colors per face

	Eigen::MatrixXd m_renderV;  // vertex positions for rendering
	Eigen::MatrixXi m_renderF;  // face indices for rendering
	Eigen::MatrixXd m_renderC;  // colors per face for rendering

	std::vector<Particle> m_particles; // individual particles for the simulation
	float m_particleRadius;

	double poly6Kernel(Vector3d distance, double h);
	double spikyKernel(Vector3d distance, double h);
	Vector3d gradSpikyKernel(Vector3d r, double h);
};
