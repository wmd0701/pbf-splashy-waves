#include "Simulation.h"
#include <vector>

class Particle {
public:
	Particle(Eigen::Vector3d);

	// intrinsic properties
	Eigen::Vector3d m_position;
	Eigen::Vector3d m_velocity;
	Eigen::Vector3d m_acceleration;
	Eigen::Vector3d m_force; // eg. due to collision

	float m_density;

	// external forces
	Eigen::Vector3d m_gravitationForce;

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
};
