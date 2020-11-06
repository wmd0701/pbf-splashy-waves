#include "Simulation.h"
#include "InstancedViewer.h"

/*
 * Example simulation that changes the colors of a cube.
 */
class DummySim : public Simulation {
public:
	DummySim();
	~DummySim() {}

	virtual void init() override;
	virtual void resetMembers() override;
	virtual void updateRenderGeometry() override;
	virtual bool advance() override;
	virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) override;

	InstancedViewer* p_iviewer;
private:

	bool initializedInstancedViewer = false;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* updatePositions;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* renderPositions;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> positions1;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> positions2;
	Eigen::VectorXf* updateColors;
	Eigen::VectorXf* renderColors;
	Eigen::VectorXf colors1;
	Eigen::VectorXf colors2;

	// only for static geometry!
	Eigen::MatrixXd m_renderV;  // vertex positions for rendering
	Eigen::MatrixXi m_renderF;  // face indices for rendering
	Eigen::MatrixXd m_renderC;  // colors per face for rendering
};