#include "DummySim.h"

using namespace std;
/*
 * Example simulation that changes the colors of a cube.
 */
DummySim::DummySim() : Simulation() {
	init();
}

void DummySim::init() {
	// just draw a point to satisfy igl... with no geometry it seems to crash!
	// vertices
	m_renderV.resize(3, 3);
	m_renderV << 0.f, 0.f, 0.f, 0.1f, 0.0f, 0.0f, 0.2f, 0.f, 0.f;
	// faces
	m_renderF.resize(1, 3);
	m_renderF << 0, 1, 2;
	// face colors
	m_renderC.resize(3, 3);
	m_renderC << 0.0f, 0.0f, 0.0f;

	int n = 100;
	positions1.resize(n * n * n, 3);
	positions2.resize(n * n * n, 3);
	updatePositions = &positions1;
	renderPositions = &positions2;

	colors1.resize(n * n * n);
	colors2.resize(n * n * n);
	updateColors = &colors1;
	renderColors = &colors2;

	float x = -n/2.0f;
	for (int i = 0; i < n; i++)
	{
		x += 1.0f;
		float y = -n/2.0f;
		for (int j = 0; j < n; j++)
		{
			y += 1.0f;
			float z = -n/2.0f;
			for (int k = 0; k < n; k++)
			{
				z += 1.0f;
				positions1.block(i * n * n + j * n + k, 0, 1, 3) << x, y, z;
				colors1[i * n * n + j * n + k] = (float)(i * n * n + j * n) / (float)(n * n * n);
			}
		}
	}

	colors2 = colors1;
	positions2 = positions1;

	initializedInstancedViewer = false;
	p_iviewer = new InstancedViewer(renderPositions, renderColors);
}

void DummySim::resetMembers() {
	int n = 100;
	positions1.resize(n * n * n, 3);
	positions2.resize(n * n * n, 3);
	updatePositions = &positions1;
	renderPositions = &positions2;

	colors1.resize(n * n * n);
	colors2.resize(n * n * n);
	updateColors = &colors1;
	renderColors = &colors2;

	float x = -n / 2.0f;
	for (int i = 0; i < n; i++)
	{
		x += 1.0f;
		float y = -n / 2.0f;
		for (int j = 0; j < n; j++)
		{
			y += 1.0f;
			float z = -n / 2.0f;
			for (int k = 0; k < n; k++)
			{
				z += 1.0f;
				positions1.block(i * n * n + j * n + k, 0, 1, 3) << x, y, z;
				colors1[i * n * n + j * n + k] = (float)(i * n * n + j * n) / (float)(n * n * n);
			}
		}
	}
}

void DummySim::updateRenderGeometry() {
	// just swap the pointers, no copying!
	auto temp = renderPositions;
	auto temp2 = renderColors;
	renderPositions = updatePositions;
	updatePositions = temp;
	renderColors = updateColors;
	updateColors = temp2;
}

bool DummySim::advance() {
	// do next step of some color animation
	int speed = 60;
	int decColor = (m_step / speed) % 3;
	int incColor = (decColor + 1) % 3;

	
	for (int i = 0; i < updatePositions->rows(); i++)
	{
		(*updatePositions)(i, 0) = (*renderPositions)(i,0) + 0.1f; // example update
	}
	
	m_step++;
	return false;
}

void DummySim::renderRenderGeometry(
	igl::opengl::glfw::Viewer &viewer) {
	//prevent the basic renderer from crashing by adding just anything.
	viewer.data().set_mesh(m_renderV, m_renderF);
	viewer.data().set_colors(m_renderC);

	if (!initializedInstancedViewer) { // do some things once!
		//viewer.core.background_color = Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
		viewer.core.animation_max_fps = 60;
		p_iviewer->init(&viewer); 
		//p_iviewer->setPerInstanceColor(false);
		initializedInstancedViewer = true; 
	}

	p_iviewer->updatePositions(renderPositions, 1);
	p_iviewer->drawInstanced();
}