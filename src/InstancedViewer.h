#ifndef INSTANCEDVIEWER_H
#define INSTANCEDVIEWER_H

#include <string>
#include "igl/opengl/glfw/Viewer.h"

class InstancedViewer
{
public:
	// Initialize a instanced viewer with a (x,3) matrix storing positions of the particles to draw.
	// and optionally a x-valued vector storing a color value (as float). x must be the same
	// the given pointers must remain valid and the values may not change during rendering!
	InstancedViewer(Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positions, Eigen::VectorXf* colors = NULL);
	~InstancedViewer() {
		glDeleteShader(m_vertexShader);
		glDeleteShader(m_fragmentShader);
	};

	// Set up shaders and basic sphere shaped mesh for particles.
	void init(igl::opengl::glfw::Viewer* v);

	// Draws all particles in a single draw call.
	void drawInstanced();

	// Change the positions of the particles.
	void updatePositions(Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positions, unsigned int toprows);

	// coloring can be uniform for all particles or per Instance.
	void updateColors(Eigen::VectorXf* colors);

	// change the amount of particles to be drawn with this function.
	void changeParticleCount(Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positions, Eigen::VectorXf* colors = NULL);

	// sets the radius for all particles
	bool setPerInstanceColor(bool isPerInstance);
	void setUniformColor(Eigen::Vector3f color);
	void setParticleSize(float radius);

private:
	void initShaders();
	void initParticleMesh();
	
	bool perInstanceColor = false;
	
	// reference to the igl viewer created by the GUI class.
	// it is used as a hook into the "normal rendering" for retrieving projection matrices and such.
	igl::opengl::glfw::Viewer* p_iglViewer;

	// Vertices, Normals and Index data for a single particle. Populated in initParticleMesh.
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> m_particleMeshData;
	Eigen::Matrix<unsigned short, -1, -1, Eigen::RowMajor> m_particleMeshIndices;

	// Pointer to the vector containing a single float color gradient value for each particle.
	// Is required to have as many entries as there are particles in the simulation.
	Eigen::VectorXf* p_particleColors;

	// Pointer to the (x,3) matrix containing 3 float values per position of each particle.
	// Its size defines how many particles are drawn and at which position.
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* p_particlePositions;
	unsigned int m_particleCount = 0;

	// Fallback matrices for when the pointers do not comply
	Eigen::VectorXf dummyColor;
	Eigen::Matrix<float, -1, -1, Eigen::RowMajor> dummyPosition;

	// Indices to the opengl buffers used by the shaders.
	unsigned int m_VAO;

	unsigned int m_colorPerInstance;
	unsigned int m_uniformParticleColor;
	unsigned int m_particleSize;

	// these are exactly the ones as in the igl viewer, but copied. 
	// (therefore inefficient as they should already be somewhere on the GPU anyways)
	unsigned int m_viewMatrix;
	unsigned int m_projMatrix;

	// these are the instanced buffer objects.
	unsigned int m_positionsVBO;
	unsigned int m_colorsVBO;

	unsigned int m_shaderProgram;
	unsigned int m_vertexShader;
	unsigned int m_fragmentShader;
};

#endif
