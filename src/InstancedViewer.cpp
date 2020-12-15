#include "InstancedViewer.h"

// parameter toprows: only render the first toprows particles. The particles left are boundary(ghost) particles
InstancedViewer::InstancedViewer(Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positions, Eigen::VectorXf* colors, unsigned int toprows)
{
	dummyColor.resize(1);
	dummyColor[0] = 0.5f;

	dummyPosition.resize(1, 3);
	dummyPosition << 0.0f, 0.0f, 0.0f;

	if (positions->cols() != 3) 
	{
		std::cout << "Instanced Viewer requires positional data in a (3,x) matrix. Using dummy particle.\n";
		perInstanceColor = false;
		m_particleCount = 1;
		p_particlePositions = &dummyPosition;
		p_particleColors = &dummyColor;
	}
	else 
	{
		m_particleCount = toprows == 0? positions->rows() : toprows;
		p_particlePositions = positions;

		if (colors != nullptr)
		{
			if (colors->size() >= m_particleCount)
			{
				perInstanceColor = true;
				p_particleColors = colors;
			}
			else
			{
				//std::cout << "Instanced Viewer requires exactly one float value as color per particle.\n";
				std::cout << "Less colors than needed. Setting needed values to dummy color.\n";
				perInstanceColor = false;
				dummyColor = Eigen::VectorXf::Ones(m_particleCount) * 0.0f; // dummy color of 0
				dummyColor.head(colors->size()) = *colors;
				p_particleColors = &dummyColor;
			}
		}
		else 
		{
			perInstanceColor = false;
			dummyColor = Eigen::VectorXf::Zero(m_particleCount);
			p_particleColors = &dummyColor;
		}
	}

}

void InstancedViewer::init(igl::opengl::glfw::Viewer* v)
{
	p_iglViewer = v;

	initParticleMesh();
	
	initShaders();
	
	glUseProgram(m_shaderProgram);
	
	glGenVertexArrays(1, &m_VAO);
	glBindVertexArray(m_VAO);
	

	// Create and Bind Particle Mesh data Buffers (Vertices, Normals and Indices)
	unsigned int VBO;
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * m_particleMeshData.size(), m_particleMeshData.data(), GL_STATIC_DRAW);

	glEnableVertexAttribArray(0); // Vertices
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);

	glEnableVertexAttribArray(1); // Normals
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3*sizeof(float)));

	unsigned int EBO;
	glGenBuffers(1, &EBO); // Indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * m_particleMeshIndices.size(), m_particleMeshIndices.data(), GL_STATIC_DRAW);


	// Create Instanced Buffers for position and color, only one float per color!
	glGenBuffers(1, &m_positionsVBO);
	glBindBuffer(GL_ARRAY_BUFFER, m_positionsVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * p_particlePositions->size(), p_particlePositions->data(), GL_DYNAMIC_DRAW);

	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glVertexAttribDivisor(2, 1);

	glGenBuffers(1, &m_colorsVBO);
	glBindBuffer(GL_ARRAY_BUFFER, m_colorsVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * p_particleColors->size(), p_particleColors->data(), GL_DYNAMIC_DRAW);

	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
	glVertexAttribDivisor(3, 1);


	m_uniformParticleColor = glGetUniformLocation(m_shaderProgram, "particleColor");
	glUniform3f(m_uniformParticleColor, 0.0f, 0.0f, 1.0f); // basic is bright blue
	m_particleSize = glGetUniformLocation(m_shaderProgram, "particleSize");
	glUniform1f(m_particleSize, 0.5f); // default is radius of 0.5
	m_colorPerInstance = glGetUniformLocation(m_shaderProgram, "colorPerInstance");
	if (perInstanceColor) glUniform1f(m_colorPerInstance, 1.0f);
	else glUniform1f(m_colorPerInstance, 0.0f);

	m_viewMatrix = glGetUniformLocation(m_shaderProgram, "view");
	glUniformMatrix4fv(m_viewMatrix, 1, GL_FALSE, p_iglViewer->core.view.data());
	m_projMatrix = glGetUniformLocation(m_shaderProgram, "proj");
	glUniformMatrix4fv(m_projMatrix, 1, GL_FALSE, p_iglViewer->core.proj.data());

	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void InstancedViewer::initShaders()
{
	std::string vshader_source =
		R"(
			#version 330 core
			layout (location = 0) in vec3 vertexPos;
			layout (location = 1) in vec3 vertexNormal;
			layout (location = 2) in vec3 particlePosition;
			layout (location = 3) in float InstanceColor;
			
			uniform float particleSize;
			uniform mat4 view;
			uniform mat4 proj;

			out vec3 position_world;
			out vec3 normal_world;
			out vec3 color;

			void main()
			{
				normal_world = vertexNormal;
				position_world = vertexPos;
				// color gradient from blue to white
				color = vec3(1.0,1.0,1.0) * InstanceColor + vec3(0.0,0.0,1.0) * (1.0-InstanceColor);
				gl_Position = proj * view * vec4(particleSize * vertexPos + particlePosition, 1.0);
			}
		)";

	std::string fshader_source =
		R"(
			#version 330 core
			out vec4 FragColor;

			in vec3 normal_world;
			in vec3 position_world;
			in vec3 color;
			
			uniform vec3 particleColor;
			uniform float colorPerInstance;

			void main()
			{
				vec3 lightPosition = vec3(100.0,100.0,100.0);

				vec3 lightDirection = normalize(lightPosition - position_world);  
				float cosineAngle = max(dot(normal_world, lightDirection), 0.0);
				vec3 diffuse = cosineAngle * vec3(1.0, 1.0, 0.95); //slightly yellow light
				
				vec3 ambient = vec3(0.3, 0.3, 0.3);

				// choose wheter to use uniform particle color or instanced color
				vec3 actualColor = colorPerInstance * color + (1.0-colorPerInstance) * particleColor;
				//vec3 actualColor = particleColor;

				FragColor = vec4((ambient + diffuse) * actualColor, 1.0);
			} 
		)";

	m_vertexShader = glCreateShader(GL_VERTEX_SHADER);
	const char* vs = vshader_source.c_str();
	glShaderSource(m_vertexShader, 1, &vs, NULL);
	glCompileShader(m_vertexShader);

	int  success;
	char infoLog[512];
	glGetShaderiv(m_vertexShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(m_vertexShader, 512, NULL, infoLog);
		std::cout << "ERROR: Vertex shader compilation failed\n" << infoLog << std::endl;
	}

	m_fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	const char* fs = fshader_source.c_str();
	glShaderSource(m_fragmentShader, 1, &fs, NULL);
	glCompileShader(m_fragmentShader);

	glGetShaderiv(m_fragmentShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(m_fragmentShader, 512, NULL, infoLog);
		std::cout << "ERROR: Fragment shader compilation failed\n" << infoLog << std::endl;
	}


	m_shaderProgram = glCreateProgram();
	glAttachShader(m_shaderProgram, m_vertexShader);
	glAttachShader(m_shaderProgram, m_fragmentShader);
	glLinkProgram(m_shaderProgram);

	glGetProgramiv(m_shaderProgram, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(m_shaderProgram, 512, NULL, infoLog);
		std::cout << "ERROR: shader program compilation failed\n" << infoLog << std::endl;
	}
}

void InstancedViewer::drawInstanced()
{
	//glEnable(GL_DEPTH_TEST);
	glUseProgram(m_shaderProgram);
	glBindVertexArray(m_VAO);
	
	glUniformMatrix4fv(m_viewMatrix, 1, GL_FALSE, p_iglViewer->core.view.data());
	glUniformMatrix4fv(m_projMatrix, 1, GL_FALSE, p_iglViewer->core.proj.data());

	glDrawElementsInstanced(GL_TRIANGLES, m_particleMeshIndices.size(), GL_UNSIGNED_SHORT, 0, m_particleCount);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// keeping the number of particles the same, just update their positions only.
void InstancedViewer::updatePositions(Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positions, unsigned int toprows)
{
	glBindBuffer(GL_ARRAY_BUFFER, m_positionsVBO);
	//glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * positions->size(), positions->data());
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * positions->topRows(toprows).size(), positions->topRows(toprows).data());
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// keeping the number of particles the same, just update their positions only.
void InstancedViewer::updatePositions(Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positions)
{
	glBindBuffer(GL_ARRAY_BUFFER, m_positionsVBO);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * positions->size(), positions->data());
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// keeping the particle count the same, just update their color only.
void InstancedViewer::updateColors(Eigen::VectorXf* colors)
{
	glBindBuffer(GL_ARRAY_BUFFER, m_colorsVBO);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * colors->size(), colors->data());
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// change the number of particles may only be called after init() was called.
void InstancedViewer::changeParticleCount(Eigen::Matrix<float, -1, -1, Eigen::RowMajor>* positions, Eigen::VectorXf* colors)
{
	if (positions->cols() != 3)
	{
		std::cout << "new positional data is required in a (3,x) matrix. Using dummy particle.\n";
		perInstanceColor = false;
		m_particleCount = 1;
		p_particlePositions = &dummyPosition;
		p_particleColors = &dummyColor;
	}
	else
	{
		m_particleCount = positions->rows();
		p_particlePositions = positions;

		if (colors != nullptr)
		{
			if (colors->size() == m_particleCount)
			{
				perInstanceColor = true;
				p_particleColors = colors;
			}
			else
			{
				std::cout << "Instanced Viewer requires exactly one float value as color per particle.\n";
				perInstanceColor = false;
				dummyColor = Eigen::VectorXf::Zero(m_particleCount);
				p_particleColors = &dummyColor;
			}
		}
		else
		{
			perInstanceColor = false;
			dummyColor = Eigen::VectorXf::Zero(m_particleCount);
			p_particleColors = &dummyColor;
		}
	}

	glBindBuffer(GL_ARRAY_BUFFER, m_positionsVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * positions->size(), positions->data(), GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, m_colorsVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * colors->size(), colors->data(), GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}


// set color drawing per instance. 
// Requires the p_particleColors vector to be set up correctly when activating.
// returns wheter the coloring is per instance from now on
bool InstancedViewer::setPerInstanceColor(bool isPerInstance)
{
	if (perInstanceColor) // change to uniform
	{
		if (!isPerInstance)
		{
			perInstanceColor = false;
			glUniform1f(m_colorPerInstance, 0.0f);
			return false;
		}
		return true;
	}
	else // change to instanced if requested.
	{
		if (isPerInstance)
		{
			// only make switch if there are enough colors per particle!
			if (p_particlePositions->rows() >= p_particleColors->size())
			{
				perInstanceColor = true;
				glUniform1f(m_colorPerInstance, 1.0f);
				return true;
			}
			else std::cout << "CAUTION: cannot change to instanced colors, array size mismatch.\n";
			return false;
		}
	}
}

void InstancedViewer::setUniformColor(Eigen::Vector3f color)
{
	glUniform3f(m_uniformParticleColor, color[0], color[1], color[2]);
}

void InstancedViewer::setParticleSize(float radius=0.5f)
{
	glUniform1f(m_particleSize, radius);
}

void InstancedViewer::initParticleMesh()
{
	m_particleMeshData.resize(30, 6); // position and normal
	m_particleMeshData << 0.000000f, 0.809017f, -0.587785f, 0.0000f, -1.0000f, 0.0000f,
		0.000000f, 0.309017f, -0.951057f, 0.0000f, -0.7689f, -0.6393f,
		0.000000f, -0.309017f, -0.951056f, 0.4998f, -0.7689f, -0.3986f,
		0.000000f, -0.809017f, -0.587785f, 0.0000f, 0.2864f, -0.9581f,
		0.459549f, 0.809017f, -0.366478f, 0.7490f, -0.2864f, -0.5973f,
		0.743566f, 0.309017f, -0.592974f, 0.0000f, -0.2864f, -0.9581f,
		0.743566f, -0.309017f, -0.592974f, 0.0000f, 0.7689f, -0.6393f,
		0.459549f, -0.809017f, -0.366478f, 0.0000f, 1.0000f, 0.0000f,
		0.573048f, 0.809017f, 0.130795f, 0.4998f, 0.7689f, -0.3986f,
		0.927212f, 0.309017f, 0.211630f, 0.7490f, 0.2864f, -0.5973f,
		0.927212f, -0.309017f, 0.211630f, 0.9340f, 0.2864f, 0.2132f,
		0.573048f, -0.809017f, 0.130795f, 0.6233f, -0.7689f, 0.1422f,
		0.255030f, 0.809017f, 0.529576f, 0.9340f, -0.2864f, 0.2132f,
		0.412648f, 0.309017f, 0.856873f, 0.6233f, 0.7689f, 0.1422f,
		0.412648f, -0.309017f, 0.856872f, 0.4157f, 0.2864f, 0.8632f,
		0.255030f, -0.809017f, 0.529576f, 0.2774f, -0.7689f, 0.5760f,
		-0.000000f, 1.000000f, 0.000000f, 0.4157f, -0.2864f, 0.8632f,
		-0.255031f, 0.809017f, 0.529576f, 0.2774f, 0.7689f, 0.5760f,
		-0.412648f, 0.309017f, 0.856873f, -0.4157f, 0.2864f, 0.8632f,
		-0.412648f, -0.309017f, 0.856872f, -0.2774f, -0.7689f, 0.5760f,
		-0.255031f, -0.809017f, 0.529576f, -0.4157f, -0.2864f, 0.8632f,
		-0.573048f, 0.809017f, 0.130795f, -0.2774f, 0.7689f, 0.5760f,
		-0.927212f, 0.309017f, 0.211630f, -0.9340f, 0.2864f, 0.2132f,
		-0.927212f, -0.309017f, 0.211630f, -0.6233f, -0.7689f, 0.1422f,
		-0.573048f, -0.809017f, 0.130795f, -0.9340f, -0.2864f, 0.2132f,
		-0.459549f, 0.809017f, -0.366478f, -0.6233f, 0.7689f, 0.1422f,
		-0.743566f, 0.309017f, -0.592974f, -0.7490f, 0.2864f, -0.5973f,
		-0.743566f, -0.309017f, -0.592974f, -0.4998f, -0.7689f, -0.3986f,
		-0.459549f, -0.809017f, -0.366478f, -0.7490f, -0.2864f, -0.5973f,
		0.000000f, -1.000000f, 0.000000f, -0.4998f, 0.7689f, -0.3986f;

	m_particleMeshIndices.resize(56, 3);
	m_particleMeshIndices << 29, 3, 7,
		1, 6, 2,
		0, 16, 4,
		2, 7, 3,
		0, 5, 1,
		4, 9, 5,
		29, 7, 11,
		5, 10, 6,
		4, 16, 8,
		6, 11, 7,
		8, 13, 9,
		29, 11, 15,
		9, 14, 10,
		8, 16, 12,
		10, 15, 11,
		12, 18, 13,
		29, 15, 20,
		13, 19, 14,
		12, 16, 17,
		14, 20, 15,
		17, 22, 18,
		29, 20, 24,
		18, 23, 19,
		17, 16, 21,
		19, 24, 20,
		21, 26, 22,
		29, 24, 28,
		23, 26, 27,
		21, 16, 25,
		23, 28, 24,
		25, 1, 26,
		29, 28, 3,
		27, 1, 2,
		25, 16, 0,
		27, 3, 28,
		1, 5, 6,
		2, 6, 7,
		0, 4, 5,
		4, 8, 9,
		5, 9, 10,
		6, 10, 11,
		8, 12, 13,
		9, 13, 14,
		10, 14, 15,
		12, 17, 18,
		13, 18, 19,
		14, 19, 20,
		17, 21, 22,
		18, 22, 23,
		19, 23, 24,
		21, 25, 26,
		23, 22, 26,
		23, 27, 28,
		25, 0, 1,
		27, 26, 1,
		27, 2, 3;
}
