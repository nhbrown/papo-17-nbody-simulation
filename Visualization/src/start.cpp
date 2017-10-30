#define _CRT_SECURE_NO_WARNINGS //Disable some warnings

//---Includes
// Std. includes
#include <iostream>
#include <map>
#include <string>
#include <windows.h>
#include <math.h>
#include <io.h>      /* Count files(iterations) MS-DOS , #include <unistd.h> for UNIX/LINUX */

// Glad
#include <glad/glad.h>

// GLFW
#include <GLFW/glfw3.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// FreeType
#include <ft2build.h>
#include FT_FREETYPE_H

// Shader class
#include <opengl/shader.h>

// Function Prototypes
#include "header.h" 
//---Includes---END

//---Settings
// Display settings
GLuint WIDTH = 800, HEIGHT = 600;
GLuint windowedWIDTH = WIDTH, windowedHEIGHT = HEIGHT; // Save width and height before entering fullscreen
const int antiAliasing = 2;
glm::vec3 textColor = glm::vec3(0.9, 0.9, 0.9); // Text color for data presentation
glm::vec3 scaleDistanceToCenterColor = glm::vec3(0, 1, 0);
glm::vec3 scaleDistanceOutsideColor = glm::vec3(1, 0, 0);

// Data settings
char dataFolder[128] = "run"; // Default folder for data
const int NUM_PARTICLE = countParticle(); // Set number of particles 
const int numOfIterations = countIterations(); // Set number of iterations
int particlesWithscaleDistanceToCenter = 0;
int particlesWithscaleDistanceOutside = 0;
int iterationCounter = 1; // Current iteration
int indexStar = 0; // Every iteration a new star
const float sparklingTime = 0.5;
int sparklingStar = 0;
float sinSparkling = 0;
double currentTime = 0;
double currentTimeForSparkling = 0;
double deltaTime = 0;
int waitOneFrame = 0;
int error = 0;
int simpleStars = 1;

int numberOfPoints = 93; // Number of points for star shape (number of vertices)
glm::vec3 *energy = (glm::vec3 *)malloc((numOfIterations+1) * sizeof(glm::vec3)); // Generate a new list of numOfIterations particles for energy data

//Hotkey customizable settings
int fullScreen = 0;
int active = 0;
int showDistance = 0;
int showData = 1; 
float pointSize = 0.3;
int slowMotion = 0;
int mouseLook = 0;
int numberOfMarkedParticles = 0;
int colorMode = 0;
float scaleDistanceToCenter = 1;
float scaleDistanceOutside = 5;
int funPaintMode = 0;

// Window input - detect one key press of hotkeys
static int oldStateSpace = GLFW_RELEASE;
static int oldStateS = GLFW_RELEASE;
static int oldStateM = GLFW_RELEASE;
static int oldStateD = GLFW_RELEASE;
static int oldStateF = GLFW_RELEASE;
static int oldStateC = GLFW_RELEASE;
static int oldStateP = GLFW_RELEASE;
static int oldStateAlt = GLFW_RELEASE;
static int oldStateEnter = GLFW_RELEASE;
static int oldStatePAGE_UP = GLFW_RELEASE;
static int oldStatePAGE_DOWN = GLFW_RELEASE;
static int oldStateKP_1 = GLFW_RELEASE;
static int oldStateKP_4 = GLFW_RELEASE;
static int oldStateKP_2 = GLFW_RELEASE;
static int oldStateKP_5 = GLFW_RELEASE;

// Mouse settings
bool firstMouse = true;
float yaw = -90.0f;	
float pitch = 0.0f;
float lastX = WIDTH / 2.0;
float lastY = HEIGHT / 2.0;
float fov = 45.0f;

// Camera settings
float alpha = 0.0f; // Default alpha angle
float zoom = 3.0f; // Default zoom

// Default camera position
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, zoom);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
//---Settings---END

// Holds all state information relevant to a character as loaded using FreeType
struct Character {
	GLuint TextureID;   // ID handle of the glyph texture
	glm::ivec2 Size;    // Size of glyph
	glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
	GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;
GLuint VAO, VBO;

// Set up vertex data of simple star
float defaultSimpleStar[] = { 0.0f, 0.0f, 0.0f };

// Set up vertex data of complex star
const float complexStar[] = {

	0.0000f,  0.0000f, 0.0000f,	0.0000f,  0.0001f, 0.0000f,	0.0000f,  -0.0001f, 0.0000f,	0.0001f,  0.0000f, 0.0000f,	-0.0001f,  0.0000f, 0.0000f,	0.0000f,  0.0000f, 0.0001f,
	0.0000f,  0.0000f, -0.0001f,	0.0000f,  0.0002f, 0.0000f,	0.0000f,  -0.0002f, 0.0000f,	0.0002f,  0.0000f, 0.0000f,	-0.0002f,  0.0000f, 0.0000f,	0.0000f,  0.0000f, 0.0002f,
	0.0000f,  0.0000f, -0.0002f,	0.0000f,  0.0003f, 0.0000f,	0.0000f,  -0.0003f, 0.0000f,	0.0003f,  0.0000f, 0.0000f,	-0.0003f,  0.0000f, 0.0000f,	0.0000f,  0.0000f, 0.0003f,
	0.0000f,  0.0000f, -0.0003f,	0.000000f,  -0.001f, 0.000000f,	0.000000f,  0.001f, 0.000000f,	0.001f,  0.000000f, 0.000000f,	-0.001f,  0.000000f, 0.000000f,	0.000000f,  0.000000f, 0.001f,
	0.000000f,  0.000000f, -0.001f,	0.000000f,  -0.0005f, 0.000000f,	0.000000f,  0.0005f, 0.000000f,	0.0005f,  0.000000f, 0.000000f,	-0.0005f,  0.000000f, 0.000000f,	0.000000f,  0.000000f, 0.0005f,
	0.000000f,  0.000000f, -0.0005f,	0.000000f,  -0.0001f, -0.0001f,	0.000000f,  0.0001f, -0.0001f,	0.000000f,  -0.0001f, 0.0001f,	0.000000f,  0.0001f, 0.0001f,	0.0001f,  0.000000f, -0.0001f,
	-0.0001f,  0.000000f, -0.0001f,	0.0001f,  0.000000f, 0.0001f,	-0.0001f,  0.000000f, 0.0001f,	0.0001f,  -0.0001f, 0.000000f,	-0.0001f,  -0.0001f, 0.000000f,	0.0001f,  0.0001f, 0.000000f,
	-0.0001f,  0.0001f, 0.000000f,	-0.0001f,  -0.0001f, -0.0001f,	-0.0001f,  -0.0001f, +0.0001f,	-0.0001f,  +0.0001f, +0.0001f,	+0.0001f,  +0.0001f, +0.0001f,	+0.0001f,  +0.0001f, -0.0001f,
	+0.0001f,  -0.0001f, -0.0001f,	-0.0001f,  +0.0001f, -0.0001f,	+0.0001f,  -0.0001f, +0.0001f,	0.000000f,  -0.0002f, -0.0002f,	0.000000f,  0.0002f, -0.0002f,	0.000000f,  -0.0002f, 0.0002f,
	0.000000f,  0.0002f, 0.0002f,	0.0002f,  0.000000f, -0.0002f,	-0.0002f,  0.000000f, -0.0002f,	0.0002f,  0.000000f, 0.0002f,	-0.0002f,  0.000000f, 0.0002f,	0.0002f,  -0.0002f, 0.000000f,
	-0.0002f,  -0.0002f, 0.000000f,	0.0002f,  0.0002f, 0.000000f,	-0.0002f,  0.0002f, 0.000000f,	-0.0002f,  -0.0002f, -0.0002f,	-0.0002f,  -0.0002f, +0.0002f,	-0.0002f,  +0.0002f, +0.0002f,
	+0.0002f,  +0.0002f, +0.0002f,	+0.0002f,  +0.0002f, -0.0002f,	+0.0002f,  -0.0002f, -0.0002f,	-0.0002f,  +0.0002f, -0.0002f,	+0.0002f,  -0.0002f, +0.0002f,	0.000000f,  -0.0005f, -0.0005f,
	0.000000f,  0.0005f, -0.0005f,	0.000000f,  -0.0005f, 0.0005f,	0.000000f,  0.0005f, 0.0005f,	0.0005f,  0.000000f, -0.0005f,	-0.0005f,  0.000000f, -0.0005f,	0.0005f,  0.000000f, 0.0005f,
	-0.0005f,  0.000000f, 0.0005f,	0.0005f,  -0.0005f, 0.000000f,	-0.0005f,  -0.0005f, 0.000000f,	0.0005f,  0.0005f, 0.000000f,	-0.0005f,  0.0005f, 0.000000f,	-0.0005f,  -0.0005f, -0.0005f,
	-0.0005f,  -0.0005f, +0.0005f,	-0.0005f,  +0.0005f, +0.0005f,	+0.0005f,  +0.0005f, +0.0005f,	+0.0005f,  +0.0005f, -0.0005f,	+0.0005f,  -0.0005f, -0.0005f,	-0.0005f,  +0.0005f, -0.0005f,
	+0.0005f,  -0.0005f, +0.0005f
};

//---Main
int main()
{
	// Glfw initialize and configure
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, antiAliasing); // Antialiasing

	// Glfw window creation
	GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "NBody_Visualization_2.0", nullptr, nullptr); // Windowed
	if (window == NULL)
	{
		printf("Failed to create GLFW window");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback); // Window resize
	glfwSetCursorPosCallback(window, mouse_callback);

	// Glad load all OpenGL function pointers
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		printf("Failed to initialize GLAD");
		return -1;
	}

	// Define the viewport dimensions
	glViewport(0, 0, WIDTH, HEIGHT);

	// Configure global opengl state
	glEnable(GL_MULTISAMPLE); //Enable antialiasing 
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_PROGRAM_POINT_SIZE);

	// Build and compile shaders for circle
	Shader shaderCircle("shaders/circle.vs", "shaders/circle.fs");
	
	// build and compile shaders for instancing
	Shader shaderInstance("shaders/instancing.vs", "shaders/instancing.fs");

	//---FreeType
	// Compile and setup the text shader
	Shader shaderText("shaders/text.vs", "shaders/text.fs");
	glm::mat4 projectionText = glm::ortho(0.0f, static_cast<GLfloat>(WIDTH), 0.0f, static_cast<GLfloat>(HEIGHT));
	shaderText.use();
	glUniformMatrix4fv(glGetUniformLocation(shaderText.ID, "projectionText"), 1, GL_FALSE, glm::value_ptr(projectionText));

	FT_Library ft;
	// All functions return a value different than 0 whenever an error occurred
	if (FT_Init_FreeType(&ft))
		printf("ERROR::FREETYPE: Could not init FreeType Library");

	// Load font as face
	FT_Face face;
	if (FT_New_Face(ft, "fonts/calibri.ttf", 0, &face))
		printf("ERROR::FREETYPE: Failed to load font");

	// Set size to load glyphs as
	FT_Set_Pixel_Sizes(face, 0, 48);

	// Disable byte-alignment restriction
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// Load first 128 characters of ASCII set
	for (GLubyte c = 0; c < 128; c++)
	{
		// Load character glyph 
		if (FT_Load_Char(face, c, FT_LOAD_RENDER))
		{
			printf("ERROR::FREETYTPE: Failed to load Glyph");
			continue;
		}
		// Generate texture
		GLuint texture;
		glGenTextures(1, &texture);
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexImage2D(
			GL_TEXTURE_2D,
			0,
			GL_RED,
			face->glyph->bitmap.width,
			face->glyph->bitmap.rows,
			0,
			GL_RED,
			GL_UNSIGNED_BYTE,
			face->glyph->bitmap.buffer
		);
		// Set texture options
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		// Store character for later use
		Character character = {
			texture,
			glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
			glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
			face->glyph->advance.x
		};
		Characters.insert(std::pair<GLchar, Character>(c, character));
	}
	glBindTexture(GL_TEXTURE_2D, 0);
	// Destroy FreeType once finished
	FT_Done_Face(face);
	FT_Done_FreeType(ft);

	// Configure VAO/VBO for texture quads
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	//---FreeType---END
	
	//Read energy-data from energy_diagnostics.csv
	readEnergyData(energy);

	// Render loop
	while (!glfwWindowShouldClose(window))
	{
		currentTime = glfwGetTime(); // For fps

		if((glfwGetTime() - currentTimeForSparkling) > sparklingTime) // If more than sparklingTime since last time
		{ 
			currentTimeForSparkling = glfwGetTime(); // Reset timer
			if (sparklingStar < NUM_PARTICLE)
			{
				sparklingStar++; // Mark next star for sparkling
			}
			else
			{
				sparklingStar = 0;
			}
			
		}
		else
		{
			sinSparkling = 1 / sparklingTime * (glfwGetTime() - currentTimeForSparkling) * 3.14; // SinSparkling goes from 0 - 3,14 => one star is smoothly sparkling 
		}
			
		// Hotkey input
		processInput(window);

		// render
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

		if (funPaintMode == 0 || waitOneFrame < 1)
		{
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		}
		if (funPaintMode == 1) // Otherwise data is still drawn
		{
			if (waitOneFrame == 0)
			{
				waitOneFrame = 1; 
			}	
		}

		if(showData == 1 && funPaintMode == 0)
		{ 
			RenderText(shaderText, "FPS: " + std::to_string((int)(1 / deltaTime)), 5.0f, 600 - 10.0f, 0.2f, textColor);

			if(active == 0) // Show only when active
			{ 				
				if (numberOfMarkedParticles > 0) // Show only if particles are marked
				{
					RenderText(shaderText, "Marked: " + std::to_string(numberOfMarkedParticles), 5.0f, 65.0f, 0.2f, textColor);
				}

				RenderText(shaderText, "Iterations: " + std::to_string(numOfIterations), 5.0f, 55.0f, 0.2f, textColor);
				RenderText(shaderText, "Particles:   " + std::to_string(NUM_PARTICLE), 5.0f, 45.0f, 0.2f, textColor);
			}		

			RenderText(shaderText, "==>", 5.0f, 35.0f, 0.2f, scaleDistanceToCenterColor);
			RenderText(shaderText, "Distance < " + std::to_string(scaleDistanceToCenter).substr(0, 4) + ": " + std::to_string(particlesWithscaleDistanceToCenter), 20.0f, 35.0f, 0.2f, textColor);
			RenderText(shaderText, "==>", 5.0f, 25.0f, 0.2f, scaleDistanceOutsideColor);
			RenderText(shaderText, "Distance > " + std::to_string(scaleDistanceOutside).substr(0, 4) + ": " + std::to_string(particlesWithscaleDistanceOutside), 20.0f, 25.0f, 0.2f, textColor);
			RenderText(shaderText, "Energy: " + std::to_string((double)energy[iterationCounter].x) + " , " + std::to_string((double)energy[iterationCounter].y) + " , " + std::to_string((double)energy[iterationCounter].z), 5.0f, 15.0f, 0.2f, textColor);
			RenderText(shaderText, "Iteration: " + std::to_string(iterationCounter), 5.0f, 5.0f, 0.2f, textColor);
		}

		if (error == 1) // If last iteration
		{
			RenderText(shaderText, "End of data reached - press R to reset", 5.0f, 600 - 20, 0.2f, glm::vec3(1.0, 0.0f, 0.0f));
		}
		if (error != 1 && slowMotion > 0) // If slowmotion on
		{
			RenderText(shaderText, "Slowmotion", 5.0f, 600 - 20, 0.2f, glm::vec3(0.9, 0.0f, 1.0f));
		}

		// Reset particle counter 
		particlesWithscaleDistanceToCenter = 0;
		particlesWithscaleDistanceOutside = 0;
			
		// Show and draw distance
		if (showDistance == 1)
		{
			drawCircle(shaderCircle, scaleDistanceToCenter, scaleDistanceToCenterColor);
			drawCircle(shaderCircle, scaleDistanceOutside, scaleDistanceOutsideColor);
		}

		// Draw instanced particles
		drawParticles(shaderInstance);
	
		// Glfw swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
		glfwSwapBuffers(window);
		glfwPollEvents();

		// If end of data is reached
		if (iterationCounter >= numOfIterations)
		{
			error = 1; // Print END OF DATA
		}
		else if(active == 1) // Else next iteration and slowmotion check
		{
			if (slowMotion == 0)
			{
				iterationCounter++;
			}
			else if (slowMotion == 1)
			{
				iterationCounter++;
				slowMotion = 2;
			}
			else if (slowMotion == 2)
			{
				slowMotion = 3;
			}
			else if (slowMotion == 3)
			{
				slowMotion = 4; // Slowmotion next iteration after 4 frames
			}
			else
			{
				slowMotion = 1;
			}
		}

		// Get new star for next iteration
		if (indexStar < NUM_PARTICLE)
		{
			indexStar++;
		}
		else
		{
			indexStar = 0;
		}

		deltaTime = glfwGetTime() - currentTime; // For fps
	}

	glfwTerminate();
	free(energy);
	energy = NULL;
	return 0;
}

void drawParticles(Shader &shader)
{
	//---READ DATA
	// Generate a new list of NUM_PARTICLE particles
	glm::vec3 *translations = (glm::vec3 *)malloc(NUM_PARTICLE * sizeof(glm::vec3));

	char file[128];

	if (iterationCounter == 1) // Initial conditions
	{
		// Read locations from initial_conditions.csv
		snprintf(file, sizeof(char) * 128, "%s\\initial_conditions.csv", dataFolder);
	}
	else // Iterations
	{
		// Read locations from iteration.csv
		snprintf(file, sizeof(char) * 128, "%s\\iteration_%i.csv", dataFolder, iterationCounter);
	}
	FILE *CSV;
	int index = 0; // Curent particle
	float x, y, z, m, vx, vy, vz;

	CSV = fopen(file, "r");
	if (CSV == NULL)
	{
		printf("Unable to open %c \n", file);
	}
	else
	{
		while ((fscanf(CSV, "%f,%f,%f,%f,%f,%f,%f\n", &x, &y, &z, &m, &vx, &vy, &vz)) > 0) // Each loop reads one row of the file
		{
			// Set position offsets for instance
			glm::vec3 translation;
			translation.x = x;
			translation.y = y;
			translation.z = z;
			translations[index] = translation;

			index++;
		}

		// Count particles outside and inside of specified scale  	
		for (int i = 0; i < NUM_PARTICLE; i++)
		{
			if (sqrt(pow(translations[i].x, 2) + pow(translations[i].y, 2) + pow(translations[i].z, 2)) < scaleDistanceToCenter)
			{
				particlesWithscaleDistanceToCenter++;
			}
			if (sqrt(pow(translations[i].x, 2) + pow(translations[i].y, 2) + pow(translations[i].z, 2)) > scaleDistanceOutside)
			{
				particlesWithscaleDistanceOutside++;
			}
		}

		if (index < NUM_PARTICLE) // If not all rows in file were counted => Error
		{
			printf("ERROR - %c: row %i - Press W to resume \n", file, index + 1); // Input-Data Error			
		}
		fclose(CSV);
	}

	// Store instance data in an array buffer
	unsigned int instanceVBO;
	glGenBuffers(1, &instanceVBO);
	glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * NUM_PARTICLE, &translations[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Store default vertice data in an array buffer
	unsigned int quadVAO, quadVBO;
	glGenVertexArrays(1, &quadVAO);
	glGenBuffers(1, &quadVBO);
	glBindVertexArray(quadVAO);
	glBindBuffer(GL_ARRAY_BUFFER, quadVBO);

	if (simpleStars > 0)
	{
		glBufferData(GL_ARRAY_BUFFER, sizeof(defaultSimpleStar), defaultSimpleStar, GL_STATIC_DRAW);
	}
	else
	{
		glBufferData(GL_ARRAY_BUFFER, sizeof(complexStar), complexStar, GL_STATIC_DRAW);
	}

	// Set default vertice 
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

	// Set instance data
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glVertexAttribDivisor(1, 1);

	free(translations); // Free allocated memory 
	//---READ DATA---END

	//---Instance
	shader.use();
	// Camera view instance
	if (mouseLook == 1)
	{
		glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp); // 1: Position of Camera, 2: Camera looks at center, 3: Up direction of camera 
		shader.setMat4("view", view);
	}
	else // No mouselook => camera looks at center
	{
		glm::mat4 view = glm::lookAt(cameraPos, glm::vec3(0.0f, 0.0f, 0.0f), cameraUp); // 1: Position of Camera, 2: Camera looks at center, 3: Up direction of camera 
		shader.setMat4("view", view);
	}

	// Projection
	glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)WIDTH / (float)HEIGHT, 0.1f, 999.0f); // Perspective projection, last parameter: far plane
	shader.setMat4("projection", projection);

	// Send indexStar to shader
	shader.setInt("indexStar", indexStar);

	// Send sparklingStar to shader
	shader.setInt("sparklingStar", sparklingStar);

	// Send sinSparkling to shader
	shader.setFloat("sinSparkling", sinSparkling);

	// Send NUM_PARTICLE to shader
	shader.setInt("NUM_PARTICLE", NUM_PARTICLE);

	// Send number of marked particles to shader
	shader.setInt("numberOfMarkedParticles", numberOfMarkedParticles);

	// Send point size to shader
	shader.setFloat("pointSize", pointSize);

	// Send color mode to shader
	shader.setInt("colorMode", colorMode);

	// Draw NUM_PARTICLE instanced quads
	glBindVertexArray(quadVAO);
	glDrawArraysInstanced(GL_POINTS, 0, numberOfPoints, NUM_PARTICLE); // NUM_PARTICLE points of numberOfPoints vertices each
	glBindVertexArray(0);
	//---Instance---END
}

void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
	// Activate corresponding render state	
	shader.use();
	glUniform3f(glGetUniformLocation(shader.ID, "textColor"), color.x, color.y, color.z);
	glActiveTexture(GL_TEXTURE0);
	glBindVertexArray(VAO);

	// Iterate through all characters
	std::string::const_iterator c;
	for (c = text.begin(); c != text.end(); c++)
	{
		Character ch = Characters[*c];

		GLfloat xpos = x + ch.Bearing.x * scale;
		GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

		GLfloat w = ch.Size.x * scale;
		GLfloat h = ch.Size.y * scale;
		// Update VBO for each character
		GLfloat vertices[6][4] = {
			{ xpos,     ypos + h,   0.0, 0.0 },
			{ xpos,     ypos,       0.0, 1.0 },
			{ xpos + w, ypos,       1.0, 1.0 },

			{ xpos,     ypos + h,   0.0, 0.0 },
			{ xpos + w, ypos,       1.0, 1.0 },
			{ xpos + w, ypos + h,   1.0, 0.0 }
		};
		// Render glyph texture over quad
		glBindTexture(GL_TEXTURE_2D, ch.TextureID);
		// Update content of VBO memory
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); 
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		// Render quad
		glDrawArrays(GL_TRIANGLES, 0, 6);
		// Advance cursors for next glyph (advance is number of 1/64 pixels)
		x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
	}
	glBindVertexArray(0);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void drawCircle(Shader &shader, float radius, glm::vec3 scaleColor)
{
	const float numberOfCircleSectors = 2 * 3.141592 / 0.05;

	shader.use();

	// View for circle
	if (mouseLook == 1)
	{
		glm::mat4 viewCircle = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp); // 1: Position of Camera, 2: Camera looks at center, 3: Up direction of camera
		shader.setMat4("view", viewCircle);
	}
	else
	{
		glm::mat4 viewCircle = glm::lookAt(cameraPos, glm::vec3(0.0f, 0.0f, 0.0f), cameraUp); // 1: Position of Camera, 2: Camera looks at center, 3: Up direction of camera
		shader.setMat4("view", viewCircle);
	}

	// Projection for circle
	glm::mat4 projectionCircle = glm::perspective(glm::radians(45.0f), (float)WIDTH / (float)HEIGHT, 0.1f, 999.0f); // perspective projection, last parameter: far plane
	shader.setMat4("projection", projectionCircle);

	// Send color to shader
	shader.setVec3("scaleColor", scaleColor);

	glm::vec3 *circle = (glm::vec3 *)malloc(numberOfCircleSectors * sizeof(glm::vec3));
	int index = 0;
	
	// Circle 		
	for (GLfloat angle = 0; angle < 2 * 3.141592; angle += 0.05f)
		{
			// Set position offsets for circle instance
			glm::vec3 translation;

			translation.x = radius*cos(angle);
			translation.y = radius*sin(angle);
			translation.z = 0;
			circle[index] = translation;

			index++;
		}

	// Store instance data in an array buffer
	unsigned int circleVBO;
	glGenBuffers(1, &circleVBO);
	glBindBuffer(GL_ARRAY_BUFFER, circleVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * numberOfCircleSectors, &circle[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// Set instance data
	unsigned int quad2VAO;
	glGenVertexArrays(1, &quad2VAO);
	glBindVertexArray(quad2VAO);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, circleVBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glVertexAttribDivisor(1, 1);
	
	// Draw
	glDrawArraysInstanced(GL_LINE_LOOP, 0, numberOfCircleSectors, numberOfCircleSectors); // numberOfCircleSectors points of numberOfCircleSectors vertices each
	glBindVertexArray(0);

	free(circle);
}

void readEnergyData(glm::vec3 *energy)
{
	//---READ ENERGY-DATA
	char file[128];

	// Read energy-data from energy_diagnostics.csv
	snprintf(file, sizeof(char) * 128, "%s\\energy_diagnostics.csv", dataFolder);

	FILE *CSV1;
	int index1 = 0;
	float x1, y1, z1;

	CSV1 = fopen(file, "r");
	if (CSV1 == NULL)
	{
		printf("Unable to open %s\\energy_diagnostics.csv \n", dataFolder);
	}
	else
	{
		while ((fscanf(CSV1, "%f,%f,%f\n", &x1, &y1, &z1)) > 0) // Each loop reads one row of the file
		{
			glm::vec3 translation;
			translation.x = x1;
			translation.y = y1;
			translation.z = z1;
			energy[index1] = translation;

			index1++;
		}

		if (index1 < numOfIterations) // If not all rows in file were counted => Error
		{
			printf("ERROR - %s\\energy_diagnostics.csv: row %i\n", dataFolder, index1 + 1); // Input-Data Error			
		}
		fclose(CSV1);
	}
	//---READ ENERGY-DATA---END
}

// Count particles 
int countParticle()
{
	char initialConditions[128];
	snprintf(initialConditions, sizeof(char) * 128, "%s\\initial_conditions.csv", dataFolder);

	int rows = 0;
	char c = 0;

	FILE *file = fopen(initialConditions, "r");

	if (file == NULL)
	{
		printf("This folder doesn't exist \n");
		printf("press enter to exit \n");
		getchar();
		exit(0);
	}
	else
	{
		while ((c = fgetc(file)) != EOF) // Check every char
		{
			if (c == '\n') // If char is \n => count as row
			{
				rows++;
			}
		}
		fclose(file);

		printf("Particles: %d\n", rows);
	}
	return rows;
}

// Count iterations 
int countIterations()
{
	int iterationCounter = 1;
	char buffer[128];

	snprintf(buffer, sizeof(char) * 128, "%s\\iteration_%i.csv", dataFolder, iterationCounter);

	while (_access(buffer, 00) != -1) // While file exists
	{
		iterationCounter++;
		snprintf(buffer, sizeof(char) * 128, "%s\\iteration_%i.csv", dataFolder, iterationCounter);
	}
	printf("Iterations: %d\n", iterationCounter - 1);
	return iterationCounter - 1;
}

// Get datafolder and check arguments
void argumentStatus(int argumentCount, char* argumentValue[])
{
	if (argumentCount > 2)
	{
		printf("Only 1 command line argument allowed");
	}
	if (argumentCount == 2)
	{
		printf("Data folder: %s \n", argumentValue[1]);
		strcpy(dataFolder, argumentValue[1]);
	}
	if (argumentCount == 1)
	{
		printf("Data folder: %s \n", dataFolder); // No arguments => default data folder "run"
	}
}

// Glfw whenever the mouse moves, this callback is called
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
	lastX = xpos;
	lastY = ypos;

	float sensitivity = 0.05f; // Mouse speed
	xoffset *= sensitivity;
	yoffset *= sensitivity;

	yaw += xoffset;
	pitch += yoffset;

	// When pitch is out of bounds, screen doesn't get flipped
	if (pitch > 89.0f)
		pitch = 89.0f;
	if (pitch < -89.0f)
		pitch = -89.0f;

	glm::vec3 front;
	front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
	front.y = sin(glm::radians(pitch));
	front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
	cameraFront = glm::normalize(front);
}

// Get hotkey input
void processInput(GLFWwindow *window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) // ESC pressed => close window
	{
		glfwSetWindowShouldClose(window, true);
	}

	// Fullscreen/window
	if (glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS && glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS) // Alt + enter for fullscreen
	{
		if (fullScreen == 0)
		{
			windowedWIDTH = WIDTH; // Save old window size for change from fullscreen to windowed
			windowedHEIGHT = HEIGHT;
			fullScreen = 1;
			const GLFWvidmode * mode = glfwGetVideoMode(glfwGetPrimaryMonitor()); // Get current desktop resolution
			glfwSetWindowMonitor(window, glfwGetPrimaryMonitor(), 0, 0, mode->width, mode->height, GLFW_DONT_CARE); // Fullscreen
		}
		else
		{
			fullScreen = 0;
			glfwSetWindowMonitor(window, NULL, 10, 50, windowedWIDTH, windowedHEIGHT, GLFW_DONT_CARE); // Windowed mode
		}
	}

	// Zoom
	float distanceToCenter = sqrt(pow(cameraPos.x, 2) + pow(cameraPos.y, 2) + pow(cameraPos.z, 2));

	if (mouseLook == 0)
	{
		if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
		{
			// Zoom speed depends on distance
			if (zoom > 0.001f)
			{
				zoom += 0.001f - distanceToCenter * 0.01;
			}
			cameraPos = glm::vec3(sin(alpha)*zoom, 0.0f, cos(-alpha)*zoom); // x , y , z, trigonometry
		}

		if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
		{
			if (zoom < 950.00f)
			{
				zoom -= 0.001f - distanceToCenter * 0.01;
			}
			cameraPos = glm::vec3(sin(alpha)*zoom, 0.0f, cos(-alpha)*zoom);
		}

		// Rotate
		if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
		{
			alpha -= 0.01f; // Rotate speed

			float camX = sin(alpha);
			float camZ = cos(alpha);

			cameraPos = glm::vec3(camX*zoom, 0.0f, camZ*zoom); // x , y , z, trigonometry
		}

		if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
		{
			alpha += 0.01f;

			float camX = sin(alpha);
			float camZ = cos(alpha);

			cameraPos = glm::vec3(camX*zoom, 0.0f, camZ*zoom);
		}
	}
	else // Mouselook on speed depends on distance
	{
		float cameraSpeed = 0.1 * deltaTime;

		if (slowMotion == 0)
		{
			cameraSpeed = 0.1 * deltaTime + distanceToCenter * 0.01;
		}
		else // Slower movement when slowmotion speed depends on distance
		{
			cameraSpeed = 0.01 * deltaTime + distanceToCenter * 0.0025;
		}
		
		if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
		{
			cameraPos += cameraSpeed * cameraFront;
		}
		if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
		{
			cameraPos -= cameraSpeed * cameraFront;
		}
		if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
		{
			cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
		}
		if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
		{
			cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) // Reset
	{
		iterationCounter = 1;
		if (error = 1)
		{
			error = 0;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) // Point size and simple star shape
	{
		pointSize = 0.3;
		simpleStars = 1;
	}

	if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) // Point size and complex star shape
	{
		pointSize = 0;
		simpleStars = 0;
	}

	if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) // Point size and complex star shape
	{
		pointSize = 1;
		simpleStars = 0;		
	}
	
	if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS) // Point size and complex star shape
	{
		pointSize = 2;
		simpleStars = 0;
	}

	int newStateSpace = glfwGetKey(window, GLFW_KEY_SPACE); 
	if (newStateSpace == GLFW_RELEASE && oldStateSpace == GLFW_PRESS) // Start/stop
	{
		if (active == 0)
		{
			active = 1;
		}
		else
		{
			active = 0;
		}			
	}
	oldStateSpace = newStateSpace;

	int newStateS = glfwGetKey(window, GLFW_KEY_S);
	if (newStateS == GLFW_RELEASE && oldStateS == GLFW_PRESS) // Slowmotion 25%
	{
		if (slowMotion == 0)
		{
			slowMotion = 1;
		}
		else
		{
			slowMotion = 0;
		}
	}
	oldStateS = newStateS;

	int newStateM = glfwGetKey(window, GLFW_KEY_M);
	if (newStateM == GLFW_RELEASE && oldStateM == GLFW_PRESS) // Mouselook
	{
		if (mouseLook == 0)
		{
			mouseLook = 1;
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED); // Hide cursor
		}
		else
		{
			mouseLook = 0;
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL); // Show cursor
		}
	}
	oldStateM = newStateM;

	int newStatePAGE_UP = glfwGetKey(window, GLFW_KEY_PAGE_UP);
	if (newStatePAGE_UP == GLFW_RELEASE && oldStatePAGE_UP == GLFW_PRESS) // Mark particle
	{
		if (numberOfMarkedParticles < NUM_PARTICLE)
		{
			numberOfMarkedParticles++;
		}		
	}
	oldStatePAGE_UP = newStatePAGE_UP;

	int newStatePAGE_DOWN = glfwGetKey(window, GLFW_KEY_PAGE_DOWN);
	if (newStatePAGE_DOWN == GLFW_RELEASE && oldStatePAGE_DOWN == GLFW_PRESS) // Hide particle
	{
		if (numberOfMarkedParticles > 0)
		{
			numberOfMarkedParticles--;
		}
	}
	oldStatePAGE_DOWN = newStatePAGE_DOWN;

	int newStateD = glfwGetKey(window, GLFW_KEY_D);
	if (newStateD == GLFW_RELEASE && oldStateD == GLFW_PRESS) // Show / hide data
	{
		if (showData == 0)
		{
			showData = 1;
		}
		else
		{
			showData = 0;
		}
	}
	oldStateD = newStateD;

	int newStateF = glfwGetKey(window, GLFW_KEY_F);
	if (newStateF == GLFW_RELEASE && oldStateF == GLFW_PRESS) // Show / hide scale
	{
		if (showDistance == 0)
		{
			showDistance = 1;
		}
		else
		{
			showDistance = 0;
		}
	}
	oldStateF = newStateF;

	int newStateKP_1 = glfwGetKey(window, GLFW_KEY_KP_1);
	if (newStateKP_1 == GLFW_RELEASE && oldStateKP_1 == GLFW_PRESS) // Decrease inner scale
	{
		if (scaleDistanceToCenter >= 5.00f)
		{
			scaleDistanceToCenter -= 1.00f;
		}
		else if (scaleDistanceToCenter > 1.00f)
		{
			scaleDistanceToCenter -= 0.50f;
		}
		else if (scaleDistanceToCenter >= 0.10f)
		{
			scaleDistanceToCenter -= 0.05f;
		}
	}
	oldStateKP_1 = newStateKP_1;

	int newStateKP_4 = glfwGetKey(window, GLFW_KEY_KP_4);
	if (newStateKP_4 == GLFW_RELEASE && oldStateKP_4 == GLFW_PRESS) // Enlarge inner scale
	{

		if (scaleDistanceToCenter >= 5.00f)
		{
			scaleDistanceToCenter += 1.00f;
		}
		else if (scaleDistanceToCenter >= 1.00f)
		{
			scaleDistanceToCenter += 0.50f;
		}
		else if (scaleDistanceToCenter >= 0.00f)
		{
			scaleDistanceToCenter += 0.05f;
		}
	}
	oldStateKP_4 = newStateKP_4;

	int newStateKP_2 = glfwGetKey(window, GLFW_KEY_KP_2);
	if (newStateKP_2 == GLFW_RELEASE && oldStateKP_2 == GLFW_PRESS) // Decrease outer scale
	{
		if (scaleDistanceOutside > 5.00f)
		{
			scaleDistanceOutside -= 1.00f;
		}
		else if (scaleDistanceOutside > 1.00f)
		{
			scaleDistanceOutside -= 0.50f;
		}
		else if (scaleDistanceOutside > 0.10)
		{
			scaleDistanceOutside -= 0.05f;
		}
	}
	oldStateKP_2 = newStateKP_2;

	int newStateKP_5 = glfwGetKey(window, GLFW_KEY_KP_5);
	if (newStateKP_5 == GLFW_RELEASE && oldStateKP_5 == GLFW_PRESS) // Enlarge outer scale
	{
		if (scaleDistanceOutside >= 5.00f)
		{
			scaleDistanceOutside += 1.00f;
		}
		else if (scaleDistanceOutside >= 1.00f)
		{
			scaleDistanceOutside += 0.50f;
		}
		else if (scaleDistanceOutside >= 0.00f)
		{
			scaleDistanceOutside += 0.05f;
		}
	}

	oldStateKP_5 = newStateKP_5;

	int newStateC = glfwGetKey(window, GLFW_KEY_C);
	if (newStateC == GLFW_RELEASE && oldStateC == GLFW_PRESS) // Color mode
	{
		if (colorMode == 0)
		{
			colorMode = 1;
		}
		else
		{
			colorMode = 0;
		}
	}
	oldStateC = newStateC;

	int newStateP = glfwGetKey(window, GLFW_KEY_P);
	if (newStateP == GLFW_RELEASE && oldStateP == GLFW_PRESS) // Experimental painting mode
	{
		if (funPaintMode == 0)
		{
			funPaintMode = 1;
			colorMode = 1;
		}
		else
		{
			funPaintMode = 0;
			colorMode = 0;
			waitOneFrame = 0;
		}
	}
	oldStateP = newStateP;
}

// Glfw whenever the window size changed this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
	if (width != 0 && height != 0)
	{
		WIDTH = width;
		HEIGHT = height;
	}
}
