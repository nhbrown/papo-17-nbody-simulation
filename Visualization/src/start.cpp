#define _CRT_SECURE_NO_WARNINGS //Disable some warnings
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <opengl/shader.h>

#include <io.h>      /* Count files(iterations) MS-DOS , #include <unistd.h> for UNIX/LINUX */

//OpenGL Mathematics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

void framebuffer_size_callback(GLFWwindow* window, int width, int height); //Callback when window is resized
void processInput(GLFWwindow *window); //Callback when keys are pressed

char dataFolder[128] = "run"; //Default folder for data

//counting of particles 
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
		while ((c = fgetc(file)) != EOF) //Check every char
		{
			if (c == '\n') //if char is \n => count as row
			{
				rows++;
			}
		}
		fclose(file);

		printf("Particles: %d\n", rows);
	}
	return rows;
}

//counting of iterations 
int countIterations()
{
	int iterationCounter = 1;
	char buffer[128];

	snprintf(buffer, sizeof(char) * 128, "%s\\iteration_%i.csv", dataFolder, iterationCounter);

	while (_access(buffer, 00) != -1) //While file exists
	{
		iterationCounter++;
		snprintf(buffer, sizeof(char) * 128, "%s\\iteration_%i.csv", dataFolder, iterationCounter);
	}
	printf("Iterations: %d\n", iterationCounter - 1);
	return iterationCounter - 1;
}

//---Settings
int pointSize = 1; //Size of particles
int antiAliasing = 1;
int effectStatus = 1;

//Default color
float red = 1.0f; //1.0f corresponds 255 in RGB
float green = 0.0f;
float blue = 0.0f;

//camera settings
float alpha = 0.0f; //Default alpha angle
float zoom = 3.0f; //Default zoom

//Default camera position
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, zoom);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

//Various settings
int active = 0; //Reading of iteration data on/off
unsigned int iterationCounter = 1; //Current iteration
int testCounter = 0; //Counting of Stars for Console-Status

int main(int argc, char* argv[])
{
	if (argc > 2)
	{
		printf("Only 1 command line argument allowed");
	}
	if (argc == 2)
	{
		printf("Data folder: %s \n", argv[1]);
		strcpy(dataFolder, argv[1]);
	}
	if (argc == 1)
	{
		printf("Data folder: %s \n", dataFolder); //No arguments => default data folder "run"
	}

	//Settings
	unsigned int NUM_PARTICLE = countParticle(); //set number of particles 
	unsigned const int numOfIterations = countIterations(); //set number of Iterations
	const unsigned int SCR_WIDTH = 1920; //Start screen resolution
	const unsigned int SCR_HEIGHT = 1080; //Start screen resolution

	//Default vertex
	float vertices[] =
	{
		//center
		0.000000f,  0.000000f, 0.000000f,

		0.0000f,  0.0001f, 0.0000f,
		0.0000f,  -0.0001f, 0.0000f,
		0.0001f,  0.0000f, 0.0000f,
		-0.0001f,  0.0000f, 0.0000f,
		0.0000f,  0.0000f, 0.0001f,
		0.0000f,  0.0000f, -0.0001f,

		0.0000f,  0.0002f, 0.0000f,
		0.0000f,  -0.0002f, 0.0000f,
		0.0002f,  0.0000f, 0.0000f,
		-0.0002f,  0.0000f, 0.0000f,
		0.0000f,  0.0000f, 0.0002f,
		0.0000f,  0.0000f, -0.0002f,

		0.0000f,  0.0003f, 0.0000f,
		0.0000f,  -0.0003f, 0.0000f,
		0.0003f,  0.0000f, 0.0000f,
		-0.0003f,  0.0000f, 0.0000f,
		0.0000f,  0.0000f, 0.0003f,
		0.0000f,  0.0000f, -0.0003f,

		//cross
		0.000000f,  -0.001f, 0.000000f,
		0.000000f,  0.001f, 0.000000f,
		0.001f,  0.000000f, 0.000000f,
		-0.001f,  0.000000f, 0.000000f,
		0.000000f,  0.000000f, 0.001f,
		0.000000f,  0.000000f, -0.001f,

		//inner cross
		0.000000f,  -0.0005f, 0.000000f,
		0.000000f,  0.0005f, 0.000000f,
		0.0005f,  0.000000f, 0.000000f,
		-0.0005f,  0.000000f, 0.000000f,
		0.000000f,  0.000000f, 0.0005f,
		0.000000f,  0.000000f, -0.0005f,

		//fillment center
		0.000000f,  -0.0001f, -0.0001f,
		0.000000f,  0.0001f, -0.0001f,
		0.000000f,  -0.0001f, 0.0001f,
		0.000000f,  0.0001f, 0.0001f,

		0.0001f,  0.000000f, -0.0001f,
		-0.0001f,  0.000000f, -0.0001f,
		0.0001f,  0.000000f, 0.0001f,
		-0.0001f,  0.000000f, 0.0001f,

		0.0001f,  -0.0001f, 0.000000f,
		-0.0001f,  -0.0001f, 0.000000f,
		0.0001f,  0.0001f, 0.000000f,
		-0.0001f,  0.0001f, 0.000000f,

		-0.0001f,  -0.0001f, -0.0001f,
		-0.0001f,  -0.0001f, +0.0001f,
		-0.0001f,  +0.0001f, +0.0001f,
		+0.0001f,  +0.0001f, +0.0001f,
		+0.0001f,  +0.0001f, -0.0001f,
		+0.0001f,  -0.0001f, -0.0001f,
		-0.0001f,  +0.0001f, -0.0001f,
		+0.0001f,  -0.0001f, +0.0001f,

		//fillment center2
		0.000000f,  -0.0002f, -0.0002f,
		0.000000f,  0.0002f, -0.0002f,
		0.000000f,  -0.0002f, 0.0002f,
		0.000000f,  0.0002f, 0.0002f,

		0.0002f,  0.000000f, -0.0002f,
		-0.0002f,  0.000000f, -0.0002f,
		0.0002f,  0.000000f, 0.0002f,
		-0.0002f,  0.000000f, 0.0002f,

		0.0002f,  -0.0002f, 0.000000f,
		-0.0002f,  -0.0002f, 0.000000f,
		0.0002f,  0.0002f, 0.000000f,
		-0.0002f,  0.0002f, 0.000000f,

		-0.0002f,  -0.0002f, -0.0002f,
		-0.0002f,  -0.0002f, +0.0002f,
		-0.0002f,  +0.0002f, +0.0002f,
		+0.0002f,  +0.0002f, +0.0002f,
		+0.0002f,  +0.0002f, -0.0002f,
		+0.0002f,  -0.0002f, -0.0002f,
		-0.0002f,  +0.0002f, -0.0002f,
		+0.0002f,  -0.0002f, +0.0002f,

		//fillment
		0.000000f,  -0.0005f, -0.0005f,
		0.000000f,  0.0005f, -0.0005f,
		0.000000f,  -0.0005f, 0.0005f,
		0.000000f,  0.0005f, 0.0005f,

		0.0005f,  0.000000f, -0.0005f,
		-0.0005f,  0.000000f, -0.0005f,
		0.0005f,  0.000000f, 0.0005f,
		-0.0005f,  0.000000f, 0.0005f,

		0.0005f,  -0.0005f, 0.000000f,
		-0.0005f,  -0.0005f, 0.000000f,
		0.0005f,  0.0005f, 0.000000f,
		-0.0005f,  0.0005f, 0.000000f,

		-0.0005f,  -0.0005f, -0.0005f,
		-0.0005f,  -0.0005f, +0.0005f,
		-0.0005f,  +0.0005f, +0.0005f,
		+0.0005f,  +0.0005f, +0.0005f,
		+0.0005f,  +0.0005f, -0.0005f,
		+0.0005f,  -0.0005f, -0.0005f,
		-0.0005f,  +0.0005f, -0.0005f,
		+0.0005f,  -0.0005f, +0.0005f,

		//cross outside
		0.000f, 0.001f, 0.000f,
		0.000f, -0.001f, 0.000f,
		0.001f, 0.000f, 0.000f,
		-0.001f, 0.000f, 0.000f,
		0.000f, 0.000f, 0.001f,
		0.000f, 0.000f, -0.001f,

			0.000f, 0.002f, 0.000f,
			0.000f, -0.002f, 0.000f,
			0.002f, 0.000f, 0.000f,
			-0.002f, 0.000f, 0.000f,
			0.000f, 0.000f, 0.002f,
			0.000f, 0.000f, -0.002f,



	};
	int numberOfPoints = 105; //Number of points for particle shape (number of vertices)

	//glfw initialize and configure
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, antiAliasing); // Antialiasing
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X

	//glfw window creation
	GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "NBodyVisualization", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback); //Window resize

	//glad load all OpenGL function pointers
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	//configure global opengl state
	glEnable(GL_MULTISAMPLE); //Enable antialiasing 
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_PROGRAM_POINT_SIZE);

	//build and compile shaders
	Shader shader("shader/vertex_shader.vs", "shader/fragment_shader.fs");
	
	// generate a new list of NUM_PARTICLE particles
	glm::vec3 *translations = (glm::vec3 *)malloc(NUM_PARTICLE * sizeof(glm::vec3));

	//generate a list of NUM_PARTICLE particles for colors
	glm::vec3 *colorTranslations = (glm::vec3 *)malloc(NUM_PARTICLE * sizeof(glm::vec3));

	//Read locations from initial_conditions.csv
	char initialConditions[128];
	snprintf(initialConditions, sizeof(char) * 128, "%s\\initial_conditions.csv", dataFolder);

	FILE *CSV;
	unsigned int index = 0; //Curent particle
	float x, y, z, m, vx, vy, vz;

	CSV = fopen(initialConditions, "r");
	if (CSV == NULL)
	{
		printf("Unable to open initial_conditions.csv \n");
	}
	else
	{
		while ((fscanf(CSV, "%f,%f,%f,%f,%f,%f,%f\n", &x, &y, &z, &m, &vx, &vy, &vz)) > 0) //Each loop reads one row of the initial_conditions.csv
		{
			//Set position offsets for instance
			glm::vec3 translation;
			translation.x = x;
			translation.y = y;
			translation.z = z;
			translations[index] = translation;

			//Set color offsets 
			glm::vec3 colorTranslation;
			colorTranslation.x = x;
			colorTranslation.y = y;
			colorTranslation.z = z;
			colorTranslations[index] = colorTranslation;

			index++;
		}

		if (index < NUM_PARTICLE) //If not all rows in iteration.cvs were counted => Error
		{
			printf("ERROR - initial_conditions.csv: row %i - Press W to resume \n", index+1); //Input-Data Error
			active = 0; //Stop
		}
		fclose(CSV);
	}

	//store instance data in an array buffer
	unsigned int instanceVBO;
	glGenBuffers(1, &instanceVBO);
	glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * NUM_PARTICLE, &translations[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//store instance color data in an array buffer
	unsigned int instanceVBOColor;
	glGenBuffers(1, &instanceVBOColor);
	glBindBuffer(GL_ARRAY_BUFFER, instanceVBOColor);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * NUM_PARTICLE, &colorTranslations[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	free(translations); //Free allocated memory
	free(colorTranslations);

	

	//Configure vertex attributes
	unsigned int VAO, VBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(3 * sizeof(float)));
	
	//set instance data
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, instanceVBOColor); // this attribute comes from a different color vertex buffer
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(3 * sizeof(float)));
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glVertexAttribDivisor(1, 1); // tell OpenGL this is an instanced vertex attribute.

	glEnableVertexAttribArray(2);
	glBindBuffer(GL_ARRAY_BUFFER, instanceVBO); // this attribute comes from a different vertex buffer
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(6 * sizeof(float)));
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glVertexAttribDivisor(2, 1); // tell OpenGL this is an instanced vertex attribute.

	//Print operation manual
	printf("Use W to start, S to stop, X to reset, UP/DOWN/LEFT/RIGHT to zoom and rotate, 1/2/3 to change point size and E/R/F for effects\n");

	//render loop
	while (!glfwWindowShouldClose(window))
	{
		//Key input
		processInput(window);

		//Reading of iteration data
		if (active == 1)
		{
			if (iterationCounter < numOfIterations)
			{
				printf("Iteration: %d - Particles with distance to center > 1: %d\n", iterationCounter, testCounter); // Print current iteration if not last iteration
			}
			testCounter = 0; //Reset counter 

			//Delete old vertex and buffers
			glDeleteVertexArrays(1, &VAO);
			glDeleteBuffers(1, &VBO);

			// generate a new list of NUM_PARTICLE particles
			glm::vec3 *translations = (glm::vec3 *)malloc(NUM_PARTICLE * sizeof(glm::vec3));

			//generate a list of NUM_PARTICLE particles for colors
			glm::vec3 *colorTranslations = (glm::vec3 *)malloc(NUM_PARTICLE * sizeof(glm::vec3));

			//Read locations from iterations
			FILE *CSV;
			unsigned int index = 0; //Current particle
			float x, y, z, m, vx, vy, vz;

			//Set folder
			char buffer[128];
			snprintf(buffer, sizeof(char) * 128, "%s\\iteration_%i.csv", dataFolder, iterationCounter);

			CSV = fopen(buffer, "r");
			if (CSV == NULL)
			{
				printf("Unable to open %s \n", buffer);
			}
			else
			{
					while ((fscanf(CSV, "%f,%f,%f,%f,%f,%f,%f\n", &x, &y, &z, &m, &vx, &vy, &vz)) > 0 ) //Each loop reads one row of the iteration.csv (-1 : EOF , 0 : not floats)
					{						
							//Set position offsets for instance
							glm::vec3 translation;
							translation.x = x;
							translation.y = y;
							translation.z = z;
							translations[index] = translation;					

							int randomForSparkling = sin(glfwGetTime());

							
							//Set color offsets 
							if (index % 2 == 0) //If  then sparkle => not every star will sparkle
							{
								if (randomForSparkling < 0.4)
								{
									randomForSparkling = randomForSparkling + 0.4;
								}

								glm::vec3 colorTranslation;
								colorTranslation.x = x; //sparkle Stars
								colorTranslation.y = y;
								colorTranslation.z = randomForSparkling*z;
								colorTranslations[index] = colorTranslation;
							}
							else
							{
								glm::vec3 colorTranslation;
								colorTranslation.x = x; 
								colorTranslation.y = y;
								colorTranslation.z = z;
								colorTranslations[index] = colorTranslation;
							}
							
							index++;
					}

					 //Count stars outside of inner circle with range 1.0	
					for (int i = 0; i < NUM_PARTICLE;i++)
					{
						if (abs(translations[i].x) > 1 || abs(translations[i].y) > 1 || abs(translations[i].z) > 1)
						{
							testCounter++;
						}
					}
				

					if (index < NUM_PARTICLE) //If not all rows in iteration.cvs were counted => Error
					{
						printf("ERROR - iteration_%i.csv: row %i - Press W to resume \n", iterationCounter, index); //Input-Data Error
						active = 0; //Stop
					}
					fclose(CSV);								

				//If end of data is reached
				if (iterationCounter >= numOfIterations)
				{
					printf("Iteration %d - end of data reached - press X to reset\n", iterationCounter);
				}
				else //Else next iteration
				{
					iterationCounter++;
				}
			}		

			//store instance position data in an array buffer
			unsigned int instanceVBO;
			glGenBuffers(1, &instanceVBO);
			glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
			glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * NUM_PARTICLE, &translations[0], GL_STATIC_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, 0);

			//store instance color data in an array buffer
			unsigned int instanceVBOColor;
			glGenBuffers(1, &instanceVBOColor);
			glBindBuffer(GL_ARRAY_BUFFER, instanceVBOColor);
			glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * NUM_PARTICLE, &colorTranslations[0], GL_STATIC_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, 0);

			free(translations); //Free allocated memory
			free(colorTranslations);

			
			//Configure vertex attributes
			unsigned int VAO, VBO;
			glGenVertexArrays(1, &VAO);
			glGenBuffers(1, &VBO);
			glBindVertexArray(VAO);
			glBindBuffer(GL_ARRAY_BUFFER, VBO);
			glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0); //Vertex
			// also set instance data
			glEnableVertexAttribArray(1);
			glBindBuffer(GL_ARRAY_BUFFER, instanceVBOColor); // this attribute comes from a different color vertex buffer
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(3 * sizeof(float))); //Color
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glVertexAttribDivisor(1, 1); // tell OpenGL this is an instanced vertex attribute.
	
			glEnableVertexAttribArray(2);
			glBindBuffer(GL_ARRAY_BUFFER, instanceVBO); // this attribute comes from a different vertex buffer
			glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(6 * sizeof(float))); //Position offset
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glVertexAttribDivisor(2, 1); // tell OpenGL this is an instanced vertex attribute.
		}



		//render
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f); //Clear and set color of screen every frame
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//Camera view 
		glm::mat4 view = glm::lookAt(cameraPos, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)); //1: Position of Camera, 2: Camera looks at center, 3: Up direction of camera 
		shader.setMat4("view", view);

		//Projection
		glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 999.0f); //perspective projection, last parameter: far plane
		shader.setMat4("projection", projection);

		//Color
		int vertexColorLocation = glGetUniformLocation(shader.ID, "ourColor"); //get adress of color uniform variable in vertex_shader
		glUniform3f(vertexColorLocation, red, green, blue); //set color  

		//Point size
		int pointSizeLocation = glGetUniformLocation(shader.ID, "pointSize"); //get adress of pointSize uniform variable in vertex_shader
		glUniform1i(pointSizeLocation, pointSize); //set pointSize 

		
		//Send random 0-1 to Shader for effects
		float randSinus = abs(sin(10*glfwGetTime()));
		if (randSinus > 0.5)
		{
			randSinus = randSinus - 0.4;
		}
		float timeLocation = glGetUniformLocation(shader.ID, "random");
		glUniform1f(timeLocation, randSinus);
		

		//Send effectStatus
		int effectLocation = glGetUniformLocation(shader.ID, "effects");
		glUniform1i(effectLocation, effectStatus);
		
		//draw NUM_PARTICLE particles instanced 
		shader.use();
		glBindVertexArray(VAO);
		glDrawArraysInstanced(GL_POINTS, 0, numberOfPoints, NUM_PARTICLE); // draw NUM_PARTICLE instanced points 
		glBindVertexArray(0);

		//glfw swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	//deallocate all resources
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);

	//End
	glfwTerminate();
	return 0;
}

//get window input
void processInput(GLFWwindow *window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) //ESC pressed => close window
		glfwSetWindowShouldClose(window, true);

	//Zoom
	if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
	{
		//Zoom speed depends on distance
		if (zoom > 5.0f)
		{
			zoom -= 0.1f;
		}
		else if (zoom > 0.1f)
		{
			zoom -= 0.01f;
		}
		else if (zoom >= 0.01f)
		{
			zoom -= 0.001f;
		}
		cameraPos = glm::vec3(sin(alpha)*zoom, 0.0f, cos(-alpha)*zoom); //x , y , z, trigonometry
	}

	if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
	{
		if (zoom > 5.0f)
		{
			zoom += 0.1f;
		}
		else if (zoom > 0.1f)
		{
			zoom += 0.01f;
		}
		else if (zoom > 0.0f)
		{
			zoom += 0.001f;
		}
		cameraPos = glm::vec3(sin(alpha)*zoom, 0.0f, cos(-alpha)*zoom);
	}

	//Rotate
	if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
	{
		alpha -= 0.01f; //Rotate speed

		float camX = sin(alpha);
		float camZ = cos(alpha);

		cameraPos = glm::vec3(camX*zoom, 0.0f, camZ*zoom); //x , y , z, trigonometry
	}

	if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
	{
		alpha += 0.01f;

		float camX = sin(alpha);
		float camZ = cos(alpha);

		cameraPos = glm::vec3(camX*zoom, 0.0f, camZ*zoom);
	}

	//Colors
	if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
	{
		if (red >= 0.99f)
		{
			red = 0.0f;
		}
		else
		{
			red += 0.002f;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
	{
		if (green >= 0.99f)
		{
			green = 0.0f;
		}
		else
		{
			green += 0.002f;
		}
	}

	if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS)
	{
		if (blue >= 0.99f)
		{
			blue = 0.0f;
		}
		else
		{
			blue += 0.002f;
		}
	}

	//Start reading of iteration data 
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
	{
		active = 1;
	}

	//Stop reading of iteration data 
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
	{
		active = 0;
	}

	//Reset to 1 iteration
	if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
	{
		iterationCounter = 1;

		if (active == 0)
		{
			printf("Press W to continue with 1 iteration \n");
		}
	}
	//pointSize  1 
	if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
	{
		pointSize = 1;
	}
	//pointSize  2 
	if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
	{
		pointSize = 2;
	}
	//pointSize  3 
	if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
	{
		pointSize = 3;
	}

	//Effects on 
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
	{
		effectStatus = 1;
	}

	//Effects off 
	if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
	{
		effectStatus = 0;
	}

	//Effects color 
	if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
	{
		effectStatus = 2;
	}


}

//glfw whenever the window size changed this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}
