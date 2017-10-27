#ifndef HEADER_H_
#define HEADER_H_

void framebuffer_size_callback(GLFWwindow* window, int width, int height); // Callback when window is resized
void argumentStatus(int argumentCount, char* argumentValue[]); // Argument check
void processInput(GLFWwindow *window); // Callback when keys are pressed
void mouse_callback(GLFWwindow* window, double xpos, double ypos); // Callback when mouse is moved
void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color); // Render text
void drawParticles(Shader &shader); // Draw instanced particles
void drawCircle(Shader &shader, float radius, glm::vec3 scaleColor);
void readEnergyData(glm::vec3 *energy);
int countParticle();
int countIterations();

#endif // HEADER_H_