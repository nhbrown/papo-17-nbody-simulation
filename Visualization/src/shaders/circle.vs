#version 330 core
layout (location = 0) in vec3 aPos;

out vec3 fColor;

uniform mat4 view;
uniform mat4 projection;
uniform vec3 scaleColor;

void main()
{
	fColor = vec3(scaleColor);
    gl_Position = projection * view * vec4(aPos, 1.0);		
}