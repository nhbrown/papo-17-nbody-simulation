#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec3 aOffset;

out vec3 fColor;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    fColor = aColor;
    gl_Position = projection * view * vec4(aPos + aOffset, 1.0);
}