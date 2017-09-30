#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec3 aOffset;

out vec3 fColor;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 ourColor;
uniform int pointSize;
uniform float random;
uniform int effects;

void main()
{
	if (effects == 1)
	{
		fColor = vec3(0.04*abs(1/aColor.x), 0.4*abs(1/aColor.x), abs(1/aColor.x)); 
    }
	else if (effects == 2)
	{
		fColor = vec3(aColor.x*0.8, aColor.y*0.8, aColor.z*0.8); 
	}
	else
	{
		fColor = vec3(0.04, 0.5, 1); 
	}
	
	gl_Position = projection * view * vec4(aPos + aOffset, 1.0);

	if((aColor.x + aColor.y) < 1)
	{
	gl_PointSize = pointSize;
	}
	else if((aColor.x + aColor.z) < 1)
	{
	gl_PointSize = pointSize+0.5;
	}
	else
	{
	gl_PointSize = pointSize+0.8;
	}	
}