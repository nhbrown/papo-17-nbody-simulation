#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aOffset;

out vec3 fColor;
flat out int InstanceID;
flat out int indexStarFragment;
flat out int colorModeFragment;
flat out int sparklingStarFragment;
flat out float sinSparklingFragment;

uniform mat4 view;
uniform mat4 projection;
uniform int indexStar;
uniform int sparklingStar;
uniform float sinSparkling;
uniform int NUM_PARTICLE;
uniform int numberOfMarkedParticles;
uniform float pointSize;
uniform int colorMode;

float sparkling = 0;
vec3 blueLightblue = vec3(0.04,0.65+sin(gl_InstanceID)*0.35,1.0);
vec3 lightblueWhite = vec3(0.6+sin(gl_InstanceID)*0.4,1.0,1.0);
vec3 bluePurple = vec3(0.15+sin(gl_InstanceID)*0.4,0.0,1);

void main()
{
	InstanceID = gl_InstanceID;
	indexStarFragment = indexStar;
	colorModeFragment = colorMode;
	sparklingStarFragment = sparklingStar;
	sinSparklingFragment = sinSparkling;

	if (sparklingStar == InstanceID)
	{
		sparkling = abs(sin(sinSparkling)) * 1;
	}
	else
	{
		sparkling = 0;
	}

	if(colorMode < 1)
	{
		if(gl_InstanceID % 2 == 0 || gl_InstanceID % 3 == 0)
		{			
			fColor = blueLightblue;
		}
		else if(gl_InstanceID % 5 == 0 || gl_InstanceID % 7 == 0 || gl_InstanceID % 11 == 0 || gl_InstanceID % 13 == 0)
		{
			fColor = lightblueWhite;
		}
		else
		{
			fColor = bluePurple;
		}
	}
	else
	{
		sparkling = 0;

		if(gl_InstanceID % 2 == 0)
		{
			fColor = vec3(sin(gl_InstanceID), cos(gl_InstanceID), tan(gl_InstanceID));
		}
		else if(gl_InstanceID % 3 == 0 || gl_InstanceID % 5 == 0)
		{
			fColor = vec3(tan(gl_InstanceID), cos(gl_InstanceID), sin(gl_InstanceID));
		}
		else
		{
			fColor = vec3(cos(gl_InstanceID)*2, sin(gl_InstanceID)*2, 0);
		}
	}

	if(gl_InstanceID < NUM_PARTICLE * 0.7)
	{
		gl_PointSize = 0.7 + sparkling + pointSize;
	}
	else if(gl_InstanceID > NUM_PARTICLE * 0.7 && gl_InstanceID < NUM_PARTICLE * 0.9)
	{
		gl_PointSize = 1.5 + sparkling + pointSize;
	}
	else
	{
		gl_PointSize = 2.5 + sparkling + pointSize;
	}

	if(gl_InstanceID < numberOfMarkedParticles)
	{
		fColor = vec3(1,0,0);
		gl_PointSize = 3 + pointSize;
	}
	
    gl_Position = projection * view * vec4(aPos + aOffset, 1.0);
}