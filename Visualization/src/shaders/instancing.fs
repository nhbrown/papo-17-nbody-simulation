#version 330 core
out vec4 FragColor;

in vec3 fColor;
flat in int InstanceID;
flat in int indexStarFragment;
flat in int sparklingStarFragment;
flat in float sinSparklingFragment;
flat in int colorModeFragment;

void main()
{
	if(colorModeFragment < 1)
	{
		if (sparklingStarFragment == InstanceID+1 || sparklingStarFragment == InstanceID+3 || sparklingStarFragment == InstanceID-15)
		{
			FragColor = vec4(fColor, abs(sin(sinSparklingFragment)) * 2); // Sparkling
		}
		else
		{
			FragColor = vec4(fColor, sin(InstanceID)*0.4+0.6);
		}
	}
	else
	{
		FragColor = vec4(fColor, 1);
	}
}