#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 offset;
layout (location = 2) in float size;
layout (location = 3) in vec4 color;

out vec4 particleColor;

uniform mat4 projection;
uniform mat4 view;

void main()
{
    vec4 worldPosition = vec4(aPos * size + offset, 1.0);
    gl_Position = projection * view * worldPosition;
    particleColor = color;
}
