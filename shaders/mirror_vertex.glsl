#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aOffset;
layout (location = 2) in vec4 aColor;
layout (location = 3) in float aSize;

out vec4 particleColor;
out vec3 reflectedPos;

uniform mat4 projection;
uniform mat4 view;
uniform vec3 mirrorNormal;
uniform vec3 mirrorPoint;

void main()
{
    particleColor = vec4(1.0, 1.0, 1.0, 1.0);

    // Compute the reflection of the particle position
    vec3 toParticle = aOffset - mirrorPoint;
    reflectedPos = aOffset - 2.0 * dot(toParticle, mirrorNormal) * mirrorNormal;

    // Transform the original particle position
    vec4 worldPos = vec4(aOffset + aPos * aSize, 1.0);
    gl_Position = projection * view * worldPos;

    // Transform the reflected particle position
    vec4 reflectedWorldPos = vec4(reflectedPos + aPos * aSize, 1.0);
    gl_Position = projection * view * reflectedWorldPos;
}