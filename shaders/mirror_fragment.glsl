#version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D screenTexture; // Texture containing the scene
uniform sampler2D fireTexture; // Fire texture
uniform vec3 fireLightColor; // Fire light color
uniform vec3 mirrorNormal;
uniform vec3 mirrorPoint;

void main() {
    vec3 I = normalize(TexCoords.xy - 0.5);
    vec3 N = normalize(mirrorNormal);

    // Reflect the vector
    vec3 R = reflect(I, N);

    // Calculate texture coordinates based on reflection
    vec2 reflectTexCoords = R.xy + 0.5;

    // Sample the color from the screen texture
    vec4 reflectColor = texture(screenTexture, reflectTexCoords);

    // Sample the fire texture
    vec4 fireColor = texture(fireTexture, TexCoords);

    // Combine the reflected color and fire color
    vec4 combinedColor = reflectColor + vec4(fireLightColor, 1.0) * fireColor;

    // Output the combined color
    FragColor = combinedColor;
}