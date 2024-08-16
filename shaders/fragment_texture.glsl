#version 330 core
in vec4 particleColor; // Assumed to come from vertex shader
in vec2 TexCoord;
out vec4 FragColor;

uniform sampler2D fireTexture;
uniform sampler2D smokeTexture;
uniform bool isSmoke; // Determine if the particle is smoke or fire
uniform sampler2D noiseTex;

uniform vec3 cameraPos;
uniform mat4 invProjection;
uniform mat4 invView;
uniform float iTime;
uniform vec3 fireLightColor;
uniform vec3 fireLightPos;

const int ITR = 256;
const float MAX = 10.0;
const float SRF = 0.01;

vec3 pal(vec3 a, vec3 b, vec3 c, vec3 d, float t)
{
    return a + b * cos(c * (t + d));
}

float map(vec3 sp)
{
    sp.x = mod(sp.x + 0.2, 2.5) - 1.25;
    sp *= 0.7;
    sp.y += 3.5;
    sp.y *= 0.25;
    return length(sp) - 1.0;
}

float mrch(vec3 ro, vec3 rd)
{
    float d0 = 0.0;
    for (int i = 0; i < ITR; ++i)
    {
        vec3 sp = ro + rd * d0;
        float ds = map(sp);
        d0 += ds;
        if (d0 > MAX || abs(ds) < SRF) break;
    }
    return d0;
}

vec3 dir(vec2 uv, vec3 r0, vec3 fx)
{
    vec3 w = normalize(fx - r0);
    vec3 u = normalize(cross(w, vec3(0.0, 1.0, 0.0)));
    vec3 v = normalize(cross(u, w));
    return mat3(u, v, w) * normalize(vec3(uv, 2.5));
}

vec3 nml(vec3 p)
{
    vec2 d = vec2(0.001, 0.0);
    return normalize(map(p) - vec3(map(p - d.xyy), map(p - d.yxy), map(p - d.yyx)));
}

float lighting(vec3 sp, vec3 lp)
{
    vec3 l = normalize(lp - sp);
    vec3 n = nml(sp);
    float df = clamp(dot(n, l), 0.0, 1.0);
    float ds = mrch(sp + n * 0.02, l);
    if (ds < length(lp - sp))
    {
        df *= 0.1;
    }
    return df;
}

vec4 sampleTexture() {
    vec2 texCoords = TexCoord; // Use TexCoord from vertex shader
    if (isSmoke) {
        return texture(smokeTexture, texCoords);
    } else {
        return texture(fireTexture, texCoords);
    }
}

void main()
{
    vec3 r0 = vec3(0.0, 0.5, -5.0 - (cos(iTime) * 0.5 + 0.5) * 0.1);
    vec3 fx = vec3(0.0, cos(iTime) * 0.5 + 0.5, sin(iTime));

    vec3 rd = normalize(fx - r0);

    // Sample texture based on the particle's state
    vec4 texColor = sampleTexture();

    // Modulate the texture color with the dynamic light effects
    texColor.rgb *= lighting(r0 + rd * 0.2, fireLightPos) * fireLightColor;

    // Apply noise texture to add dynamic effects
    vec2 noiseCoords = gl_FragCoord.xy / vec2(textureSize(noiseTex, 0));
    float noiseValue = texture(noiseTex, noiseCoords).r;

    // Blend the computed texture color with the particle's color
    vec4 finalColor = mix(texColor, particleColor, 0.5) * noiseValue; // Simple 50% blend modulated by noise

    FragColor = vec4(finalColor.rgb, 1.0); // Ensure full opacity
}