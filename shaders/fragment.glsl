#version 330 core
in vec4 particleColor; // Assumed to come from vertex shader
out vec4 FragColor;

uniform vec3 cameraPos;
uniform mat4 invProjection;
uniform mat4 invView;
uniform sampler3D volumeTexture;
uniform float iTime;
uniform vec2 iResolution;
uniform vec3 fireLightColor;
uniform vec3 fireLightPos;

const int ITR = 256;
const float MAX = 10.0;
const float SRF = 0.01;
const float airRefractiveIndex = 1.0;

vec3 pal(vec3 a, vec3 b, vec3 c, vec3 d, float t)
{
    return a + b * cos(c * (t + d));
}

float fresnel(vec3 I, vec3 N, float refractiveIndex) {
    float cosi = clamp(dot(I, N), -1.0, 1.0);
    float etai = airRefractiveIndex, etat = refractiveIndex;
    if (cosi > 0.0) { float tmp = etai; etai = etat; etat = tmp; }
    float sint = etai / etat * sqrt(max(0.0, 1.0 - cosi * cosi));
    if (sint >= 1.0) {
        return 1.0;
    } else {
        float cost = sqrt(max(0.0, 1.0 - sint * sint));
        cosi = abs(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) / 2.0;
    }
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

vec4 marchRay(vec3 rayOrigin, vec3 rayDir, float tMin, float tMax)
{
    float d0 = mrch(rayOrigin, rayDir);
    vec3 cl = vec3(0.0);

    if (d0 > 0.0)
    {
        vec3 hp = rayOrigin + rayDir * d0;
        vec3 a = vec3(0.2, 0.2, 0.2);
        vec3 b = vec3(0.5, 0.25, 0.23);
        vec3 c = vec3(7.0, 1.1, 0.76);
        vec3 d = vec3(1.0, 1.0, 1.75);
        cl += lighting(hp, fireLightPos) * pal(a, b, c, d, d0);
    }

    cl = pow(cl, vec3(3.0)) * 22.0;
    return vec4(cl, 1.0);
}

void main()
{
    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    uv = (uv - 0.5) * 2.0;
    uv.x *= iResolution.x / iResolution.y;

    vec3 r0 = vec3(0.0, 0.5, -5.0 - (cos(iTime) * 0.5 + 0.5) * 0.1);
    vec3 fx = vec3(0.0, cos(iTime) * 0.5 + 0.5, sin(iTime));
    vec3 rd = dir(uv, r0, fx);

    vec4 rayColor = marchRay(r0, rd, 0.0, 1.0);

    // Calculate blending factor based on particle's lifetime
    float blendFactor = particleColor.a * 2; // Assuming alpha encodes the lifetime progress from 1.0 (new) to 0.0 (dead)

    // Blend the particle color with the ray marched color
    vec4 blendedColor = mix(rayColor, particleColor, blendFactor);

    // Apply the fire light color dynamically
    vec3 fireLightEffect = fireLightColor * blendedColor.rgb;

    FragColor = vec4(fireLightEffect, blendedColor.a);
}
