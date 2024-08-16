#version 330 core
in vec3 FragPos;
in vec3 Normal;

out vec4 FragColor;

uniform vec3 lightPos;
uniform vec3 lightColor;
uniform vec3 viewPos;
uniform float refractiveIndex;

vec3 calculateRefraction(vec3 incident, vec3 normal, float eta) {
    float cosi = clamp(dot(incident, normal), -1.0, 1.0);
    float etai = 1.0, etat = eta;
    vec3 n = normal;
    if (cosi < 0.0) { cosi = -cosi; }
    else {
        float tmp = etai;
        etai = etat;
        etat = tmp;
        n = -normal;
    }
    float etaRatio = etai / etat;
    float k = 1.0 - etaRatio * etaRatio * (1.0 - cosi * cosi);
    return k < 0.0 ? vec3(0.0) : etaRatio * incident + (etaRatio * cosi - sqrt(k)) * n;
}

vec3 calculateReflection(vec3 incident, vec3 normal) {
    return incident - 2.0 * dot(incident, normal) * normal;
}

void main()
{
    // Ambient
    float ambientStrength = 0.1;
    vec3 ambient = ambientStrength * lightColor;

    // Diffuse
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // Specular
    float specularStrength = 0.5;
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = calculateReflection(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = specularStrength * spec * lightColor;

    // Reflection and Refraction
    vec3 I = normalize(FragPos - viewPos);
    vec3 R = calculateReflection(I, norm);
    vec3 refractColor = calculateRefraction(I, norm, refractiveIndex);

    vec3 color = mix(refractColor, R, 0.5) + ambient + diffuse + specular;
    FragColor = vec4(color, 1.0);
}
