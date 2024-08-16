#version 330 core
out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;

uniform vec3 cameraPos;
uniform vec3 lightPos[100];
uniform vec3 lightColor[100];
uniform int numLights;
uniform vec3 objectColor;

void main()
{
    vec3 ambient = vec3(0.0);
    vec3 diffuse = vec3(0.0);
    vec3 specular = vec3(0.0);
    vec3 viewDir = normalize(cameraPos - FragPos);
    vec3 norm = normalize(Normal);

    for (int i = 0; i < numLights; ++i) {
        vec3 lightDir = normalize(lightPos[i] - FragPos);
        vec3 reflectDir = reflect(-lightDir, norm);

        float diff = max(dot(norm, lightDir), 0.0);
        diffuse += diff * lightColor[i];

        float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
        specular += spec * lightColor[i];

        ambient += 0.1 * lightColor[i];
    }

    vec3 result = ambient + diffuse + specular;
    result *= objectColor;
    FragColor = vec4(result, 1.0);
}
