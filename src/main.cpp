#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#include "cs488.h"

#define CUBE_GLASS
//#define TEXTURE

// Error callback for GLFW
void glfwErrorCallback(int error, const char* description) {
    std::cerr << "GLFW Error: " << description << std::endl;
}

// Window size callback for GLFW
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

// Function to initialize OpenGL
GLFWwindow* initializeOpenGL() {
    glfwSetErrorCallback(glfwErrorCallback);

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return nullptr;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(800, 600, "Particle System", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return nullptr;
    }

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return nullptr;
    }

    glViewport(0, 0, 800, 600);
    return window;
}

// Function to read a file and return its contents as a string
std::string readFile(const char* filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

// Function to load, compile, and link shaders
GLuint loadShader(const char* vertexPath, const char* fragmentPath) {
    std::string vertexCode = readFile(vertexPath);
    const char* vertexSource = vertexCode.c_str();
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    glCompileShader(vertexShader);

    GLint success;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cerr << "Error: Vertex Shader compilation failed\n" << infoLog << std::endl;
    }

    std::string fragmentCode = readFile(fragmentPath);
    const char* fragmentSource = fragmentCode.c_str();
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cerr << "Error: Fragment Shader compilation failed\n" << infoLog << std::endl;
    }

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cerr << "Error: Shader Program linking failed\n" << infoLog << std::endl;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return shaderProgram;
}

GLuint cubeVAO, cubeVBO, cubeEBO;
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 2.0f);

void setupCube(float pos) {
    float cubeVertices[] = {
            // positions // normals
            -0.5f, -0.5f, -0.5f, 0.0f, 0.0f, -1.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 0.0f, -1.0f,
            0.5f, 0.5f, -0.5f, 0.0f, 0.0f, -1.0f,
            -0.5f, 0.5f, -0.5f, 0.0f, 0.0f, -1.0f,
            -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 1.0f,
            -0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 1.0f,
            -0.5f, 0.5f, 0.5f, -1.0f, 0.0f, 0.0f,
            -0.5f, 0.5f, -0.5f, -1.0f, 0.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, -1.0f, 0.0f, 0.0f,
            -0.5f, -0.5f, 0.5f, -1.0f, 0.0f, 0.0f,
            0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 0.0f,
            0.5f, 0.5f, -0.5f, 1.0f, 0.0f, 0.0f,
            0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f,
            0.5f, -0.5f, 0.5f, 1.0f, 0.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 0.0f, -1.0f, 0.0f,
            0.5f, -0.5f, -0.5f, 0.0f, -1.0f, 0.0f,
            0.5f, -0.5f, 0.5f, 0.0f, -1.0f, 0.0f,
            -0.5f, -0.5f, 0.5f, 0.0f, -1.0f, 0.0f,
            -0.5f, 0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f, 0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f,
            -0.5f, 0.5f, 0.5f, 0.0f, 1.0f, 0.0f
    };

    for (int i = 0; i < 24; i++) {
        cubeVertices[i * 6] *= 0.3f;
        cubeVertices[i * 6 + 1] *= 0.3f;
        cubeVertices[i * 6 + 2] *= 0.3f;
        cubeVertices[i * 6] += pos;
        cubeVertices[i * 6 + 1] += pos;
        cubeVertices[i * 6 + 2] += pos;
    }

    unsigned int cubeIndices[] = {
            0, 1, 2, 2, 3, 0,
            4, 5, 6, 6, 7, 4,
            8, 9, 10, 10, 11, 8,
            12, 13, 14, 14, 15, 12,
            16, 17, 18, 18, 19, 16,
            20, 21, 22, 22, 23, 20
    };

    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glGenBuffers(1, &cubeEBO);

    glBindVertexArray(cubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cubeIndices), cubeIndices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);
}

#ifdef CUBE_GLASS
void renderCube(GLuint cubeShaderProgram, const glm::mat4& projection, const glm::mat4& view, const glm::vec3& lightPos, const glm::vec3& lightColor) {
    glUseProgram(cubeShaderProgram);

    GLint projectionLoc = glGetUniformLocation(cubeShaderProgram, "projection");
    GLint viewLoc = glGetUniformLocation(cubeShaderProgram, "view");
    GLint modelLoc = glGetUniformLocation(cubeShaderProgram, "model");
    GLint lightPosLoc = glGetUniformLocation(cubeShaderProgram, "lightPos");
    GLint lightColorLoc = glGetUniformLocation(cubeShaderProgram, "lightColor");
    GLint viewPosLoc = glGetUniformLocation(cubeShaderProgram, "viewPos");
    GLint refractiveIndexLoc = glGetUniformLocation(cubeShaderProgram, "refractiveIndex");

    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, glm::vec3(-0.5f, -0.5f, -2.0f)); // Move the cube to a position to see multiple faces
    model = glm::rotate(model, glm::radians(45.0f), glm::vec3(1.0f, 1.0f, 0.0f)); // Rotate to see multiple faces

    glm::vec3 viewPos = glm::vec3(0.0f, 0.0f, 2.0f); // Position of the camera

    glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
    glUniform3fv(lightPosLoc, 1, glm::value_ptr(lightPos));
    glUniform3fv(lightColorLoc, 1, glm::value_ptr(lightColor));
    glUniform3fv(viewPosLoc, 1, glm::value_ptr(viewPos));
    glUniform1f(refractiveIndexLoc, 1.52f); // Refractive index for glass

    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glUseProgram(0);
}

#else

void renderCube(GLuint cubeShaderProgram, const glm::mat4& projection, const glm::mat4& view, const glm::vec3& lightPos, const glm::vec3& lightColor) {
    glUseProgram(cubeShaderProgram);

    GLint projectionLoc = glGetUniformLocation(cubeShaderProgram, "projection");
    GLint viewLoc = glGetUniformLocation(cubeShaderProgram, "view");
    GLint modelLoc = glGetUniformLocation(cubeShaderProgram, "model");
    GLint lightPosLoc = glGetUniformLocation(cubeShaderProgram, "lightPos");
    GLint lightColorLoc = glGetUniformLocation(cubeShaderProgram, "lightColor");
    GLint viewPosLoc = glGetUniformLocation(cubeShaderProgram, "viewPos");

    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, glm::vec3(-0.5f, -0.5f, -2.0f)); // Move the cube to a position to see multiple faces
    model = glm::rotate(model, glm::radians(45.0f), glm::vec3(1.0f, 1.0f, 0.0f)); // Rotate to see multiple faces

    glm::vec3 viewPos = glm::vec3(0.0f, 0.0f, 2.0f); // Position of the camera

    glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
    glUniform3fv(lightPosLoc, 1, glm::value_ptr(lightPos));
    glUniform3fv(lightColorLoc, 1, glm::value_ptr(lightColor));
    glUniform3fv(viewPosLoc, 1, glm::value_ptr(viewPos));

    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glUseProgram(0);
}
#endif



void renderMirror(GLuint shaderProgram, GLuint mirrorVAO, glm::mat4& projection, const glm::mat4& view, GLuint fireTexture, glm::vec3 fireLightColor) {
    glUseProgram(shaderProgram);

    GLint projectionLoc = glGetUniformLocation(shaderProgram, "projection");
    GLint viewLoc = glGetUniformLocation(shaderProgram, "view");
    GLint modelLoc = glGetUniformLocation(shaderProgram, "model");
    GLint screenTextureLoc = glGetUniformLocation(shaderProgram, "screenTexture");
    GLint fireTextureLoc = glGetUniformLocation(shaderProgram, "fireTexture");
    GLint fireLightColorLoc = glGetUniformLocation(shaderProgram, "fireLightColor");
    GLint mirrorNormalLoc = glGetUniformLocation(shaderProgram, "mirrorNormal");
    GLint mirrorPointLoc = glGetUniformLocation(shaderProgram, "mirrorPoint");

    if (projectionLoc == -1) std::cerr << "Failed to get uniform location for projection." << std::endl;
    if (viewLoc == -1) std::cerr << "Failed to get uniform location for view." << std::endl;
    if (modelLoc == -1) std::cerr << "Failed to get uniform location for model." << std::endl;
    if (screenTextureLoc == -1) std::cerr << "Failed to get uniform location for screenTexture." << std::endl;
    if (fireTextureLoc == -1) std::cerr << "Failed to get uniform location for fireTexture." << std::endl;
    if (fireLightColorLoc == -1) std::cerr << "Failed to get uniform location for fireLightColor." << std::endl;
    if (mirrorNormalLoc == -1) std::cerr << "Failed to get uniform location for mirrorNormal." << std::endl;
    if (mirrorPointLoc == -1) std::cerr << "Failed to get uniform location for mirrorPoint." << std::endl;


    if (projectionLoc == -1 || viewLoc == -1 || modelLoc == -1 || screenTextureLoc == -1 || fireTextureLoc == -1 || fireLightColorLoc == -1 || mirrorNormalLoc == -1 || mirrorPointLoc == -1) {
        std::cerr << "Failed to get uniform location for mirror shader." << std::endl;
        return;
    }

    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, glm::vec3(-1.0f, 0.0f, 0.0f)); // Move the mirror to the left
    model = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f)); // Rotate to face the scene

    glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, fireTexture); //
    glUniform1i(screenTextureLoc, 0);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, fireTexture);
    glUniform1i(fireTextureLoc, 1);

    glUniform3fv(fireLightColorLoc, 1, glm::value_ptr(fireLightColor));
    glUniform3fv(mirrorNormalLoc, 1, glm::value_ptr(glm::vec3(1.0f, 0.0f, 0.0f))); // Example normal pointing to the left
    glUniform3fv(mirrorPointLoc, 1, glm::value_ptr(glm::vec3(-1.5f, 0.0f, 0.0f))); // Example point on the mirror plane

    glBindVertexArray(mirrorVAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glUseProgram(0);
}

void setupMirrorBuffers(GLuint &mirrorVAO, GLuint &mirrorVBO, GLuint &mirrorEBO) {
    float mirrorVertices[] = {
            // positions // texture coords
            -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
            1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
            1.0f, 1.0f, 0.0f, 1.0f, 1.0f,
            -1.0f, 1.0f, 0.0f, 0.0f, 1.0f
    };

    unsigned int mirrorIndices[] = {
            0, 1, 2,
            2, 3, 0
    };

    glGenVertexArrays(1, &mirrorVAO);
    glGenBuffers(1, &mirrorVBO);
    glGenBuffers(1, &mirrorEBO);

    glBindVertexArray(mirrorVAO);

    glBindBuffer(GL_ARRAY_BUFFER, mirrorVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(mirrorVertices), mirrorVertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mirrorEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(mirrorIndices), mirrorIndices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);
}


int main() {
    std::cout << "Initializing OpenGL..." << std::endl;
    GLFWwindow* window = initializeOpenGL();
    if (!window) {
        return -1;
    }

    glfwSetKeyCallback(window, keyFunc);
    std::cout << "Loading shaders..." << std::endl;

//    GLuint mirrorShaderProgram = loadShader("../shaders/mirror_vertex.glsl", "../shaders/mirror_fragment.glsl");

#ifdef CUBE_GLASS
    GLuint cubeShaderProgram = loadShader("../shaders/cube_vertex.glsl", "../shaders/cube_glass_fragment.glsl");
#else
    GLuint cubeShaderProgram = loadShader("../shaders/cube_vertex.glsl", "../shaders/cube_fragment.glsl");
#endif

    std::cout << "Shader Programs loaded" << std::endl;

    std::cout << "Initializing particle system..." << std::endl;
    ParticleSystem particleSystem(4000); // Assuming 7000 particles

#ifdef TEXTURE
    GLuint particleShaderProgram = loadShader("../shaders/vertex.glsl", "../shaders/fragment_texture.glsl");
    particleSystem.loadTextures();
#else
    GLuint particleShaderProgram = loadShader("../shaders/vertex.glsl", "../shaders/fragment.glsl");
#endif

    std::cout << "Particle System initialized with 7000 particles" << std::endl;

    setupCube(0.5f);

    glm::mat4 projection = glm::mat4(
            1.81066f, 0.0f, 0.0f, 0.0f,
            0.0f, 2.41421f, 0.0f, 0.0f,
            0.0f, 0.0f, -1.002f, -1.0f,
            0.0f, 0.0f, -0.2002f, 0.0f
    );

    glm::mat4 view = glm::lookAt(glm::vec3(0.0f, 0.0f, 2.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 invProjection = glm::inverse(projection);
    glm::mat4 invView = glm::inverse(view);

//    // Ensure the shader program is valid before proceeding
//    if (particleShaderProgram == 0 || cubeShaderProgram == 0 || mirrorShaderProgram == 0) {
//        std::cerr << "Shader program is not valid." << std::endl;
//        return -1;
//    }

    // Check if uniform locations are valid
    GLint projectionLoc = glGetUniformLocation(particleShaderProgram, "projection");
    GLint viewLoc = glGetUniformLocation(particleShaderProgram, "view");
    GLint iTimeLoc = glGetUniformLocation(particleShaderProgram, "iTime");
//    GLint iResolutionLoc = glGetUniformLocation(particleShaderProgram, "iResolution");

    if (projectionLoc == -1) std::cerr << "Failed to get uniform location for projection." << std::endl;
    if (viewLoc == -1) std::cerr << "Failed to get uniform location for view." << std::endl;
    if (iTimeLoc == -1) std::cerr << "Failed to get uniform location for iTime." << std::endl;
//    if (iResolutionLoc == -1) std::cerr << "Failed to get uniform location for iResolution." << std::endl;

    if (projectionLoc == -1 || viewLoc == -1 || iTimeLoc == -1) {
        return -1;
    }

    std::cout << "Uniform locations obtained." << std::endl;

//    GLuint mirrorVAO, mirrorVBO, mirrorEBO;
//    setupMirrorBuffers(mirrorVAO, mirrorVBO, mirrorEBO);

    // Set clear color to black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    glm::vec3 fireLightPos = glm::vec3(0.0f, 0.5f, -2.0f); // Example fire light position
    glm::vec3 fireLightColor = glm::vec3(1.0f, 0.5f, 0.0f); // Example fire light color

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Update and render particles
        particleSystem.updateParticles(0.016f); // Update particles with a fixed time step
//        particleSystem.renderParticles(particleShaderProgram, mirrorShaderProgram, projection, view);
        particleSystem.renderParticles(particleShaderProgram, cubeShaderProgram, projection, view);

        // Render cube
//        renderCube(cubeShaderProgram, projection, view, fireLightPos, fireLightColor);
#ifdef CUBE_GLASS
        // Render cube
//        renderCube(cubeShaderProgram, projection, view, fireLightPos, fireLightColor);
#endif
//        renderMirror(mirrorShaderProgram, mirrorVAO, projection, view, particleSystem.fireTexture, fireLightColor);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
