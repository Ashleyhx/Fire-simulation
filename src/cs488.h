// =======================================
// CS488/688 base code
// (written by Toshiya Hachisuka)
// =======================================
#pragma once
#define _CRT_SECURE_NO_WARNINGS
#define NOMINMAX


// OpenGL
//#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>


// image loader and writer
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"


// linear algebra 
#include "linalg.h"
using namespace linalg::aliases;


// animated GIF writer
#include "gif.h"


// misc
#include <iostream>
#include <vector>
#include <cfloat>

#include <cmath>
#include <random>
#include <fstream>
#include <sstream>


#include <cstdint>
#include <algorithm>
#include <thread>
#include <mutex>
#include <glm/fwd.hpp>
#include <queue>


//#define FUEL_REDUCTION
//#define TEXTURE
#define SMOKE
#define WIND
//#define CUBE

// window size and resolution
// (do not make it too large - will be slow!)
constexpr int globalWidth = 512;
constexpr int globalHeight = 384;


static bool globalRecording = false;
static GifWriter globalGIFfile;
constexpr int globalGIFdelay = 1;


const float Tignition = 800.0f;
const float Tmax = 1500.0f;



namespace PCG32 {
    static uint64_t mcg_state = 0xcafef00dd15ea5e5u;	// must be odd
    static uint64_t const multiplier = 6364136223846793005u;
    uint32_t pcg32_fast(void) {
        uint64_t x = mcg_state;
        const unsigned count = (unsigned)(x >> 61);
        mcg_state = x * multiplier;
        x ^= x >> 22;
        return (uint32_t)(x >> (22 + count));
    }
    float rand() {
        return float(double(pcg32_fast()) / 4294967296.0);
    }
}


void keyFunc(GLFWwindow* window, int key, int scancode, int action, int mods) {
    printf("111");
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        switch (key) {


            case GLFW_KEY_F: {
                printf("F key pressed\n");
                if (!globalRecording) {
                    char fileName[1024];
                    sprintf(fileName, "output%d.gif", int(1000.0 * PCG32::rand()));
                    printf("Saving \"%s\"...\n", fileName);
                    GifBegin(&globalGIFfile, fileName, globalWidth, globalHeight, globalGIFdelay);
                    globalRecording = true;
                    printf("(Recording started)\n");
                } else {
                    GifEnd(&globalGIFfile);
                    globalRecording = false;
                    printf("(Recording done)\n");
                }
                break;
            }
        }
    }
}

static float randomFloat() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    return dist(gen);
}


struct Color {
    float r, g, b, a;
    Color(float r = 1.0f, float g = 1.0f, float b = 1.0f, float a = 1.0f)
            : r(r), g(g), b(b), a(a) {}

    static Color lerp(const Color& start, const Color& end, float t) {
        return Color(
                start.r + t * (end.r - start.r),
                start.g + t * (end.g - start.g),
                start.b + t * (end.b - start.b),
                start.a + t * (end.a - start.a)
        );
    }
};


class Particle {
public:
    glm::vec3 position;
    glm::vec3 velocity;
    float temperature;
    float fuel;
    Color startColor;
    Color endColor;
    float lifetime;
    float maxLifetime;
    Color currentColor;
    float size;
    float reactionCoordinate;
    float emissionIntensity;
    bool isSmoke;
    bool onAnotherObject = false;
    glm::vec3 windForce = glm::vec3(2.0f, 0.0f, 0.0f);
    glm::vec3 vorticity;
    glm::vec3 confinementForce;

    Particle(const glm::vec3& firePosition = glm::vec3(0.0f)) {
        reset(false, firePosition);
    }

    void reset(bool smoke, const glm::vec3& firePosition = glm::vec3(0.0f)) {
        if (smoke) {
            if (randomFloat() <= 0.7f) {
                return;
            }
            // Initialize as smoke particle
            position = glm::vec3(firePosition.x, firePosition.y, firePosition.z );
            velocity = glm::vec3(randomFloat() * 0.01f, randomFloat() * 0.01f, randomFloat() * 0.01f);
            temperature = 300.0f; // Smoke temperature
            fuel = 0.0f;
            startColor = Color(0.2f, 0.2f, 0.2f, 0.8f); // Dark gray for smoke
            endColor = Color(0.1f, 0.1f, 0.1f, 0.2f); // Fade out to almost black
            maxLifetime = 2.0f + randomFloat() * 1.0f; // Shorter lifetime for smoke
            size = 0.03f;
            vorticity = glm::vec3(0.0f);
            confinementForce = glm::vec3(0.0f);
        } else {
            float horizontalRadius = 0.23f; // Wider left to right
            float verticalRadius = 0.07f; // Narrower bottom to top

            // Random angle within the oval
            float theta = 2.0f * M_PI * randomFloat();

            float x = horizontalRadius * cos(theta);
            float z = verticalRadius * sin(theta);
            float y = -0.2f;

            float x_vel = x > 0 ? -0.3f : 0.3f;
            float z_vel = z > 0 ? -0.5f : 0.5f;
            float cur_rand = randomFloat();
            cur_rand *= cur_rand;
            position = glm::vec3(x, z, y);
            velocity = glm::vec3((cur_rand + 0.5f) * x_vel, cur_rand + 0.8f, (cur_rand + 0.5f) * z_vel * 0.2f);

            temperature = Tignition; // Start at ignition temperature
            fuel = 1.0f;
            startColor = Color(1.0f, 1.0f, 0.0f, 1.0f);
            endColor = Color(1.0f, 0.0f, 0.0f, 0.0f);
            maxLifetime = cur_rand;
            size = 0.015f;
            vorticity = glm::vec3(0.0f);
            confinementForce = glm::vec3(0.0f);
        }

        lifetime = 0.0f;
        currentColor = startColor;

        // Initialize reaction coordinate Y
        reactionCoordinate = 1.0f;
        emissionIntensity = temperature / Tmax;
        isSmoke = smoke;
    }

    void step(float deltaTime, const std::vector<float>& phi, const std::vector<float>& uf, const std::vector<float>& vf, const std::vector<float>& wf, int gridSize, float h) {
        lifetime += deltaTime;
        if (lifetime >= maxLifetime) {
            reset(isSmoke);
        } else {
            if (!std::isnan(confinementForce.x) && !std::isnan(confinementForce.y) && !std::isnan(confinementForce.z)) {
                velocity += confinementForce * deltaTime;
            }

            applyForces(deltaTime, phi, uf, vf, wf, gridSize, h);
            position += velocity * deltaTime;
            updateColorAndAlpha(deltaTime);

            float convectionTerm = glm::dot(velocity, glm::vec3(phi[position.x], phi[position.y], phi[position.z])) * deltaTime;
            reactionCoordinate = std::max(0.0f, reactionCoordinate - deltaTime - convectionTerm);

            updateTemperature();
            emissionIntensity = temperature / Tmax;
            updateSize();

#ifdef SMOKE
            // Transition to smoke when fuel is depleted
            if (!isSmoke && reactionCoordinate < 0.05f) {
                isSmoke = true;
                reset(true, position); // Reset as smoke particle, starting from the current fire position
            }
#endif

#ifdef WIND
            applyWindForce(deltaTime);
#endif
        }
    }

    void stepOnFire(float deltaTime, const std::vector<float>& phi, const std::vector<float>& uf, const std::vector<float>& vf, const std::vector<float>& wf, int gridSize, float h) {
        lifetime += deltaTime;
        if (lifetime >= maxLifetime) {
            reset(isSmoke);
        } else {
            position += velocity * deltaTime;
            updateColorAndAlpha(deltaTime);
            float convectionTerm = glm::dot(velocity, glm::vec3(phi[position.x], phi[position.y], phi[position.z])) * deltaTime;
            reactionCoordinate = std::max(0.0f, reactionCoordinate - deltaTime - convectionTerm);
            updateSize();
#ifdef SMOKE
            if (!isSmoke && reactionCoordinate < 0.05f) {
                isSmoke = true;
                reset(true, position); // Reset as smoke particle, starting from the current fire position
            }
#endif

#ifdef WIND
            applyWindForce(deltaTime);
#endif
            velocity += glm::vec3(0.01f, 0.0f, 0.0f);
        }
    }

    void updateTemperature() {
        if (reactionCoordinate >= 0.9f) {
            temperature = Tignition + (Tmax - Tignition) * (1.0f - reactionCoordinate);
        } else {
            temperature = Tmax;
        }
    }

    void updateColorAndAlpha(float deltaTime) {
        float t = lifetime / maxLifetime;
        if (!isSmoke) {
            currentColor = Color::lerp(startColor, endColor, t);
            currentColor.a = 1.0f - t; // Fade out over lifetime
        } else {
            currentColor = Color::lerp(startColor, endColor, t);
//            currentColor.a = 1.0f - t; // Fade out over lifetime
        }
    }

    void updateSize() {
        // Reduce size as fuel burns out
        size = std::min(0.015f, reactionCoordinate * 0.015f + 0.007f);
    }

    void applyForces(float deltaTime, const std::vector<float>& phi, const std::vector<float>& uf, const std::vector<float>& vf, const std::vector<float>& wf, int gridSize, float h) {
        int i = static_cast<int>((position.x + 1.0f) * 0.5f * (gridSize - 1));
        int j = static_cast<int>((position.y + 1.0f) * 0.5f * (gridSize - 1));
        int k = static_cast<int>((position.z + 1.0f) * 0.5f * (gridSize - 1));

        if (i >= 0 && i < gridSize && j >= 0 && j < gridSize && k >= 0 && k < gridSize) {
            int index = i + gridSize * (j + gridSize * k);

            // Get the velocity field at the particle's position
            float u = (uf[index] + uf[index + 1]) * 0.5f;
            float v = (vf[index] + vf[index + gridSize]) * 0.5f;
            float w = (wf[index] + wf[index + gridSize * gridSize]) * 0.5f;
            glm::vec3 velocityField = glm::vec3(u, v, w);

            // Calculate the gradient (normal) of phi at the particle's position
            float phiX = (phi[index + 1] - phi[index - 1]) / (2 * h);
            float phiY = (phi[index + gridSize] - phi[index - gridSize]) / (2 * h);
            float phiZ = (phi[index + gridSize * gridSize] - phi[index - gridSize * gridSize]) / (2 * h);
            glm::vec3 normal = glm::normalize(glm::vec3(phiX, phiY, phiZ));

            velocity += velocityField * deltaTime;
            if (isSmoke) {
                velocity.y += deltaTime * 0.3f; // Smoke rises
            } else {
                velocity.y -= 9.81f * deltaTime * 0.1f; // Fire affected by gravity
            }

            // Apply buoyancy force
            float Tair = 300.0f; // Ambient temperature
            float alpha = 0.5f; // Buoyancy constant
            glm::vec3 buoyancy = alpha * (temperature - Tair) * glm::vec3(0.0f, 1.0f, 0.0f);
//            velocity += buoyancy * deltaTime;
        }
    }

    void applyWindForce(float deltaTime) {
        velocity += windForce * deltaTime;
    }

    static Color blackbodyColor(float temperature) {
        // Implement the blackbody radiation color calculation based on temperature
        float t = glm::clamp((temperature - Tignition) / (Tmax - Tignition), 0.0f, 1.0f);
        return Color(t * 0.6f, t * 0.8f, t * 0.4f, 0.3f);
    }
};


struct ParticleInstance {
    glm::vec3 offset;
    float size;
    glm::vec4 color;
    glm::vec2 texCoords;
};


class ParticleSystem {
public:
    std::vector<Particle> particles;
    std::vector<Particle> cubeParticles; // Separate vector for particles on the cube
    GLuint quadVAO, quadVBO, instanceVBO, EBO;
    GLuint volumeTexture;
    const int volumeSize = 64; // Size of the volume texture
    const int gridSize = 64; // Grid size for level set
    float h; // Grid spacing
    std::vector<float> phi; // Level set grid
    std::vector<float> uf, vf, wf; // Velocity fields for fuel
    float fuelAmount; // Total fuel amount
    float initialFuelAmount; // Initial fuel amount for resetting
    GLuint fireTexture, smokeTexture;
    GLuint noiseTexture;

    std::vector<std::vector<std::vector<std::vector<Particle*>>>> grid;
    int gridResolution = 128;
    float gridCellSize;

    glm::vec3 cubePosition;
    glm::vec3 cubeVelocity = glm::vec3(0.2f, 0.0f, 0.0f);
    bool cubeOnFire;

    ParticleSystem(int numParticles, float initialFuel = 1000.0f) {
        particles.resize(numParticles);
        phi.resize(gridSize * gridSize * gridSize, 1.0f);
        uf.resize(gridSize * gridSize * gridSize, 0.0f);
        vf.resize(gridSize * gridSize * gridSize, 0.0f);
        wf.resize(gridSize * gridSize * gridSize, 0.0f);
        h = 2.0f / (gridSize - 1);
        initializeParticles();
        initializeGrid();
        setupOpenGLBuffers();
        setupVolumeTexture();

        initialFuelAmount = initialFuel;
        fuelAmount = initialFuel;
        noiseTexture = loadTexture("../shaders/img.png");
        cubePosition = glm::vec3(-1.5f, 0.0f, -1.5f);
        cubeVelocity = glm::vec3(0.2f, 0.0f, 0.0f);
        cubeOnFire = false;
    }

    void initializeParticles() {
        for (auto& particle : particles) {
            particle.reset(false); // Start as fire particles
        }
    }

    void initializeGrid() {
        gridCellSize = 2.0f / gridResolution;
        grid.resize(gridResolution, std::vector<std::vector<std::vector<Particle*>>>(gridResolution, std::vector<std::vector<Particle*>>(gridResolution)));
    }

    glm::ivec3 getGridIndex(const glm::vec3& position) {
        return glm::ivec3((position + 1.0f) * 0.5f * float (gridResolution));
    }

    void setupOpenGLBuffers() {
        float quadVertices[] = {
                -0.5f, -0.5f, 0.0f,
                0.5f, -0.5f, 0.0f,
                0.5f, 0.5f, 0.0f,
                -0.5f, 0.5f, 0.0f,
        };
        unsigned int quadIndices[] = {
                0, 1, 2, 2, 3, 0
        };

        glGenVertexArrays(1, &quadVAO);
        glGenBuffers(1, &quadVBO);
        glGenBuffers(1, &EBO);
        glGenBuffers(1, &instanceVBO);

        glBindVertexArray(quadVAO);

        glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIndices), quadIndices, GL_STATIC_DRAW);

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferData(GL_ARRAY_BUFFER, (particles.size() + cubeParticles.size()) * sizeof(ParticleInstance), nullptr, GL_DYNAMIC_DRAW);

        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleInstance), (void*)offsetof(ParticleInstance, offset));
        glVertexAttribDivisor(1, 1);

        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(ParticleInstance), (void*)offsetof(ParticleInstance, size));
        glVertexAttribDivisor(2, 1);

        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(ParticleInstance), (void*)offsetof(ParticleInstance, color));
        glVertexAttribDivisor(3, 1);

        glEnableVertexAttribArray(4);
        glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(ParticleInstance), (void*)offsetof(ParticleInstance, color));
        glVertexAttribDivisor(4, 1);

        glEnableVertexAttribArray(5);
        glVertexAttribPointer(5, 2, GL_FLOAT, GL_FALSE, sizeof(ParticleInstance), (void*)offsetof(ParticleInstance, texCoords));
        glVertexAttribDivisor(5, 1);

        glBindVertexArray(0);
    }

    void setupVolumeTexture() {
        glGenTextures(1, &volumeTexture);
        glBindTexture(GL_TEXTURE_3D, volumeTexture);
        glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, volumeSize, volumeSize, volumeSize, 0, GL_RED, GL_FLOAT, nullptr);

        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        glBindTexture(GL_TEXTURE_3D, 0);
    }

    void updateParticles(float deltaTime) {
        for (auto& particle : particles) {
            particle.step(deltaTime, phi, uf, vf, wf, gridSize, h);
        }

        for (auto& particle : cubeParticles) {
            particle.stepOnFire(deltaTime, phi, uf, vf, wf, gridSize, h);
            particle.position += cubeVelocity * deltaTime; // Move with the cube
        }

        handleCollisions();
        updateLevelSet(deltaTime);

        computeVorticity();
        computeConfinementForce();

        updateVolumeTexture();
        updateInstanceBuffer();

        cubePosition += cubeVelocity * deltaTime;
        checkCubeFireCollision();

        if (cubeOnFire) {
            addParticlesAtCubePosition();
        }

        // Decrease fuel amount and reduce particle count if fuel is low
#ifdef FUEL_REDUCTION
        fuelAmount -= deltaTime * 5.0f; // Decrease fuel over time
        if (fuelAmount <= 0.0f) {
            // If fuel is exhausted, gradually reduce the number of particles
            if (!particles.empty()) {
                particles.pop_back();
            }
        } else {
            // Adjust the number of particles based on remaining fuel
            int targetParticleCount = static_cast<int>((fuelAmount / initialFuelAmount) * particles.size());
            while (particles.size() > targetParticleCount && !particles.empty()) {
                particles.pop_back();
            }
        }
#endif
    }

    void checkCubeFireCollision() {
        if (cubePosition.x > 0.0f) {
            cubeOnFire = true;
        }
    }

    void addParticlesAtCubePosition() {
        glm::vec3 newPosition;
        float x;

        for (int i = 0; i < 30; ++i) {
            if (randomFloat() < 0.5f) {
                x = -1.0f;
            } else {
                x = 1.0f;
            }
            newPosition = cubePosition + glm::vec3(randomFloat() * 0.07f * x, 0.04f, randomFloat() * 0.1f);
            Particle newParticle;
            newParticle.reset(false, newPosition);
            newParticle.position = newPosition;
            newParticle.onAnotherObject = true;
            newParticle.size = 0.010f;
            newParticle.maxLifetime = randomFloat() * 0.4;
            particles.push_back(newParticle);
        }
    }

    void computeVorticity() {
        // Clear grid
        for (auto& plane : grid) {
            for (auto& row : plane) {
                for (auto& cell : row) {
                    cell.clear();
                }
            }
        }

        for (auto& particle : particles) {
            glm::ivec3 idx = getGridIndex(particle.position);
            if (idx.x >= 0 && idx.x < gridResolution && idx.y >= 0 && idx.y < gridResolution && idx.z >= 0 && idx.z < gridResolution) {
                grid[idx.x][idx.y][idx.z].push_back(&particle);
            }
        }

        // Compute vorticity for each particle
        for (auto& plane : grid) {
            for (auto& row : plane) {
                for (auto& cell : row) {
                    for (auto& particle : cell) {
                        glm::vec3 vorticity(0.0f);
                        for (auto& neighbor : cell) {
                            if (particle != neighbor) {
                                glm::vec3 r = neighbor->position - particle->position;
                                glm::vec3 v = neighbor->velocity - particle->velocity;
                                if (glm::length(r) < neighbor->size && glm::dot(r, r) != 0) {
                                    vorticity += glm::cross(v, r) / glm::dot(r, r);
                                }
                            }
                        }
                        particle->vorticity = vorticity;
                    }
                }
            }
        }
    }

    void computeConfinementForce() {
        float epsilon = 0.05f;
        for (auto& particle : particles) {
            glm::vec3 gradient(0.0f);
            glm::ivec3 idx = getGridIndex(particle.position);

            // Compute gradient of vorticity magnitude
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        glm::ivec3 neighborIdx = idx + glm::ivec3(dx, dy, dz);
                        if (neighborIdx.x >= 0 && neighborIdx.x < gridResolution && neighborIdx.y >= 0 && neighborIdx.y < gridResolution && neighborIdx.z >= 0 && neighborIdx.z < gridResolution) {
                            for (auto& neighbor : grid[neighborIdx.x][neighborIdx.y][neighborIdx.z]) {
                                glm::vec3 r = neighbor->position - particle.position;
                                gradient += (neighbor->vorticity - particle.vorticity) / (glm::dot(r, r) + 1e-5f); // Avoid division by zero
                            }
                        }
                    }
                }
            }

            glm::vec3 n = glm::normalize(gradient);
            glm::vec3 confinementForce = epsilon * glm::cross(n, particle.vorticity);

            if (!std::isnan(confinementForce.x) && !std::isnan(confinementForce.y) && !std::isnan(confinementForce.z)) {
                particle.confinementForce = confinementForce;
            }
        }
    }


    void handleCollisions() {
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                if (checkCollision(particles[i], particles[j])) {
                    resolveCollision(particles[i], particles[j]);
                }
            }
        }
    }

    static bool checkCollision(const Particle& a, const Particle& b) {
        float distance = glm::distance(a.position, b.position);
        return distance < (a.size + b.size);
    }

    static void resolveCollision(Particle& a, Particle& b) {
        glm::vec3 normal = glm::normalize(b.position - a.position);
        glm::vec3 relativeVelocity = b.velocity - a.velocity;
        float velocityAlongNormal = glm::dot(relativeVelocity, normal);

        if (velocityAlongNormal > 0) {
            return;
        }

        float restitution = 0.2f; // Coefficient of restitution (bounciness)
        float impulseMagnitude = -(1 + restitution) * velocityAlongNormal / 2.0f;
        glm::vec3 impulse = impulseMagnitude * normal;

        a.velocity -= impulse;
        b.velocity += impulse;
    }

    void updateInstanceBuffer() {
        std::vector<ParticleInstance> instances;
        instances.reserve(particles.size());
        for (auto& particle : particles) {
            ParticleInstance instance;
            instance.offset = particle.position;
            instance.size = particle.size;
            instance.color = glm::vec4(particle.currentColor.r, particle.currentColor.g, particle.currentColor.b, particle.currentColor.a);
            instance.texCoords = glm::vec2(0.0f, 0.0f); // Replace with actual texture coordinates if needed
            instances.push_back(instance);
        }

        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(ParticleInstance), instances.data(), GL_DYNAMIC_DRAW);
    }

    void updateVolumeTexture() {
        std::vector<GLubyte> volumeData(volumeSize * volumeSize * volumeSize, 0);

        // Update volume data based on particle positions and colors
        for (const auto& particle : particles) {
            int x = static_cast<int>((particle.position.x + 1.0f) * 0.5f * volumeSize);
            int y = static_cast<int>((particle.position.y + 1.0f) * 0.5f * volumeSize);
            int z = static_cast<int>((particle.position.z + 1.0f) * 0.5f * volumeSize);

            if (x >= 0 && x < volumeSize && y >= 0 && y < volumeSize && z >= 0 && z < volumeSize) {
                int index = x + volumeSize * (y + volumeSize * z);
                volumeData[index] = static_cast<GLubyte>(glm::clamp(particle.currentColor.r * 255.0f, 0.0f, 255.0f));
            }
        }

        glBindTexture(GL_TEXTURE_3D, volumeTexture);
        glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, volumeSize, volumeSize, volumeSize, GL_RED, GL_UNSIGNED_BYTE, volumeData.data());
        glBindTexture(GL_TEXTURE_3D, 0);
    }

    void updateLevelSet(float deltaTime) {
        // Calculate gradients and velocities
        std::vector<float> gradPhiX(gridSize * gridSize * gridSize, 0.0f);
        std::vector<float> gradPhiY(gridSize * gridSize * gridSize, 0.0f);
        std::vector<float> gradPhiZ(gridSize * gridSize * gridSize, 0.0f);

        for (int i = 1; i < gridSize - 1; ++i) {
            for (int j = 1; j < gridSize - 1; ++j) {
                for (int k = 1; k < gridSize - 1; ++k) {
                    int index = i + gridSize * (j + gridSize * k);
                    gradPhiX[index] = (phi[index + 1] - phi[index - 1]) / (2 * h);
                    gradPhiY[index] = (phi[index + gridSize] - phi[index - gridSize]) / (2 * h);
                    gradPhiZ[index] = (phi[index + gridSize * gridSize] - phi[index - gridSize * gridSize]) / (2 * h);
                }
            }
        }

        // Update level set φ using upwind differencing
        std::vector<float> phiNew = phi;

        for (int i = 1; i < gridSize - 1; ++i) {
            for (int j = 1; j < gridSize - 1; ++j) {
                for (int k = 1; k < gridSize - 1; ++k) {
                    int index = i + gridSize * (j + gridSize * k);
                    glm::vec3 w = glm::vec3(uf[index], vf[index], wf[index]) + glm::vec3(gradPhiX[index], gradPhiY[index], gradPhiZ[index]);
                    float phiX = (w.x > 0) ? (phi[index] - phi[index - 1]) / h : (phi[index + 1] - phi[index]) / h;
                    float phiY = (w.y > 0) ? (phi[index] - phi[index - gridSize]) / h : (phi[index + gridSize] - phi[index]) / h;
                    float phiZ = (w.z > 0) ? (phi[index] - phi[index - gridSize * gridSize]) / h : (phi[index + gridSize * gridSize] - phi[index]) / h;

                    phiNew[index] = phi[index] - deltaTime * (w.x * phiX + w.y * phiY + w.z * phiZ);
                }
            }
        }

        phi = phiNew;

        // Reinitialize φ to maintain it as a signed distance function
        reinitializePhi();
    }

    void reinitializePhi() {
        std::vector<float> phiNew = phi;

        std::priority_queue<std::tuple<float, int, int, int>> heap;
        for (int i = 1; i < gridSize - 1; ++i) {
            for (int j = 1; j < gridSize - 1; ++j) {
                for (int k = 1; k < gridSize - 1; ++k) {
                    int index = i + gridSize * (j + gridSize * k);
                    if (std::abs(phi[index]) < h) {
                        heap.push(std::make_tuple(std::abs(phi[index]), i, j, k));
                    }
                }
            }
        }

        while (!heap.empty()) {
            auto [dist, i, j, k] = heap.top();
            heap.pop();

            int index = i + gridSize * (j + gridSize * k);
            if (i > 0) updatePhiDistance(i - 1, j, k, index, phiNew, heap);
            if (i < gridSize - 1) updatePhiDistance(i + 1, j, k, index, phiNew, heap);
            if (j > 0) updatePhiDistance(i, j - 1, k, index, phiNew, heap);
            if (j < gridSize - 1) updatePhiDistance(i, j + 1, k, index, phiNew, heap);
            if (k > 0) updatePhiDistance(i, j, k - 1, index, phiNew, heap);
            if (k < gridSize - 1) updatePhiDistance(i, j, k + 1, index, phiNew, heap);
        }

        phi = phiNew;
    }

    void updatePhiDistance(int i, int j, int k, int centerIndex, std::vector<float>& phiNew, std::priority_queue<std::tuple<float, int, int, int>>& heap) const {
        int index = i + gridSize * (j + gridSize * k);
        if (phiNew[index] > 0) {
            float dist = std::min(phiNew[centerIndex] + h, phiNew[index]);
            if (dist < phiNew[index]) {
                phiNew[index] = dist;
                heap.push(std::make_tuple(dist, i, j, k));
            }
        }
    }

    void getLightData(std::vector<glm::vec3>& lightPositions, std::vector<glm::vec3>& lightColors) const {
        for (const auto& particle : particles) {
            lightPositions.push_back(particle.position);
            lightColors.push_back(glm::vec3(particle.currentColor.r, particle.currentColor.g, particle.currentColor.b));
        }
    }

    GLuint generateNoiseTexture() {
        const int texSize = 256;
        std::vector<unsigned char> data(texSize * texSize);

        for (int y = 0; y < texSize; ++y) {
            for (int x = 0; x < texSize; ++x) {
                data[y * texSize + x] = rand() % 256;
            }
        }

        GLuint textureID;
        glGenTextures(1, &textureID);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, texSize, texSize, 0, GL_RED, GL_UNSIGNED_BYTE, data.data());
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        return textureID;
    }

    GLuint loadTexture(const char* path) {
        int width, height, nrChannels;
        auto data = stbi_load("../shaders/s1.png", &width, &height, &nrChannels, 0);

        if (stbi_failure_reason()) {
            std::cerr << "Failed to load texture: " << stbi_failure_reason() << std::endl;
        }

        GLenum format = (nrChannels == 4) ? GL_RGBA : GL_RGB;

        GLuint textureID;
        glGenTextures(1, &textureID);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);

        glGenerateMipmap(GL_TEXTURE_2D);
        stbi_image_free(data);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        return textureID;
    }

    // Load textures
    void loadTextures() {
        fireTexture = loadTexture("../shaders/fire_particle.png");
        smokeTexture = loadTexture("../shaders/fire_particle.png");
        if (fireTexture == 0 || smokeTexture == 0) {
            std::cerr << "Failed to load textures" << std::endl;
        } else {
            std::cout << "Loaded textures" << std::endl;
        }
    }

    GLuint cubeVAO, cubeVBO, cubeEBO;
    glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 2.0f);

    void setupCube(float pos) {
        float cubeVertices[] = {
                // positions // normals
                -0.15f, -0.15f, -0.15f, 0.0f, 0.0f, -1.0f,
                0.15f, -0.15f, -0.15f, 0.0f, 0.0f, -1.0f,
                0.15f, 0.15f, -0.15f, 0.0f, 0.0f, -1.0f,
                -0.15f, 0.15f, -0.15f, 0.0f, 0.0f, -1.0f,
                -0.15f, -0.15f, 0.15f, 0.0f, 0.0f, 1.0f,
                0.15f, -0.15f, 0.15f, 0.0f, 0.0f, 1.0f,
                0.15f, 0.15f, 0.15f, 0.0f, 0.0f, 1.0f,
                -0.15f, 0.15f, 0.15f, 0.0f, 0.0f, 1.0f,
                -0.15f, 0.15f, 0.15f, -1.0f, 0.0f, 0.0f,
                -0.15f, 0.15f, -0.15f, -1.0f, 0.0f, 0.0f,
                -0.15f, -0.15f, -0.15f, -1.0f, 0.0f, 0.0f,
                -0.15f, -0.15f, 0.15f, -1.0f, 0.0f, 0.0f,
                0.15f, 0.15f, 0.15f, 1.0f, 0.0f, 0.0f,
                0.15f, 0.15f, -0.15f, 1.0f, 0.0f, 0.0f,
                0.15f, -0.15f, -0.15f, 1.0f, 0.0f, 0.0f,
                0.15f, -0.15f, 0.15f, 1.0f, 0.0f, 0.0f,
                -0.15f, -0.15f, -0.15f, 0.0f, -1.0f, 0.0f,
                0.15f, -0.15f, -0.15f, 0.0f, -1.0f, 0.0f,
                0.15f, -0.15f, 0.15f, 0.0f, -1.0f, 0.0f,
                -0.15f, -0.15f, 0.15f, 0.0f, -1.0f, 0.0f,
                -0.15f, 0.15f, -0.15f, 0.0f, 1.0f, 0.0f,
                0.15f, 0.15f, -0.15f, 0.0f, 1.0f, 0.0f,
                0.15f, 0.15f, 0.15f, 0.0f, 1.0f, 0.0f,
                -0.15f, 0.15f, 0.15f, 0.0f, 1.0f, 0.0f
        };

        for (int i = 0; i < 24; i++) {
            cubeVertices[i * 6] *= 0.5f;
            cubeVertices[i * 6 + 1] *= 0.5f;
            cubeVertices[i * 6 + 2] *= 0.5f;
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

    void renderCube(GLuint cubeShaderProgram, const glm::mat4& projection, const glm::mat4& view, const glm::mat4& model) {
        glUseProgram(cubeShaderProgram);

        GLint projectionLoc = glGetUniformLocation(cubeShaderProgram, "projection");
        GLint viewLoc = glGetUniformLocation(cubeShaderProgram, "view");
        GLint modelLoc = glGetUniformLocation(cubeShaderProgram, "model");
        GLint lightPosLoc = glGetUniformLocation(cubeShaderProgram, "lightPos");
        GLint viewPosLoc = glGetUniformLocation(cubeShaderProgram, "viewPos");
        GLint colorLoc = glGetUniformLocation(cubeShaderProgram, "color");

        glm::vec3 lightPos = glm::vec3(0.0f, 1.0f, 2.0f);
        glm::vec3 viewPos = glm::vec3(0.0f, 0.0f, 2.0f);
        glm::vec3 color = glm::vec3(1.0f, 0.0f, 0.0f);

        glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        glUniform3fv(lightPosLoc, 1, glm::value_ptr(lightPos));
        glUniform3fv(viewPosLoc, 1, glm::value_ptr(viewPos));
        glUniform3fv(colorLoc, 1, glm::value_ptr(color));

        glBindVertexArray(cubeVAO);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        glUseProgram(0);
    }

    void renderParticles(GLuint shaderProgram, GLuint cubeShaderProgram, const glm::mat4& projection, const glm::mat4& view) {
#ifdef CUBE
        static float cubeTime = 0.0f;
        cubeTime += 0.01f;
        cubePosition = glm::vec3(-0.5f + cubeTime, 0.1f, 0.0f);
        glm::mat4 model = glm::translate(glm::mat4(1.0f), cubePosition);
        setupCube(0.0f);
        glm::vec3 fireLightPos = glm::vec3(0.0f, 0.5f, -2.0f);
        glm::vec3 fireLightColor2 = glm::vec3(1.0f, 0.5f, 0.0f);

        renderCube(cubeShaderProgram, projection, view, model);
#endif
        glUseProgram(shaderProgram);

        GLint projectionLoc = glGetUniformLocation(shaderProgram, "projection");
        GLint viewLoc = glGetUniformLocation(shaderProgram, "view");
        GLint fireLightColorLoc = glGetUniformLocation(shaderProgram, "fireLightColor");

        if (projectionLoc == -1) std::cerr << "Failed to get uniform location for projection." << std::endl;
        if (viewLoc == -1) std::cerr << "Failed to get uniform location for view." << std::endl;
        if (fireLightColorLoc == -1) std::cerr << "Failed to get uniform location for fireLightColor." << std::endl;

        glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 2.0f);
        glm::vec3 fireLightColor = glm::vec3(1.0f, 0.5f, 0.0f);

        glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
        glUniform3fv(fireLightColorLoc, 1, glm::value_ptr(fireLightColor));

#ifdef TEXTURE
        // Bind noise texture
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, noiseTexture);
//        glUniform1i(noiseTexLoc, 0);
#endif

        // Draw particles
        glBindVertexArray(quadVAO);
        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, particles.size());
        glBindVertexArray(0);

        glUseProgram(0);
    }
};



