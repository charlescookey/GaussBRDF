#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <memory.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <random>
#include <iomanip>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#define SQ(x) (x * x)


struct Ray {
    glm::vec3 o; // origin
    glm::vec3 dir; // direction

    Ray(const glm::vec3& origin, const glm::vec3& direction) : o(origin), dir(direction) {}

    glm::vec3 at(float t) const {
        return o + dir * t;
    }
};

struct Gaussian {
    glm::vec3 pos;
    glm::mat3 covariance;
    float amplitude;

    Gaussian(const glm::vec3& position, const glm::mat3& cov, float amp = 1.0f)
        : pos(position), covariance(cov), amplitude(amp) {}

    // Evaluate Gaussian at point p
    float evaluate(const glm::vec3& p) const {
        glm::vec3 diff = p - pos;
        glm::vec3 temp = covariance * diff;
        float exponent = -0.5f * glm::dot(temp ,diff);


        float det = glm::determinant(covariance);
        float norm = 1.0f / std::sqrt(std::pow(2.0f * M_PI, 3.0f) * det);

        return amplitude * norm * std::exp(exponent);
    }
};

// Your analytical solution (with bug fix)
float computeGaussianIntegral(Ray ray, const Gaussian& gaussian) {
    glm::vec3 r = gaussian.pos - ray.o;

    ray.dir = glm::normalize(ray.dir);
    glm::vec3 _a = gaussian.covariance * ray.dir;
    float a = glm::dot(_a , ray.dir);
    glm::vec3 _c = gaussian.covariance * r;
    float c = glm::dot(_c , r);
    float b = glm::dot(_c , ray.dir);

    float firstPart = std::exp((SQ(b) / (2 * a) - (c / 2)));
    float secondPart = std::sqrt(M_PI / (2 * a));

    return firstPart * secondPart; // Fixed: was firstPart2 * secondPart2
}

// Monte Carlo integration along the ray
float monteCarloIntegral(const Ray& ray, const Gaussian& gaussian,
    float tMin, float tMax, int numSamples) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(tMin, tMax);

    float sum = 0.0f;
    for (int i = 0; i < numSamples; i++) {
        float t = dis(gen);
        glm::vec3 point = ray.at(t);
        sum += gaussian.evaluate(point);
    }

    float average = sum / numSamples;
    float integral = average * (tMax - tMin);

    return integral;
}