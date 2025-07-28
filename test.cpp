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
        float exponent = -0.5f * glm::dot(temp, diff);
        float det = glm::determinant(covariance);
        float norm = 1.0f / std::sqrt(std::pow(2.0f * M_PI, 3.0f) * det);
        return amplitude * norm * std::exp(exponent);
    }
};

// Corrected analytical solution - the original formula appears to be missing a factor of 2
float computeGaussianIntegral(Ray ray, const Gaussian& gaussian, bool debug = false) {
    glm::vec3 r = gaussian.pos - ray.o;
    ray.dir = glm::normalize(ray.dir);
    glm::vec3 _a = gaussian.covariance * ray.dir;
    float a = glm::dot(_a, ray.dir);
    glm::vec3 _c = gaussian.covariance * r;
    float c = glm::dot(_c, r);
    float b = glm::dot(_c, ray.dir);

    if (debug) {
        std::cout << "Debug analytical computation:" << std::endl;
        std::cout << "  r = (" << r.x << ", " << r.y << ", " << r.z << ")" << std::endl;
        std::cout << "  ray.dir = (" << ray.dir.x << ", " << ray.dir.y << ", " << ray.dir.z << ")" << std::endl;
        std::cout << "  a = " << a << std::endl;
        std::cout << "  b = " << b << std::endl;
        std::cout << "  c = " << c << std::endl;
    }

    float exponent = (SQ(b) / (2 * a) - (c / 2));
    float firstPart = std::exp(exponent);
    float secondPart = std::sqrt(M_PI / (2 * a));

    // Add the normalization factor and amplitude to match Monte Carlo
    float det = glm::determinant(gaussian.covariance);
    float norm = gaussian.amplitude / std::sqrt(std::pow(2.0f * M_PI, 3.0f) * det);

    if (debug) {
        std::cout << "  exponent = " << exponent << std::endl;
        std::cout << "  firstPart = " << firstPart << std::endl;
        std::cout << "  secondPart = " << secondPart << std::endl;
        std::cout << "  norm = " << norm << std::endl;
        std::cout << "  raw result = " << firstPart * secondPart << std::endl;
        std::cout << "  final result (with factor 2) = " << 2.0f * norm * firstPart * secondPart << std::endl;
    }

    // The original formula appears to be missing a factor of 2

    return 2.0f * norm * firstPart * secondPart;
}

// Original analytical solution for comparison
float computeGaussianIntegralOriginal(Ray ray, const Gaussian& gaussian) {
    glm::vec3 r = gaussian.pos - ray.o;
    ray.dir = glm::normalize(ray.dir);
    glm::vec3 _a = gaussian.covariance * ray.dir;
    float a = glm::dot(_a, ray.dir);
    glm::vec3 _c = gaussian.covariance * r;
    float c = glm::dot(_c, r);
    float b = glm::dot(_c, ray.dir);

    float firstPart = std::exp((SQ(b) / (2 * a) - (c / 2)));
    float secondPart = std::sqrt(M_PI / (2 * a));

    float det = glm::determinant(gaussian.covariance);
    float norm = gaussian.amplitude / std::sqrt(std::pow(2.0f * M_PI, 3.0f) * det);

    return norm * firstPart * secondPart;
}

// Alternative implementation: Direct numerical integration for comparison
float directNumericalIntegral(const Ray& ray, const Gaussian& gaussian,
    float tMin, float tMax, int numSamples) {
    float dt = (tMax - tMin) / numSamples;
    float sum = 0.0f;

    for (int i = 0; i < numSamples; i++) {
        float t = tMin + (i + 0.5f) * dt;  // midpoint rule
        glm::vec3 point = ray.at(t);
        sum += gaussian.evaluate(point);
    }

    return sum * dt;
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

// Debug function to check if the Gaussian makes sense
void debugGaussian(const Gaussian& gaussian) {
    float det = glm::determinant(gaussian.covariance);
    std::cout << "Gaussian debug:" << std::endl;
    std::cout << "  Position: (" << gaussian.pos.x << ", " << gaussian.pos.y << ", " << gaussian.pos.z << ")" << std::endl;
    std::cout << "  Amplitude: " << gaussian.amplitude << std::endl;
    std::cout << "  Determinant: " << det << std::endl;
    std::cout << "  Value at center: " << gaussian.evaluate(gaussian.pos) << std::endl;
    std::cout << std::endl;
}

int main2() {
    std::cout << std::fixed << std::setprecision(8);

    // Test case 1: Simple case with identity-like covariance
    std::cout << "=== Test Case 1: Simple Gaussian ===" << std::endl;
    Ray ray1(glm::vec3(0, 0, 0), glm::vec3(1, 0, 0));
    glm::mat3 cov1(1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 1.0f);
    Gaussian gaussian1(glm::vec3(2, 0, 0), cov1, 1.0f);

    debugGaussian(gaussian1);

    float analyticalOriginal1 = computeGaussianIntegralOriginal(ray1, gaussian1);
    float analyticalCorrected1 = computeGaussianIntegral(ray1, gaussian1, true);
    float monteCarlo1 = monteCarloIntegral(ray1, gaussian1, -10.0f, 15.0f, 1000000);
    float numerical1 = directNumericalIntegral(ray1, gaussian1, -10.0f, 15.0f, 100000);

    std::cout << "Original analytical result: " << analyticalOriginal1 << std::endl;
    std::cout << "Corrected analytical result: " << analyticalCorrected1 << std::endl;
    std::cout << "Monte Carlo result: " << monteCarlo1 << std::endl;
    std::cout << "Numerical integration: " << numerical1 << std::endl;
    std::cout << "Original vs Numerical error: " << std::abs(analyticalOriginal1 - numerical1) / numerical1 * 100 << "%" << std::endl;
    std::cout << "Corrected vs Numerical error: " << std::abs(analyticalCorrected1 - numerical1) / numerical1 * 100 << "%" << std::endl;

    // Test case 2: Off-axis ray
    std::cout << "\n=== Test Case 2: Off-axis Ray ===" << std::endl;
    Ray ray2(glm::vec3(0, 1, 0), glm::vec3(1, 0, 0));
    Gaussian gaussian2(glm::vec3(2, 1, 0), cov1, 1.0f);

    float analytical2 = computeGaussianIntegral(ray2, gaussian2);
    float monteCarlo2 = monteCarloIntegral(ray2, gaussian2, -10.0f, 15.0f, 1000000);
    float numerical2 = directNumericalIntegral(ray2, gaussian2, -10.0f, 15.0f, 100000);

    std::cout << "Analytical result: " << analytical2 << std::endl;
    std::cout << "Monte Carlo result: " << monteCarlo2 << std::endl;
    std::cout << "Numerical integration: " << numerical2 << std::endl;
    std::cout << "Analytical vs Numerical error: " << std::abs(analytical2 - numerical2) / numerical2 * 100 << "%" << std::endl;

    // Test case 3: Anisotropic Gaussian
    std::cout << "\n=== Test Case 3: Anisotropic Gaussian ===" << std::endl;
    Ray ray3(glm::vec3(-1, 0, 0), glm::vec3(1, 0, 0));
    glm::mat3 cov3(0.5f, 0.0f, 0.0f,
        0.0f, 2.0f, 0.0f,
        0.0f, 0.0f, 1.0f);
    Gaussian gaussian3(glm::vec3(1, 0, 0), cov3, 1.0f);

    float analytical3 = computeGaussianIntegral(ray3, gaussian3);
    float monteCarlo3 = monteCarloIntegral(ray3, gaussian3, -10.0f, 15.0f, 1000000);

    std::cout << "Analytical result: " << analytical3 << std::endl;
    std::cout << "Monte Carlo result: " << monteCarlo3 << std::endl;
    std::cout << "Relative error: " << std::abs(analytical3 - monteCarlo3) / std::max(analytical3, monteCarlo3) * 100 << "%" << std::endl;

    // Test case 4: Diagonal ray
    std::cout << "\n=== Test Case 4: Diagonal Ray ===" << std::endl;
    Ray ray4(glm::vec3(0, 0, 0), glm::normalize(glm::vec3(1, 1, 0)));
    Gaussian gaussian4(glm::vec3(1, 1, 0), cov1, 1.0f);

    float analytical4 = computeGaussianIntegral(ray4, gaussian4);
    float monteCarlo4 = monteCarloIntegral(ray4, gaussian4, -10.0f, 15.0f, 1000000);

    std::cout << "Analytical result: " << analytical4 << std::endl;
    std::cout << "Monte Carlo result: " << monteCarlo4 << std::endl;
    std::cout << "Relative error: " << std::abs(analytical4 - monteCarlo4) / std::max(analytical4, monteCarlo4) * 100 << "%" << std::endl;

    // Additional test: Simple sanity check with a very peaked Gaussian
    std::cout << "\n=== Test Case 5: Sanity Check - Peaked Gaussian ===" << std::endl;
    Ray ray5(glm::vec3(0, 0, 0), glm::vec3(1, 0, 0));
    glm::mat3 cov5(0.01f, 0.0f, 0.0f,
        0.0f, 0.01f, 0.0f,
        0.0f, 0.0f, 0.01f);
    Gaussian gaussian5(glm::vec3(1, 0, 0), cov5, 1.0f);

    debugGaussian(gaussian5);

    float analytical5 = computeGaussianIntegral(ray5, gaussian5);
    float monteCarlo5 = monteCarloIntegral(ray5, gaussian5, -5.0f, 10.0f, 2000000);

    std::cout << "Analytical result: " << analytical5 << std::endl;
    std::cout << "Monte Carlo result: " << monteCarlo5 << std::endl;
    std::cout << "Relative error: " << std::abs(analytical5 - monteCarlo5) / std::max(analytical5, monteCarlo5) * 100 << "%" << std::endl;

    return 0;
}