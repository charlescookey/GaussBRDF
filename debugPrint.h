#pragma once

#include "happly.h"
#include "Math.h"
#include "Imaging.h"



void printAllGaussianDetails(std::vector<Gaussian>& gaussians) {

	std::ofstream outFile("gaussians.txt");

	std::cout << "Number of Gaussians: " << gaussians.size() << std::endl;
	for (size_t i = 0; i < gaussians.size(); i++) {
		outFile << "Gaussian " << i << ": Position: ("
			<< gaussians[i]._position.x << ", "
			<< gaussians[i]._position.y << ", "
			<< gaussians[i]._position.z << "), "
			<< "Opacity: " << gaussians[i].opacity
			<< ", Scale: (" << gaussians[i]._scale.x << ", " << gaussians[i]._scale.y << ", " << gaussians[i]._scale.z
			<< "), Rotation: (" << gaussians[i]._rotation.x << ", " << gaussians[i]._rotation.y << ", " << gaussians[i]._rotation.z << ", " << gaussians[i]._rotation.w
			<< "), ZeroSH: (" << gaussians[i].ZeroSH.x << ", " << gaussians[i].ZeroSH.y << ", " << gaussians[i].ZeroSH.z << ")"
			<< ", HigherSH: [";
		for (const float& sh : gaussians[i].higherSH) {
			outFile << sh << ", ";
		}
		outFile << "]\n" << std::endl;
	}
	outFile.close();
}

void testGauss() {
	//std::vector<Gaussian> testGaussians = {
	//{ { 0.0f, 0.0f,  2.0f }, { 0, 0, 1 }, { 1.0f, 0.0f, 0.0f }, {0, 0, 0, 1}, {0.2f, 0.2f, 0.2f}, 0.9f },
	//{ { 0.5f, 0.3f, 2.5f }, { 0, 1, 0 }, { 0.0f, 1.0f, 0.0f }, {0, 0.707f, 0, 0.707f}, {0.3f, 0.2f, 0.1f}, 0.8f },
	//{ { -0.4f, -0.3f, 2.2f }, { 1, 0, 0 }, { 0.0f, 0.0f, 1.0f }, {0.707f, 0, 0, 0.707f}, {0.1f, 0.3f, 0.2f}, 0.7f },
	//{ { 0.2f, -0.5f, 2.8f }, { 0, 1, 0 }, { 1.0f, 1.0f, 0.0f }, {0.5f, 0.5f, 0.5f, 0.5f}, {0.15f, 0.15f, 0.3f}, 0.85f },
	//{ { -0.6f, 0.4f, 3.0f }, { 0, 0, 1 }, { 0.0f, 1.0f, 1.0f }, {0, 1, 0, 0}, {0.25f, 0.1f, 0.25f}, 0.75f },
	//{ { 0.3f, 0.2f, 1.8f }, { 1, 0, 1 }, { 1.0f, 0.5f, 0.0f }, {0.707f, 0.707f, 0, 0}, {0.2f, 0.2f, 0.2f}, 1.0f },
	//{ { -0.3f, -0.2f, 2.1f }, { 1, 1, 0 }, { 0.5f, 0.5f, 1.0f }, {0.923f, 0, 0.383f, 0}, {0.3f, 0.1f, 0.2f}, 0.6f },
	//{ { 0.1f, 0.1f, 1.9f }, { 0, 0, 1 }, { 0.7f, 0.7f, 0.7f }, {1, 0, 0, 0}, {0.25f, 0.25f, 0.25f}, 0.95f },
	//{ { -0.1f, 0.3f, 2.4f }, { 0, 1, 0 }, { 0.3f, 0.6f, 0.9f }, {0, 0, 0, 1}, {0.2f, 0.3f, 0.2f}, 0.65f },
	//{ { 0.0f, -0.1f, 2.6f }, { 0, 1, 1 }, { 1.0f, 0.2f, 0.5f }, {0.5f, 0, 0.5f, 0.707f}, {0.3f, 0.3f, 0.15f}, 0.9f }
	//};

	//Ray r(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f));

	//Colour color = GaussianColor(r, testGaussians);

	//std::cout << "Computed Color: ("
	//	<< color.r << ", "
	//	<< color.g << ", "
	//	<< color.b << ")" << std::endl;
}

glm::vec3 testDensityandColorPrint(std::ofstream& outFile, Ray& ray, std::vector<Gaussian>& gaussians) {
	outFile << "Testing Density and Color for Ray: (" << ray.o.x << ", " << ray.o.y << ", " << ray.o.z << ") Direction: ("
		<< ray.dir.x << ", " << ray.dir.y << ", " << ray.dir.z << ")" << std::endl;
	float tr = 1.0f;
	glm::vec3 color(0.0f, 0.0f, 0.0f);

	std::vector<Gaussian> sorted = gaussians;
	std::sort(sorted.begin(), sorted.end(), [&](const Gaussian& a, const Gaussian& b) {
		return glm::length(a._position - ray.o) < glm::length(b._position - ray.o);
		});

	for (Gaussian& gaussian : sorted) {
		float density = computeGaussianIntegral(ray, gaussian);

		float alpha = density * gaussian.opacity;

		alpha = 1.0f - exp(-density);


		glm::vec3 viewDir = glm::normalize(ray.o - gaussian._position);
		glm::vec3 SHColor = evaluateSphericalHarmonics(viewDir, gaussian);

		bool intersect = gaussian.aabb.rayAABB(ray);


		if (tr < 0.001f) {
			break; // Stop if transmittance is very low
		}

		if (density < 1e-6f) {
			continue; // Skip if density is negligible
		}

		color = color + (SHColor * alpha * tr);
		//color = color.normalize();
		correct(color);
		tr *= (1.0f - alpha);

		outFile << "Gaussian Index: " << gaussian.index << " with density: " << density
			<< " at position: (" << gaussian._position.x << ", " << gaussian._position.y << ", " << gaussian._position.z << ")"
			<< " and color: (" << SHColor.r << ", " << SHColor.g << ", " << SHColor.b << ")"
			<< " intersect: " << (intersect ? "true" : "false") << " alpha: " << alpha << " Transmittance: " << tr << "\n\n";
	}
	return color;
}

void testDensityandColor(std::vector<Gaussian>& gaussians, Camera& camera, RTCamera& viewCamera, BVHNode& bvh) {
	std::ofstream outFile("IntersectedGaussians.txt");

	outFile << "Camera Origin: (" << camera.origin.x << ", " << camera.origin.y << ", " << camera.origin.z << ")"
		<< " View Direction: (" << camera.viewDirection.x << ", " << camera.viewDirection.y << ", " << camera.viewDirection.z << ")"
		<< " From: (" << viewCamera.from.x << ", " << viewCamera.from.y << ", " << viewCamera.from.z << ")"
		<< " To: (" << viewCamera.to.x << ", " << viewCamera.to.y << ", " << viewCamera.to.z << ")"
		<< " Up: (" << viewCamera.up.x << ", " << viewCamera.up.y << ", " << viewCamera.up.z << ")\n\n"
		<< std::endl;


	Ray ray = camera.generateRay(0.5, 0.5);
	bvh.traverse(ray, gaussians);
	glm::vec3 color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
	outFile << "Final Color at screen pos 0.5,0.5 :(" << color.r << ", " << color.g << ", " << color.b << ")with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;

	ray = camera.generateRay(50.5, 50.5);
	bvh.traverse(ray, gaussians);
	color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
	outFile << "Final Color at screen pos 50.5,50.5 :(" << color.r << ", " << color.g << ", " << color.b << ")with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;

	ray = camera.generateRay(100.5, 100.5);
	bvh.traverse(ray, gaussians);
	color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
	outFile << "Final Color at screen pos 100.5,100.5 :(" << color.r << ", " << color.g << ", " << color.b << ") with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;

	ray = camera.generateRay(150.5, 150.5);
	bvh.traverse(ray, gaussians);
	color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
	outFile << "Final Color at screen pos 150.5,150.5 :(" << color.r << ", " << color.g << ", " << color.b << ")with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;

	ray = camera.generateRay(199.5, 199.5);
	bvh.traverse(ray, gaussians);
	color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
	outFile << "Final Color at screen pos 199.5,199.5 :(" << color.r << ", " << color.g << ", " << color.b << ")with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;

	outFile.close();


}


void empty(Camera& camera, RTCamera& viewcamera) {
	std::ofstream outFile("CameraPosition.txt");

	for (int i = 0; i < 15; i++) {
		viewcamera.left();
	}
	for (int i = 0; i < 5; i++) {
		viewcamera.back();
	}

	outFile << "Camera Origin: (" << camera.origin.x << ", " << camera.origin.y << ", " << camera.origin.z << ")"
		<< " View Direction: (" << camera.viewDirection.x << ", " << camera.viewDirection.y << ", " << camera.viewDirection.z << ")"
		<< " From: (" << viewcamera.from.x << ", " << viewcamera.from.y << ", " << viewcamera.from.z << ")"
		<< " To: (" << viewcamera.to.x << ", " << viewcamera.to.y << ", " << viewcamera.to.z << ")"
		<< " Up: (" << viewcamera.up.x << ", " << viewcamera.up.y << ", " << viewcamera.up.z << ")\n\n"
		<< std::endl;

}