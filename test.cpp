//#include <iostream>
//
//
//#include <fstream>
//
//#include "happly.h"
//#include "Math.h"
//
//Colour GaussianColor(Ray& ray, std::vector<Gaussian> gaussians) {}; 
//float computeGaussianIntegral(Ray& ray, Gaussian& gaussian) {};
//Colour evaluateSphericalHarmonics(const Vec3& viewDir, Gaussian& gaussian) {};
//
//
//void printAllGaussianDetails(std::vector<Gaussian>& gaussians) {
//
//	std::ofstream outFile("gaussians.txt");
//
//	std::cout << "Number of Gaussians: " << gaussians.size() << std::endl;
//	for (size_t i = 0; i < gaussians.size(); i++) {
//		outFile << "Gaussian " << i << ": Position: ("
//			<< gaussians[i].pos.x << ", "
//			<< gaussians[i].pos.y << ", "
//			<< gaussians[i].pos.z << "), "
//			<< "Opacity: " << gaussians[i].opacity
//			<< ", Scale: (" << gaussians[i].scale.x << ", " << gaussians[i].scale.y << ", " << gaussians[i].scale.z
//			<< "), Rotation: (" << gaussians[i].rotation.x << ", " << gaussians[i].rotation.y << ", " << gaussians[i].rotation.z << ", " << gaussians[i].rotation.w
//			<< "), ZeroSH: (" << gaussians[i].ZeroSH.x << ", " << gaussians[i].ZeroSH.y << ", " << gaussians[i].ZeroSH.z << ")"
//			<< ", HigherSH: [";
//		for (const float& sh : gaussians[i].higherSH) {
//			outFile << sh << ", ";
//		}
//		outFile << "]\n" << std::endl;
//	}
//	outFile.close();
//}
//
//
//void testScaleActivation() {
//	happly::PLYData plyIn("point_cloud.ply");
//	std::vector<float> elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("scale_0");
//	std::vector<float> elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("scale_1");
//	std::vector<float> elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("scale_2");
//
//	float maxScale = 0;
//	float maxExpScale = 0;
//	float maxSigScale = 0;
//
//
//	//for (size_t i = 0; i < elementA_prop1.size(); i++) {
//	for (size_t i = 0; i < 50; i++) {
//		Vec3 scale = Vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
//
//		Vec3 expScale = scale.exponent();
//		expScale = expScale * 10.f;
//		Vec3 sigScale = scale.Sigmoid();
//		sigScale = sigScale * 10.f;
//		std::cout << "Scale: " << scale.x << ", " << scale.y << ", " << scale.z << " ExScale: " << expScale.x << ", " << expScale.y << ", " << expScale.z << " Sigmoid Scale: " << sigScale.x << ", " << sigScale.y << ", " << sigScale.z << std::endl;
//		std::cout << std::endl;
//
//		maxScale = (std::max)(maxScale, scale._max());
//		maxExpScale = (std::max)(maxExpScale, expScale._max());
//		maxSigScale = (std::max)(maxSigScale, sigScale._max());
//	}
//	std::cout << "Max Scale: " << maxScale << std::endl;
//	std::cout << "Max Exponential Scale: " << maxExpScale << std::endl;
//	std::cout << "Max Sigmoid Scale: " << maxSigScale << std::endl;
//
//}
//
//void testGauss() {
//	//	Vec3 pos				 normal				ZeroSH			rotation		scale			opacity;
//	std::vector<Gaussian> testGaussians = {
//	{ { 0.0f, 0.0f,  2.0f }, { 0, 0, 1 }, { 1.0f, 0.0f, 0.0f }, {0, 0, 0, 1}, {0.2f, 0.2f, 0.2f}, 0.9f },
//	{ { 0.5f, 0.3f, 2.5f }, { 0, 1, 0 }, { 0.0f, 1.0f, 0.0f }, {0, 0.707f, 0, 0.707f}, {0.3f, 0.2f, 0.1f}, 0.8f },
//	{ { -0.4f, -0.3f, 2.2f }, { 1, 0, 0 }, { 0.0f, 0.0f, 1.0f }, {0.707f, 0, 0, 0.707f}, {0.1f, 0.3f, 0.2f}, 0.7f },
//	{ { 0.2f, -0.5f, 2.8f }, { 0, 1, 0 }, { 1.0f, 1.0f, 0.0f }, {0.5f, 0.5f, 0.5f, 0.5f}, {0.15f, 0.15f, 0.3f}, 0.85f },
//	{ { -0.6f, 0.4f, 3.0f }, { 0, 0, 1 }, { 0.0f, 1.0f, 1.0f }, {0, 1, 0, 0}, {0.25f, 0.1f, 0.25f}, 0.75f },
//	{ { 0.3f, 0.2f, 1.8f }, { 1, 0, 1 }, { 1.0f, 0.5f, 0.0f }, {0.707f, 0.707f, 0, 0}, {0.2f, 0.2f, 0.2f}, 1.0f },
//	{ { -0.3f, -0.2f, 2.1f }, { 1, 1, 0 }, { 0.5f, 0.5f, 1.0f }, {0.923f, 0, 0.383f, 0}, {0.3f, 0.1f, 0.2f}, 0.6f },
//	{ { 0.1f, 0.1f, 1.9f }, { 0, 0, 1 }, { 0.7f, 0.7f, 0.7f }, {1, 0, 0, 0}, {0.25f, 0.25f, 0.25f}, 0.95f },
//	{ { -0.1f, 0.3f, 2.4f }, { 0, 1, 0 }, { 0.3f, 0.6f, 0.9f }, {0, 0, 0, 1}, {0.2f, 0.3f, 0.2f}, 0.65f },
//	{ { 0.0f, -0.1f, 2.6f }, { 0, 1, 1 }, { 1.0f, 0.2f, 0.5f }, {0.5f, 0, 0.5f, 0.707f}, {0.3f, 0.3f, 0.15f}, 0.9f }
//	};
//
//	Ray r(Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f));
//
//	Colour color = GaussianColor(r, testGaussians);
//
//	std::cout << "Computed Color: ("
//		<< color.r << ", "
//		<< color.g << ", "
//		<< color.b << ")" << std::endl;
//}
//
//
//Colour testDensityandColorPrint(std::ofstream& outFile, Ray& ray, std::vector<Gaussian>& gaussians) {
//	outFile << "Testing Density and Color for Ray: (" << ray.o.x << ", " << ray.o.y << ", " << ray.o.z << ") Direction: ("
//		<< ray.dir.x << ", " << ray.dir.y << ", " << ray.dir.z << ")" << std::endl;
//	float tr = 1.0f;
//	Colour color(0.0f, 0.0f, 0.0f);
//
//	std::vector<Gaussian> sorted = gaussians;
//	std::sort(sorted.begin(), sorted.end(), [&](const Gaussian& a, const Gaussian& b) {
//		return (a.pos - ray.o).lengthSq() < (b.pos - ray.o).lengthSq();
//		});
//
//	for (Gaussian& gaussian : sorted) {
//		float density = computeGaussianIntegral(ray, gaussian);
//
//		float alpha = density * gaussian.opacity;
//
//		alpha = 1.0f - exp(-density);
//
//
//		Vec3 viewDir = (ray.o - gaussian.pos).normalize();
//		Colour SHColor = evaluateSphericalHarmonics(viewDir, gaussian);
//
//		bool intersect = gaussian.aabb.rayAABB(ray);
//
//
//		if (tr < 0.001f) {
//			break; // Stop if transmittance is very low
//		}
//
//		if (density < 1e-6f) {
//			continue; // Skip if density is negligible
//		}
//
//		color = color + (SHColor * alpha * tr);
//		//color = color.normalize();
//		color.correct();
//		tr *= (1.0f - alpha);
//
//		outFile << "Gaussian Index: " << gaussian.index << " with density: " << density
//			<< " at position: (" << gaussian.pos.x << ", " << gaussian.pos.y << ", " << gaussian.pos.z << ")"
//			<< " and color: (" << SHColor.r << ", " << SHColor.g << ", " << SHColor.b << ")"
//			<< " intersect: " << (intersect ? "true" : "false") << " alpha: " << alpha << " Transmittance: " << tr << "\n\n";
//	}
//	return color;
//}
//
//void testDensityandColor(std::vector<Gaussian>& gaussians, Camera& camera, RTCamera& viewCamera, BVHNode& bvh) {
//	std::ofstream outFile("IntersectedGaussians.txt");
//
//	outFile << "Camera Origin: (" << camera.origin.x << ", " << camera.origin.y << ", " << camera.origin.z << ")"
//		<< " View Direction: (" << camera.viewDirection.x << ", " << camera.viewDirection.y << ", " << camera.viewDirection.z << ")"
//		<< " From: (" << viewCamera.from.x << ", " << viewCamera.from.y << ", " << viewCamera.from.z << ")"
//		<< " To: (" << viewCamera.to.x << ", " << viewCamera.to.y << ", " << viewCamera.to.z << ")"
//		<< " Up: (" << viewCamera.up.x << ", " << viewCamera.up.y << ", " << viewCamera.up.z << ")\n\n"
//		<< std::endl;
//
//
//	Ray ray = camera.generateRay(0.5, 0.5);
//	bvh.traverse(ray, gaussians);
//	Colour color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
//	outFile << "Final Color at screen pos 0.5,0.5 :(" << color.r << ", " << color.g << ", " << color.b << ")with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;
//
//	ray = camera.generateRay(50.5, 50.5);
//	bvh.traverse(ray, gaussians);
//	color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
//	outFile << "Final Color at screen pos 50.5,50.5 :(" << color.r << ", " << color.g << ", " << color.b << ")with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;
//
//	ray = camera.generateRay(100.5, 100.5);
//	bvh.traverse(ray, gaussians);
//	color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
//	outFile << "Final Color at screen pos 100.5,100.5 :(" << color.r << ", " << color.g << ", " << color.b << ") with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;
//
//	ray = camera.generateRay(150.5, 150.5);
//	bvh.traverse(ray, gaussians);
//	color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
//	outFile << "Final Color at screen pos 150.5,150.5 :(" << color.r << ", " << color.g << ", " << color.b << ")with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;
//
//	ray = camera.generateRay(199.5, 199.5);
//	bvh.traverse(ray, gaussians);
//	color = testDensityandColorPrint(outFile, ray, bvh.getIntersectedGaussians());
//	outFile << "Final Color at screen pos 199.5,199.5 :(" << color.r << ", " << color.g << ", " << color.b << ")with " << bvh.getIntersectedGaussians().size() << " intersected gaussians\n\n" << std::endl;
//
//	outFile.close();
//
//
//}
//
//
//void empty(Camera& camera, RTCamera& viewcamera) {
//	std::ofstream outFile("CameraPosition.txt");
//
//	for (int i = 0; i < 15; i++) {
//		viewcamera.left();
//	}
//	for (int i = 0; i < 5; i++) {
//		viewcamera.back();
//	}
//
//	outFile << "Camera Origin: (" << camera.origin.x << ", " << camera.origin.y << ", " << camera.origin.z << ")"
//		<< " View Direction: (" << camera.viewDirection.x << ", " << camera.viewDirection.y << ", " << camera.viewDirection.z << ")"
//		<< " From: (" << viewcamera.from.x << ", " << viewcamera.from.y << ", " << viewcamera.from.z << ")"
//		<< " To: (" << viewcamera.to.x << ", " << viewcamera.to.y << ", " << viewcamera.to.z << ")"
//		<< " Up: (" << viewcamera.up.x << ", " << viewcamera.up.y << ", " << viewcamera.up.z << ")\n\n"
//		<< std::endl;
//
//}