#include <iostream>
#include "happly.h"
#include "Math.h"
#include "GamesEngineeringBase.h"




float GaussianDensity(Ray& ray, const Gaussian& gaussian)
{
	Vec3 diff = gaussian.pos - ray.o;
	float t = diff.dot(ray.dir);
	Vec3 p = ray.at(t) - gaussian.pos;

	//Mahalanobis distance

	float dist = p.x * p.x / gaussian.covariance.a[0][0] +
		p.y * p.y / gaussian.covariance.a[1][1] +
		p.z * p.z / gaussian.covariance.a[2][2];

	return expf(-0.5f * dist) * gaussian.opacity;
}

float computeGaussianIntegral(Ray& ray , Gaussian& gaussian) {
	Vec3 r = gaussian.pos - ray.o;
	
	ray.dir = ray.dir.normalize();

	Vec3 _a = gaussian.covariance.mulRowVec(ray.dir);
	float a = _a.dot(ray.dir);

	Vec3 _c = gaussian.covariance.mulRowVec(r);
	float c = _c.dot(r);


	float b = _c.dot(ray.dir) ;

	float firstPart = std::expf((SQ(b) / (2 * a) - (c / 2)));

	float secondPart = std::sqrtf(M_PI / (2*a));

	std::cout << "Gaussian Integral: " << firstPart << " * " << secondPart << " * " << gaussian.opacity ;

	return firstPart * secondPart * gaussian.opacity;
}

Colour GaussianColor(Ray& ray, std::vector<Gaussian> gaussians)
{
	Colour color(0.0f, 0.0f, 0.0f);
	float tr = 1.0f; // Transmittance, assuming full opacity for simplicity

	std::vector<Gaussian> sorted = gaussians;
	std::sort(sorted.begin(), sorted.end(), [&](const Gaussian& a, const Gaussian& b) {
		return (a.pos - ray.o).lengthSq()< (b.pos - ray.o).lengthSq();
	});

	for (Gaussian& gaussian : sorted) {
		float density = computeGaussianIntegral(ray, gaussian);
		//float density = GaussianDensity(ray, gaussian);

		float alpha = 1.0f - exp(-density);
		

		if (tr < 0.001f) {
			break; // Stop if transmittance is very low
		}

		if (density < 1e-6f){
			continue; // Skip if density is negligible
		}

		color = color + ( gaussian.color * alpha * tr);
		//color = color.normalize();
		color.correct();
		tr *= (1.0f - alpha); // Update transmittance
		std::cout << " density: " << density <<" alpha: " << alpha<< " Transmittance: " << tr << std::endl;
		std::cout << "Color: (" << color.r << ", " << color.g << ", " << color.b << ")" << " gaussian color: (" << gaussian.color.r << ", " << gaussian.color.g << ", " << gaussian.color.b << ")" << std::endl;

	}
	return color;
}

void evaluateSphericalHarmonics(const Vec3& pos, const Vec3& normal, Colour& shColor) {
	// Placeholder for SH evaluation logic
	// This function should compute the spherical harmonics based on the position and normal
	// For now, we will just set a dummy color
	shColor = Colour(1.0f, 0.5f, 0.5f); // Example color
}

void rendersplats() {
	//read splats#
	//gen ray
	//get all inetsected gaussian
	//compute color
	//render color
}

void parsePLY(std::string filename, std::vector<Gaussian>& gaussians) {
	happly::PLYData plyIn("point_cloud.ply");
	std::vector<float> elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("x");
	std::vector<float> elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("y");
	std::vector<float> elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("z");

	int size = elementA_prop1.size();

	gaussians = std::vector<Gaussian>(size);

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		//std::cout << "x: " << elementA_prop1[i] << ", y: " << elementA_prop2[i] << ", z: " << elementA_prop3[i] << std::endl;
		gaussians[i].pos = Vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
		gaussians[i].index = i;
	}

	/*
	elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("nx");
	elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("ny");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("nz");

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i].normal = Vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
	}
	*/

	elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("f_dc_0");
	elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("f_dc_1");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("f_dc_2");


	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i].ZeroSH = Vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
		gaussians[i].color += gaussians[i].ZeroSH;
		gaussians[i].color = gaussians[i].color * 0.28f;
		gaussians[i].color += Vec3(0.5f, 0.5f, 0.5f); // Adding a base color

		std::cout << "ZeroSH: " << gaussians[i].ZeroSH.x << ", " << gaussians[i].ZeroSH.y << ", " << gaussians[i].ZeroSH.z << std::endl;
		std::cout << "Color: " << gaussians[i].color.r << ", " << gaussians[i].color.g << ", " << gaussians[i].color.b << std::endl;
	}

	elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("scale_0");
	elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("scale_1");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("scale_2");

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i].scale = Vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
		gaussians[i].compute_gaussian_aabb();
	}

	std::vector<float> elementA_prop4 = plyIn.getElement("vertex").getProperty<float>("opacity");

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i].opacity =	sigmoid(elementA_prop4[i]);
		//std::cout << "Opacity: " << gaussians[i].opacity << std::endl;
	}

	elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("rot_0");
	elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("rot_1");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("rot_2");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("rot_3");

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i].rotation = Vec4(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i], elementA_prop4[i]);
		gaussians[i].compute_gaussian_covariance();
	}
}

void testGauss() {
	//	Vec3 pos				 normal				ZeroSH			rotation		scale			opacity;
	std::vector<Gaussian> testGaussians = {
	{ { 0.0f, 0.0f,  2.0f }, { 0, 0, 1 }, { 1.0f, 0.0f, 0.0f }, {0, 0, 0, 1}, {0.2f, 0.2f, 0.2f}, 0.9f },
	{ { 0.5f, 0.3f, 2.5f }, { 0, 1, 0 }, { 0.0f, 1.0f, 0.0f }, {0, 0.707f, 0, 0.707f}, {0.3f, 0.2f, 0.1f}, 0.8f },
	{ { -0.4f, -0.3f, 2.2f }, { 1, 0, 0 }, { 0.0f, 0.0f, 1.0f }, {0.707f, 0, 0, 0.707f}, {0.1f, 0.3f, 0.2f}, 0.7f },
	{ { 0.2f, -0.5f, 2.8f }, { 0, 1, 0 }, { 1.0f, 1.0f, 0.0f }, {0.5f, 0.5f, 0.5f, 0.5f}, {0.15f, 0.15f, 0.3f}, 0.85f },
	{ { -0.6f, 0.4f, 3.0f }, { 0, 0, 1 }, { 0.0f, 1.0f, 1.0f }, {0, 1, 0, 0}, {0.25f, 0.1f, 0.25f}, 0.75f },
	{ { 0.3f, 0.2f, 1.8f }, { 1, 0, 1 }, { 1.0f, 0.5f, 0.0f }, {0.707f, 0.707f, 0, 0}, {0.2f, 0.2f, 0.2f}, 1.0f },
	{ { -0.3f, -0.2f, 2.1f }, { 1, 1, 0 }, { 0.5f, 0.5f, 1.0f }, {0.923f, 0, 0.383f, 0}, {0.3f, 0.1f, 0.2f}, 0.6f },
	{ { 0.1f, 0.1f, 1.9f }, { 0, 0, 1 }, { 0.7f, 0.7f, 0.7f }, {1, 0, 0, 0}, {0.25f, 0.25f, 0.25f}, 0.95f },
	{ { -0.1f, 0.3f, 2.4f }, { 0, 1, 0 }, { 0.3f, 0.6f, 0.9f }, {0, 0, 0, 1}, {0.2f, 0.3f, 0.2f}, 0.65f },
	{ { 0.0f, -0.1f, 2.6f }, { 0, 1, 1 }, { 1.0f, 0.2f, 0.5f }, {0.5f, 0, 0.5f, 0.707f}, {0.3f, 0.3f, 0.15f}, 0.9f }
	};

	Ray r(Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 1.0f));

	Colour color = GaussianColor(r, testGaussians);

	std::cout << "Computed Color: ("
		<< color.r << ", "
		<< color.g << ", "
		<< color.b << ")" << std::endl;
}

void setCamera(Camera &camera, RTCamera& viewCamera) {
	Vec3 from(4.63008f, 2.00661f, 2.06466f);
	viewCamera.from = from;

	float pitch = -21.5944f * (M_PI / 180.0f);  // X rotation
	float yaw = 65.9668f * (M_PI / 180.0f);     // Y rotation
	float roll = 0.0f * (M_PI / 180.0f);        // Z rotation

	Vec3 forward(
		cos(pitch) * cos(yaw),
		sin(pitch),
		cos(pitch) * sin(yaw)
	);

	Vec3 to = from + forward;
	viewCamera.to = to;

	Vec3 up(0.0f, 1.0f, 0.0f);
	viewCamera.up = up;

	viewCamera.setCamera(&camera);
}

int main() {
	GamesEngineeringBase::Window canvas;
	GamesEngineeringBase::Timer tim;

	//can be called inside camera
	//int width = 1024;
	int width = 200;
	//int height = 768;
	int height = 200;
	canvas.create(width, height, "Charles GE");
	float fov = 45;


	Matrix P = Matrix::perspective(0.001f, 10000.0f, (float)width / (float)height, fov);

	Camera camera;
	camera.init(P, width, height);

	RTCamera viewcamera;
	
	setCamera(camera, viewcamera);


	std::cout << "Parsing PLY file...\n";
	std::vector<Gaussian> gaussians{};
	parsePLY("point_cloud.ply", gaussians);
	//parsePLY("perf.ply", gaussians);
	std::cout << "Done PLY file...\n";


	std::cout << "Building BVH...\n";
	BVHNode bvh;
	bvh.build(gaussians);
	std::cout << "Done building BVH...\n";

	std::cout << "BVH Statistics:\n";
	bvh.checkTraverse();
	bvh.printStat();
	GamesEngineeringBase::Timer timer;
	

	bool running = true;
	while (running)
	{
		canvas.checkInput();
		canvas.clear();

		timer.reset();

		std::cout << "Rendering Gaussians...\n";

		for (unsigned int y = 0; y < height; y++)
		{
			for (unsigned int x = 0; x < width; x++)
			{
				float px = x + 0.5f;
				float py = y + 0.5f;
				Ray ray = camera.generateRay(px, py);

				bvh.traverse(ray, gaussians);

				Colour color = GaussianColor(ray, bvh.getIntersectedGaussians());

				canvas.draw(x, y, color.r * 255.0f, color.g * 255.0f, color.b * 255.0f);

				//std::cout << "\nRendering pixel (" << x << ", " << y << ") - Color: ("<< color.r << ", "<< color.g << ", "<< color.b << ")\n";
			}
		}
		float t = timer.dt();
		std::cout << "rendeiring time: " << t << std::endl;
		canvas.present();


		if (canvas.keyPressed(VK_ESCAPE))
		{
			break;
		}
		if (canvas.keyPressed('W'))
		{
			viewcamera.forward();
		}
		if (canvas.keyPressed('S'))
		{
			viewcamera.back();
		}
		if (canvas.keyPressed('A'))
		{
			viewcamera.left();
		}
		if (canvas.keyPressed('D'))
		{
			viewcamera.right();
		}
		if (canvas.keyPressed('Q'))
		{
			break;
		}

	}


	return 0;
}