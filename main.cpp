#include <iostream>


#include <fstream>

#include "happly.h"
#include "Math.h"
#include "Imaging.h"


float GaussianDensity(Ray& ray, const Gaussian& gaussian)
{
	glm::vec3 diff = gaussian._position - ray.o;
	float t = glm::dot(diff , ray.dir);
	glm::vec3 p = ray.at(t) - gaussian._position;

	//Mahalanobis distance

	float dist = p.x * p.x / gaussian.covariance3D[0][0] +
		p.y * p.y / gaussian.covariance3D[1][1] +
		p.z * p.z / gaussian.covariance3D[2][2];

	return expf(-0.5f * dist);
}
float GaussianDensity2(Ray& ray, Gaussian& gaussian)
{
	glm::vec3 diff = gaussian._position - ray.o;
	float t = glm::dot(diff , ray.dir); // projection onto ray
	glm::vec3 p = ray.at(t) - gaussian._position;

	// Diagonal Mahalanobis distance
	float dist = (p.x * p.x) / gaussian.covariance3D[0][0] +
		(p.y * p.y) / gaussian.covariance3D[1][1] +
		(p.z * p.z) / gaussian.covariance3D[2][2];

	// Determinant for normalization (diagonal case)
	float det = gaussian.covariance3D[0][0] *
		gaussian.covariance3D[1][1] *
		gaussian.covariance3D[2][2];

	float norm = 1.0f / (powf(2.0f * M_PI, 1.5f) * sqrtf(det));

	float density = norm * expf(-0.5f * dist);

	return density;
}



float GaussianDensity3(Ray& ray, Gaussian& g) {
	glm::vec3 diff = g._position - ray.o;
	float t = glm::dot(diff, ray.dir);  // projection of diff onto ray.dir
	glm::vec3 x = ray.at(t) - g._position;    // vector from Gaussian center to ray point

	glm::vec3 scaled = glm::vec3(
		glm::dot(x, g.R[0]) / g._scaleMuld.x,
		glm::dot(x, g.R[1]) / g._scaleMuld.y,
		glm::dot(x, g.R[2]) / g._scaleMuld.z
	);


	// Use double to avoid exp underflow
	double dist2 = glm::dot(glm::dvec3(scaled), glm::dvec3(scaled));
	double det = static_cast<double>(g._scaleMuld.x * g._scaleMuld.y * g._scaleMuld.z);
	double norm = 1.0 / (pow(2.0 * M_PI, 1.5) * det);

	float density =  static_cast<float>(norm * exp(-0.5 * dist2));


	// Debug output
	//std::cout << "Gaussian Density3: " << norm << " * exp(-0.5 * " << dist2
	//	<< ") = " << density
	//	<< " scaled: (" << scaled.x << ", " << scaled.y << ", " << scaled.z << ")"
	//	<< " det: " << det
	//	<< " Gaussian scale: (" << g._scaleMuld.x << ", " << g._scaleMuld.y << ", " << g._scaleMuld.z << ")"
	//	<< std::endl;

	return density;
}

float evaluateGaussian(Gaussian& splat, Ray& ray) {
	glm::vec3 to_splat = splat._position - ray.o;
	float closest_approach_t = glm::dot(to_splat, ray.dir);
	glm::vec3 _point = ray.at(closest_approach_t);
	glm::vec3 point = glm::vec3(_point.x, _point.y, _point.z);

	glm::vec3 diff = point - splat._position;

	// Compute inverse of covariance matrix for numerical stability
	//glm::mat3 cov_inv = glm::inverse(splat.RawCovariance3D + 0.001f * glm::mat3(1.0f));

	float mahal_dist = glm::dot(diff, splat.RawCovariance3D_inv * diff);

	// Gaussian density (unnormalized)
	return exp(-0.5f * mahal_dist);
}

float _evaluateGaussian(Gaussian& splat, glm::vec3& point) {
	glm::vec3 diff = point - splat._position;

	// Compute inverse of covariance matrix for numerical stability
	glm::mat3 cov_inv = glm::inverse(splat.covariance3D + 0.001f * glm::mat3(1.0f));

	float mahal_dist = glm::dot(diff, cov_inv * diff);

	// Gaussian density (unnormalized)
	return exp(-0.5f * mahal_dist);
}

float computeGaussianIntegral(Ray& ray, Gaussian& gaussian) {
	glm::vec3 r = gaussian._position - ray.o;

	ray.dir = glm::normalize(ray.dir);
	glm::vec3 _a = gaussian.inverse_covariance3D * ray.dir;
	float a = glm::dot(_a, ray.dir);

	glm::vec3 _c = gaussian.inverse_covariance3D * r;
	float c = glm::dot(_c , r);


	float b = glm::dot(_c , ray.dir);

	float firstPart = std::expf((SQ(b) / (2 * a) - (c / 2)));

	float secondPart = std::sqrtf(M_PI / (2 * a));

	//std::cout << "Gaussian Integral: " << firstPart << " * " << secondPart << " * " << gaussian.opacity ;

	return firstPart * secondPart;
}

glm::vec3 evaluateSphericalHarmonics(const glm::vec3& viewDir, Gaussian& gaussian) {

	glm::vec3 dir = glm::normalize(viewDir);
	float x = dir.x;
	float y = dir.y;
	float z = dir.z;

	float xx = x * x;
	float yy = y * y;
	float zz = z * z;
	float xy = x * y;
	float xz = x * z;
	float yz = y * z;

	glm::vec3 color =
		gaussian.ZeroSH * SH_C0;

	color = color

		- (glm::vec3(gaussian.higherSH[0], gaussian.higherSH[1], gaussian.higherSH[2]) * SH_C1 * y)
		+ (glm::vec3(gaussian.higherSH[3], gaussian.higherSH[4], gaussian.higherSH[5]) * SH_C1 * z)
		- (glm::vec3(gaussian.higherSH[6], gaussian.higherSH[7], gaussian.higherSH[8]) * SH_C1 * x);

	color = color


		+ (glm::vec3(gaussian.higherSH[9], gaussian.higherSH[10], gaussian.higherSH[11]) * SH_C2_0 * xy)
		+ (glm::vec3(gaussian.higherSH[12], gaussian.higherSH[13], gaussian.higherSH[14]) * SH_C2_1 * yz)
		+ (glm::vec3(gaussian.higherSH[15], gaussian.higherSH[16], gaussian.higherSH[17]) * SH_C2_2 * (2.0f * zz - xx - yy))
		+ (glm::vec3(gaussian.higherSH[18], gaussian.higherSH[19], gaussian.higherSH[20]) * SH_C2_3 * xz)
		+ (glm::vec3(gaussian.higherSH[21], gaussian.higherSH[22], gaussian.higherSH[23]) * SH_C2_4 * (xx - yy));

	color = color

		+ (glm::vec3(gaussian.higherSH[24], gaussian.higherSH[25], gaussian.higherSH[26]) * SH_C3_0 * y * (3.0f * xx - yy))
		+ (glm::vec3(gaussian.higherSH[27], gaussian.higherSH[28], gaussian.higherSH[29]) * SH_C3_1 * xy * z)
		+ (glm::vec3(gaussian.higherSH[30], gaussian.higherSH[31], gaussian.higherSH[32]) * SH_C3_2 * y * (4.0f * zz - xx - yy))
		+ (glm::vec3(gaussian.higherSH[33], gaussian.higherSH[34], gaussian.higherSH[35]) * SH_C3_3 * z * (2.0f * zz - 3.0f * xx - 3.0f * yy))
		+ (glm::vec3(gaussian.higherSH[36], gaussian.higherSH[37], gaussian.higherSH[38]) * SH_C3_4 * x * (4.0f * zz - xx - yy))
		+ (glm::vec3(gaussian.higherSH[39], gaussian.higherSH[40], gaussian.higherSH[41]) * SH_C3_5 * z * (xx - yy))
		+ (glm::vec3(gaussian.higherSH[42], gaussian.higherSH[43], gaussian.higherSH[44]) * SH_C3_6 * x * (xx - 3.0f * yy));


	//std::cout<<"Gaussian color : (" << gaussian.color.r << ", " << gaussian.color.g << ", " << gaussian.color.b << ")" << std::endl;
	return color;
}



glm::vec3 GaussianColor(Ray& ray, std::vector<Gaussian> gaussians)
{
	glm::vec3 color(0.0f, 0.0f, 0.0f);
	float tr = 1.0f; // Transmittance, assuming full opacity for simplicity

	std::vector<Gaussian> sorted = gaussians;
	std::sort(sorted.begin(), sorted.end(), [&](const Gaussian& a, const Gaussian& b) {
		return glm::pow(glm::length(a._position - ray.o),2) < glm::pow(glm::length(b._position - ray.o),2);
		});

	for (Gaussian& gaussian : sorted) {
		float density = computeGaussianIntegral(ray, gaussian);
		//float density = GaussianDensity3(ray, gaussian) * 1e4;
		//float density = GaussianDensity(ray, gaussian);
		//float density = evaluateGaussian(gaussian, ray);


		//std::cout << "Density: " << density << " Density2: " << density2 << " Density3 : " << density3 << " Density4 : " << density4 << std::endl;

		density *= gaussian.opacity;


		if (density < 1e-6f) {
			continue; // Skip if density is negligible
		}

		float alpha = 1.0f - expf(-density);
		//alpha *=gaussian.opacity;

		if (tr < 0.001f) {
			break; // Stop if transmittance is very low
		}

		glm::vec3 viewDir = glm::normalize(ray.o - gaussian._position);
		//Colour SHColor = evaluateSphericalHarmonics(viewDir, gaussian);

		color = color + (gaussian._color * alpha * tr);
		//color = color.normalize();
		correct(color);
		tr *= (1.0f - alpha); // Update transmittance
	}
	toneMap(color);
	return color;
}

glm::vec3 GaussianColor3(Ray& ray, std::vector<Gaussian>& gaussians)
{
	glm::vec3 color(0.0f, 0.0f, 0.0f);
	float tr = 1.0f; // transmittance

	// Optional: sort by distance to ray origin
	std::vector<Gaussian> sorted = gaussians;
	std::sort(sorted.begin(), sorted.end(), [&](const Gaussian& a, const Gaussian& b) {
		return glm::pow(glm::length(a._position - ray.o), 2) < glm::pow(glm::length(b._position - ray.o), 2);
		});

	for (Gaussian& g : sorted)
	{
		float density = evaluateGaussian(g , ray);

		float alpha = density * g.opacity;


		if (density < 1e-6f || tr < 1e-3f)
			continue;
		// Use view-independent color or SH-evaluated if desired

		color += (g._color * alpha * tr);
		tr *= (1.0f - alpha);
	}
	toneMap(color);
	//std::cout << "Final Color: (" << color.r << ", " << color.g << ", " << color.b << ")" << std::endl;
	return color;
}


void parsePLY(std::string filename, std::vector<Gaussian>& gaussians) {
	happly::PLYData plyIn(filename.c_str());
	std::vector<float> elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("x");
	std::vector<float> elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("y");
	std::vector<float> elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("z");

	int size = elementA_prop1.size();

	gaussians = std::vector<Gaussian>(size);

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i]._position = glm::vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
		gaussians[i].index = i;
	}

	elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("f_dc_0");
	elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("f_dc_1");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("f_dc_2");


	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i].ZeroSH = glm::vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
		gaussians[i].viewIndependent();

		//std::cout << "ZeroSH: " << gaussians[i].ZeroSH.x << ", " << gaussians[i].ZeroSH.y << ", " << gaussians[i].ZeroSH.z << std::endl;
		//std::cout << "Color: " << gaussians[i].color.r << ", " << gaussians[i].color.g << ", " << gaussians[i].color.b << std::endl;
	}

	elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("scale_0");
	elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("scale_1");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("scale_2");


	float average_scale = 0.0f;

	float average_raw = 0.0f;
	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		glm::vec3 raw = glm::vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
		average_raw += raw.length();
	}
	average_raw /= elementA_prop1.size();

	float raw_multiplier = 1.0f / average_raw;

	// Step 2: Apply normalization and exponentiation
	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		glm::vec3 raw = glm::vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
		gaussians[i]._scaleRaw = raw;                  // Save raw if needed
		gaussians[i]._scale = glm::exp(raw);
		raw *= raw_multiplier;                // Normalize raw latent scale
		gaussians[i]._scaleMuld = glm::exp(raw);      // 
	}

	float pos_scale = 3.0f;  // Optional: increase to bring Gaussians into camera view
	//for (auto& g : gaussians) {
		//g.SclaedPos = g.pos * pos_scale;
		//g.scaled = g.scaled * pos_scale;  // Maintain consistent spatial unit scaling

		//std::cout << "Scale: " << g.scale.x << ", " << g.scale.y << ", " << g.scale.z << std::endl;
		//std::cout << "after multiplier: " << g.scaled.x << ", " << g.scaled.y << ", " << g.scaled.z << std::endl;
		//std::cout << "Pos: " << g.pos.x << ", " << g.pos.y << ", " << g.pos.z << std::endl;
		//std::cout << "SclaedPos: " << g.SclaedPos.x << ", " << g.SclaedPos.y << ", " << g.SclaedPos.z << std::endl;
	//}

	for (auto& g : gaussians) {
		g._scaleMuld = g._scaleMuld * pos_scale;
		//g.SclaedPos = g.pos * scale_multiplier;

	}

	std::vector<float> elementA_prop4 = plyIn.getElement("vertex").getProperty<float>("opacity");

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i].compute_gaussian_aabb();
		gaussians[i].opacity = sigmoid(elementA_prop4[i]);
		//std::cout << "Opacity: " << gaussians[i].opacity << std::endl;
	}

	elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("rot_0");
	elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("rot_1");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("rot_2");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("rot_3");

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i]._rotation = glm::quat(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i], elementA_prop4[i]);
		gaussians[i].compute_gaussian_covariance_glm();
	}


	for (int i = 0; i < 45; i++) {
		std::string name = "f_rest_" + std::to_string(i);

		elementA_prop1 = plyIn.getElement("vertex").getProperty<float>(name);

		for (size_t i = 0; i < elementA_prop1.size(); i++) {
			gaussians[i].higherSH.push_back(elementA_prop1[i]);
		}
	}
}



void setCamera(Camera& camera, RTCamera& viewCamera) {
	Vec3 from(4.63008f, -2.00661f, 20.06466f);
	viewCamera.from = from;

	//float pitch = -21.5944f * (M_PI / 180.0f);  // X rotation
	float pitch = -21.5944f * (M_PI / 180.0f);  // X rotation
	//float yaw = 65.9668f * (M_PI / 180.0f);     // Y rotation
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
	int width = 300;
	//int height = 768;
	int height = 300;
	canvas.create(width, height, "Charles GE");
	float fov = 45;


	//testScaleActivation();
	//return 0;

	Matrix P = Matrix::perspective(0.001f, 10000.0f, (float)width / (float)height, fov);

	Camera camera;
	camera.init(P, width, height);

	RTCamera viewcamera;

	setCamera(camera, viewcamera);

	for (int i = 0; i < 3; i++) {
		//viewcamera.right();
	}
	for (int i = 0; i < 2; i++) {
		//viewcamera.flyDown();
	}
	for (int i = 0; i < 50; i++) {
		//viewcamera.forward();
		//viewcamera.back();
	}

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
	int loopCount = 0;
	std::string filename = "s3_p1_vd_for30_fullres_inv_diff_";
	std::string fullFilename;
	std::cout << "Rendering Gaussians...\n";

	while (running)
	{
		canvas.checkInput();
		canvas.clear();

		if (canvas.keyPressed(VK_ESCAPE))
		{
			break;
		}
		if (canvas.keyPressed('W'))
		{
			viewcamera.forward();
			std::cout << "w pressed, moving forward\n";
		}
		if (canvas.keyPressed('S'))
		{
			viewcamera.back();
			std::cout << "s pressed, moving back\n";
		}
		if (canvas.keyPressed('A'))
		{
			viewcamera.left();
			std::cout << "a pressed, moving left\n";
		}
		if (canvas.keyPressed('D'))
		{
			viewcamera.right();
			std::cout << "d pressed, moving right\n";
		}
		if (canvas.keyPressed('Q'))
		{
			std::cout << "q pressed, Quit\n";
			break;
		}

		timer.reset();


		for (unsigned int y = 0; y < height; y++)
		{
			for (unsigned int x = 0; x < width; x++)
			{
				float px = x + 0.5f;
				float py = y + 0.5f;
				Ray ray = camera.generateRay(px, py);

				bvh.traverse(ray, gaussians);

				glm::vec3 color = GaussianColor(ray, bvh.getIntersectedGaussians());

				canvas.draw(x, y, color.r * 255.0f, color.g * 255.0f, color.b * 255.0f);

				//std::cout << "\nRendering pixel (" << x << ", " << y << ") - Color: ("<< color.r << ", "<< color.g << ", "<< color.b << ")\n";
			}
			//print percentage at intervals of 10%
			if (y % (height / 10) == 0) {
				std::cout << "Rendering progress: " << (y * 100 / height) << "%\r";
			}

		}
		float t = timer.dt();
		std::cout << "\nRendering time count " << loopCount << ": " << t << std::endl;

		fullFilename = filename + "_" + std::to_string(loopCount++) + ".png";
		savePNG(fullFilename, &canvas);
		canvas.present();




		if (loopCount % 2 == 0) {
			//viewcamera.forward();
		}
		else {
			//viewcamera.forward();
		}
		viewcamera.forward();


	}
	std::cout << "Camera Origin: (" << camera.origin.x << ", " << camera.origin.y << ", " << camera.origin.z << ")"
		<< " View Direction: (" << camera.viewDirection.x << ", " << camera.viewDirection.y << ", " << camera.viewDirection.z << ")"
		<< " From: (" << viewcamera.from.x << ", " << viewcamera.from.y << ", " << viewcamera.from.z << ")"
		<< " To: (" << viewcamera.to.x << ", " << viewcamera.to.y << ", " << viewcamera.to.z << ")"
		<< " Up: (" << viewcamera.up.x << ", " << viewcamera.up.y << ", " << viewcamera.up.z << ")\n\n"
		<< std::endl;

	return 0;
}