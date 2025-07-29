#include <iostream>


#include <fstream>

#include "happly.h"
#include "Math.h"
#include "Imaging.h"


#define SH_C0 0.28209479177387814f
#define SH_C1 0.4886025119029199f

#define SH_C2_0 1.0925484305920792f
#define SH_C2_1 -1.0925484305920792f
#define SH_C2_2 0.31539156525252005f
#define SH_C2_3 -1.0925484305920792f
#define SH_C2_4 0.5462742152960396f

#define SH_C3_0 -0.5900435899266435f
#define SH_C3_1 2.890611442640554f
#define SH_C3_2 -0.4570457994644658f
#define SH_C3_3 0.3731763325901154f
#define SH_C3_4 -0.4570457994644658f
#define SH_C3_5 1.445305721320277f
#define SH_C3_6 -0.5900435899266435f

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

	glm::vec3 _r = r.ToGlm();
	glm::vec3 _dir = ray.dir.ToGlm();

	glm::vec3 __a = gaussian.covariance3D * _dir;
	float a2 = glm::dot(__a, _dir);

	Vec3 _a = gaussian.covariance.mulRowVec(ray.dir);
	float a = _a.dot(ray.dir);

	Vec3 _c = gaussian.covariance.mulRowVec(r);
	float c = _c.dot(r);

	glm::vec3 __c = gaussian.covariance3D * _r;
	float c2 = glm::dot(__c, _r);


	float b = _c.dot(ray.dir) ;
	float b2 = glm::dot(__c, _dir);

	float firstPart = std::expf((SQ(b) / (2 * a) - (c / 2)));
	float firstPart2 = std::expf((SQ(b2) / (2 * a2) - (c2 / 2)));

	float secondPart = std::sqrtf(M_PI / (2*a));
	float secondPart2 = std::sqrtf(M_PI / (2 * a2));

	//std::cout << "Gaussian Integral: " << firstPart << " * " << secondPart << " * " << gaussian.opacity ;

	//return firstPart * secondPart;
	return firstPart2 * secondPart2 * 0.5;
}

void evaluateSphericalHarmonics(const Vec3& viewDir, std::vector<Gaussian>& gaussians) {
	Vec3 dir = viewDir.normalize();
	float x = dir.x;
	float y = dir.y;
	float z = dir.z;

	float xx = x * x;
	float yy = y * y;
	float zz = z * z;
	float xy = x * y;
	float xz = x * z;
	float yz = y * z;

	for (Gaussian& gaussian : gaussians) {
		Vec3 color =
			Vec3(gaussian.higherSH[0], gaussian.higherSH[1], gaussian.higherSH[2]) * SH_C0;

		color = color

			- (Vec3(gaussian.higherSH[0], gaussian.higherSH[1], gaussian.higherSH[2]) * SH_C1 * y)
			+ (Vec3(gaussian.higherSH[3], gaussian.higherSH[4], gaussian.higherSH[5]) * SH_C1 * z)
			- (Vec3(gaussian.higherSH[6], gaussian.higherSH[7], gaussian.higherSH[8]) * SH_C1 * x);

		color = color


			+ (Vec3(gaussian.higherSH[9], gaussian.higherSH[10], gaussian.higherSH[11]) * SH_C2_0 * xy)
			+ (Vec3(gaussian.higherSH[12], gaussian.higherSH[13], gaussian.higherSH[14]) * SH_C2_1 * yz)
			+ (Vec3(gaussian.higherSH[15], gaussian.higherSH[16], gaussian.higherSH[17]) * SH_C2_2 * (2.0f * zz - xx - yy))
			+ (Vec3(gaussian.higherSH[18], gaussian.higherSH[19], gaussian.higherSH[20]) * SH_C2_3 * xz)
			+ (Vec3(gaussian.higherSH[21], gaussian.higherSH[22], gaussian.higherSH[23]) * SH_C2_4 * (xx - yy));

		color = color

			+ (Vec3(gaussian.higherSH[24], gaussian.higherSH[25], gaussian.higherSH[26]) * SH_C3_0 * y * (3.0f * xx - yy))
			+ (Vec3(gaussian.higherSH[27], gaussian.higherSH[28], gaussian.higherSH[29]) * SH_C3_1 * xy * z)
			+ (Vec3(gaussian.higherSH[30], gaussian.higherSH[31], gaussian.higherSH[32]) * SH_C3_2 * y * (4.0f * zz - xx - yy))
			+ (Vec3(gaussian.higherSH[33], gaussian.higherSH[34], gaussian.higherSH[35]) * SH_C3_3 * z * (2.0f * zz - 3.0f * xx - 3.0f * yy))
			+ (Vec3(gaussian.higherSH[36], gaussian.higherSH[37], gaussian.higherSH[38]) * SH_C3_4 * x * (4.0f * zz - xx - yy))
			+ (Vec3(gaussian.higherSH[39], gaussian.higherSH[40], gaussian.higherSH[41]) * SH_C3_5 * z * (xx - yy))
			+ (Vec3(gaussian.higherSH[42], gaussian.higherSH[43], gaussian.higherSH[44]) * SH_C3_6 * x * (xx - 3.0f * yy) );

		gaussian.color += color;
		//gaussian.color += Vec3(0.5f);


		//std::cout<<"Gaussian color : (" << gaussian.color.r << ", " << gaussian.color.g << ", " << gaussian.color.b << ")" << std::endl;
	}
}

Colour evaluateSphericalHarmonics(const Vec3& viewDir, Gaussian& gaussian) {
	Vec3 dir = viewDir.normalize();
	float x = dir.x;
	float y = dir.y;
	float z = dir.z;

	float xx = x * x;
	float yy = y * y;
	float zz = z * z;
	float xy = x * y;
	float xz = x * z;
	float yz = y * z;

	Vec3 color =
		gaussian.ZeroSH * SH_C0;

	color = color

		- (Vec3(gaussian.higherSH[0], gaussian.higherSH[1], gaussian.higherSH[2]) * SH_C1 * y)
		+ (Vec3(gaussian.higherSH[3], gaussian.higherSH[4], gaussian.higherSH[5]) * SH_C1 * z)
		- (Vec3(gaussian.higherSH[6], gaussian.higherSH[7], gaussian.higherSH[8]) * SH_C1 * x);

	color = color


		+ (Vec3(gaussian.higherSH[9], gaussian.higherSH[10], gaussian.higherSH[11]) * SH_C2_0 * xy)
		+ (Vec3(gaussian.higherSH[12], gaussian.higherSH[13], gaussian.higherSH[14]) * SH_C2_1 * yz)
		+ (Vec3(gaussian.higherSH[15], gaussian.higherSH[16], gaussian.higherSH[17]) * SH_C2_2 * (2.0f * zz - xx - yy))
		+ (Vec3(gaussian.higherSH[18], gaussian.higherSH[19], gaussian.higherSH[20]) * SH_C2_3 * xz)
		+ (Vec3(gaussian.higherSH[21], gaussian.higherSH[22], gaussian.higherSH[23]) * SH_C2_4 * (xx - yy));

	color = color

		+ (Vec3(gaussian.higherSH[24], gaussian.higherSH[25], gaussian.higherSH[26]) * SH_C3_0 * y * (3.0f * xx - yy))
		+ (Vec3(gaussian.higherSH[27], gaussian.higherSH[28], gaussian.higherSH[29]) * SH_C3_1 * xy * z)
		+ (Vec3(gaussian.higherSH[30], gaussian.higherSH[31], gaussian.higherSH[32]) * SH_C3_2 * y * (4.0f * zz - xx - yy))
		+ (Vec3(gaussian.higherSH[33], gaussian.higherSH[34], gaussian.higherSH[35]) * SH_C3_3 * z * (2.0f * zz - 3.0f * xx - 3.0f * yy))
		+ (Vec3(gaussian.higherSH[36], gaussian.higherSH[37], gaussian.higherSH[38]) * SH_C3_4 * x * (4.0f * zz - xx - yy))
		+ (Vec3(gaussian.higherSH[39], gaussian.higherSH[40], gaussian.higherSH[41]) * SH_C3_5 * z * (xx - yy))
		+ (Vec3(gaussian.higherSH[42], gaussian.higherSH[43], gaussian.higherSH[44]) * SH_C3_6 * x * (xx - 3.0f * yy));

	Colour c;
	c += color;
	//c += Vec3(0.5f);

	//std::cout<<"Gaussian color : (" << gaussian.color.r << ", " << gaussian.color.g << ", " << gaussian.color.b << ")" << std::endl;
	return c;
}

Colour viewIndependent(Gaussian& gaussian) {

	Vec3 color = gaussian.ZeroSH * SH_C0;
	Colour c;
	c += color;
	c += Vec3(0.5f);
	return c;
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

		density *= gaussian.opacity;

		float alpha = 1.0f - exp(-density);
		

		if (tr < 0.001f) {
			break; // Stop if transmittance is very low
		}

		if (density < 1e-6f){
			continue; // Skip if density is negligible
		}
		Vec3 viewDir = (ray.o - gaussian.pos).normalize();
		Colour SHColor = evaluateSphericalHarmonics(viewDir, gaussian);
		//Colour SHColor = viewIndependent(gaussian);

		color = color + (SHColor * alpha * tr);
		//color = color.normalize();
		color.correct();
		tr *= (1.0f - alpha); // Update transmittance
		//std::cout << " density: " << density <<" alpha: " << alpha<< " Transmittance: " << tr << std::endl;
		//std::cout << "Color: (" << color.r << ", " << color.g << ", " << color.b << ")" << " gaussian color: (" << gaussian.color.r << ", " << gaussian.color.g << ", " << gaussian.color.b << ")" << std::endl;

	}
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
		//gaussians[i].color += gaussians[i].ZeroSH;
		//gaussians[i].color = gaussians[i].color * 0.28f;
		//gaussians[i].color += Vec3(0.5f, 0.5f, 0.5f); // Adding a base color

		//std::cout << "ZeroSH: " << gaussians[i].ZeroSH.x << ", " << gaussians[i].ZeroSH.y << ", " << gaussians[i].ZeroSH.z << std::endl;
		//std::cout << "Color: " << gaussians[i].color.r << ", " << gaussians[i].color.g << ", " << gaussians[i].color.b << std::endl;
	}

	elementA_prop1 = plyIn.getElement("vertex").getProperty<float>("scale_0");
	elementA_prop2 = plyIn.getElement("vertex").getProperty<float>("scale_1");
	elementA_prop3 = plyIn.getElement("vertex").getProperty<float>("scale_2");

	for (size_t i = 0; i < elementA_prop1.size(); i++) {
		gaussians[i].scale = Vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i]);
		//std::cout << "Scale: " << gaussians[i].scale.x << ", " << gaussians[i].scale.y << ", " << gaussians[i].scale.z << std::endl;
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
		gaussians[i].rotation = Vec3(elementA_prop1[i], elementA_prop2[i], elementA_prop3[i], elementA_prop4[i]);
		gaussians[i].compute_gaussian_covariance();
	}


	for (int i = 0; i < 45; i++) {
		std::string name = "f_rest_" + std::to_string(i);

		elementA_prop1 = plyIn.getElement("vertex").getProperty<float>(name);

		for (size_t i = 0; i < elementA_prop1.size(); i++) {
			gaussians[i].higherSH.push_back(elementA_prop1[i]);
		}
	}


	//deleete

	/*for (size_t i = 0; i < elementA_prop1.size(); i++) {
		std::cout << "Higher SH: ";
		std::cout << gaussians[i].ZeroSH.x << ", " << gaussians[i].ZeroSH.y << ", " << gaussians[i].ZeroSH.z << std::endl;

		for (float a : gaussians[i].higherSH) {
			
			std::cout << a << ", ";
		}
		std::cout << std::endl;

	}*/


	//printing to file

}


void setCamera(Camera &camera, RTCamera& viewCamera) {
	//Vec3 from(4.63008f, 2.00661f, 2.06466f);
	Vec3 from(4.63008f, 2.00661f, 20.06466f);
	viewCamera.from = from;

	//float pitch = -21.5944f * (M_PI / 180.0f);  // X rotation
	float pitch = -5.5944f * (M_PI / 180.0f);  // X rotation
	//float yaw = 65.9668f * (M_PI / 180.0f);     // Y rotation
	float yaw = 11.9668f * (M_PI / 180.0f);     // Y rotation
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
	int width = 500;
	//int height = 768;
	int height = 500;
	canvas.create(width, height, "Charles GE");
	float fov = 45;


	//testScaleActivation();
	//return 0;

	Matrix P = Matrix::perspective(0.001f, 10000.0f, (float)width / (float)height, fov);

	Camera camera;
	camera.init(P, width, height);

	RTCamera viewcamera;
	
	setCamera(camera, viewcamera);

	for (int i = 0; i < 19; i++) {
		viewcamera.left();
	}
	for (int i = 0; i < 10; i++) {
		viewcamera.forward();
	}
	for (int i = 0; i < 1; i++) {
		viewcamera.flyDown();
	}
	
	for (int i = 0; i < 3; i++) {
		viewcamera.strafeLeft();
	}
	for (int i = 0; i < 3; i++) {
		viewcamera.left();
	}

	for (int i = 0; i < 4; i++) {
		viewcamera.forward();
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


	//std::cout << "Evaluating Spherical Harmonics:\n";
	//evaluateSphericalHarmonics(camera.viewDirection, gaussians);


	//testDensityandColor(gaussians , camera, viewcamera,bvh);
	

	bool running = true;
	int loopCount = 0;
	std::string filename = "19l_f10_1d_sl3_l3_4f_SH_half";
	std::string fullFilename;
	std::cout << "Rendering Gaussians...\n";

	std::ofstream outFile("CameraPositionsLog.txt");

	while (running)
	{
		canvas.checkInput();
		canvas.clear();

		timer.reset();


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
			//print percentage at intervals of 10%
			if (y % (height / 10) == 0) {
				std::cout << "Rendering progress: " << (y * 100 / height) << "%\r";
			}

		}
		float t = timer.dt();
		std::cout << "\nRendering time count "<<loopCount<<": " << t << std::endl;

		fullFilename = filename + "_" + std::to_string(loopCount++) + ".png";
		savePNG(fullFilename, &canvas);
		canvas.present();

		

		//print froma and to
		outFile << loopCount <<" - Camera From: (" << viewcamera.from.x << ", " << viewcamera.from.y << ", " << viewcamera.from.z << ")"
			<< " To: (" << viewcamera.to.x << ", " << viewcamera.to.y << ", " << viewcamera.to.z << ")"
			<< " Up: (" << viewcamera.up.x << ", " << viewcamera.up.y << ", " << viewcamera.up.z << ")"
			<< std::endl;


		if (loopCount > 30)
		{
			std::cout << "Stopping after 30 loops." << std::endl;
			break;
		}


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
		
		
		viewcamera.strafeRight();


	}
	outFile.close();

	return 0;
}