/***********************************************************
	 Starter code for Assignment 3

	 This code was originally written by Jack Wang for
			CSC418, SPRING 2005

		Implementations of functions in raytracer.h,
		and the main function which specifies the
		scene to be rendered.

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <stdlib.h>

#include <string>
#include <limits>

#include <time.h>

//Set some default image names
std::string TEXTURE_IMG2 = "Earth-Clouds.bmp";
std::string TEXTURE_IMG = "negy.bmp";
std::string TEXTURE_IMG3 = "nebula.bmp";


bool TEXTURE_MAP_FLAG = true;
bool REFRACTION_FLAG = true;
bool REFRACTION_FLAG_DEBUG = false;

//Constants to easily change for material reflectance and refractance
double LARGE_SPH_REFLECT = 0.1;
double  LARGE_SPH_REFRACT = 0.6;
double LARGE_SPH_REFRAC_INDX = 1.3;


//For easy change of light position in the scene for the defaul scene
bool LIGHT_DEFAULT = true;
bool LIGHT_TEST_REFRAC = false;
Point3D LIGHT_POS_TEST = Point3D( 0, 5, 0);

//Multiple scenes
void refraction_scene_1();
void wmonkey_scene_2();

//Different materials
Material RED( Colour(0, 0, 0), Colour(0.9, 0.05, 0.05),
             Colour(0.4, 0.2, 0.2),
             30.2, 0.0, 1.5, 0.9);
Material GREEN_TRANSP( Colour(0, 0, 0), Colour(0.2, 0.9, 0.2),
             Colour(0.5, 0.6, 0.5),
             30.2, 0.0, 1.3, 0.9);


//Following defaults are for depth of field
bool DEPTH_OF_FIELD_FLAG = true;
int NUM_RAYS = 5;       //number of random rays to shoot at the scene to get an approximation of depth of field


Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();

	// defaults
	setAAMode(NONE);
	setShadingMode(SCENE_MODE_SIGNATURE);
	setShadows(NONE);
	setReflDepth(0);
	setColorSpaceMode(NONE);

	_num_ss_samples = 4;

	_refractionMode = false;
}

Raytracer::~Raytracer() {
	delete _root;
}

void Raytracer::setShadingMode(int mode) {
	switch ( mode ) {
	case SCENE_MODE_SIGNATURE:
		printf("setShadingMode: SIGNATURE\n");
		break;

	case SCENE_MODE_DIFFUSE:
		printf("setShadingMode: DIFFUSE\n");
		break;

	case SCENE_MODE_PHONG:
		printf("setShadingMode: PHONG\n");
		break;
	}

	_shadingMode = mode;
}
int Raytracer::getShadingMode() {
	return _shadingMode;
}

void Raytracer::setAAMode(int mode) {
	switch ( mode ) {
	case AA_SUPER_SAMPLING:
		printf("setAAMode: AA_SUPER_SAMPLING\n");
		break;

	case NONE:
		printf("setAAMode: NONE\n");
		break;
	}

	_aaMode = mode;
}
int Raytracer::getAAMode() {
	return _aaMode;
}

void Raytracer::setShadows(int mode) {
	switch ( mode ) {
	case SHADOW_CAST:
		printf("setShadows: SHADOW_CAST\n");
		break;

	case NONE:
		printf("setShadows: NONE\n");
		break;
	}

	_shadows = mode;
}
int Raytracer::getShadows() {
	return _shadows;
}

void Raytracer::setReflDepth(int depth) {
	printf("setReflDepth: %d\n", depth);

	_MAX_DEPTH = depth;
}
int Raytracer::getReflDepth() {
	return _MAX_DEPTH;
}

void Raytracer::setEnvMap(EnvMap map) {
	_envmap = map;
}
EnvMap Raytracer::getEnvMap() {
	return _envmap;
}

// Set texture map for all the objects (for now, it's only for the sphere)
void Raytracer::setTextureMap(TextureMap map){
    _txtMapSphere = map;
}

TextureMap Raytracer::getTextureMapSph(){
    return _txtMapSphere;
}

// Refraction mode
void Raytracer::setRefractionMode(bool mode){
    printf("setting Refraction mode\n");
    _refractionMode = mode;
}
bool Raytracer::getRefractionMode(){
    return _refractionMode;
}

void Raytracer::setColorSpaceMode(int mode) {
	switch ( mode ) {
	case COLOR_ENC_SRGB_GAMMA_CORRECT:
		printf("setColorSpaceMode: COLOR_ENC_SRGB_GAMMA_CORRECT\n");
		break;

	case NONE:
		printf("setColorSpaceMode: NONE\n");
		break;
	}

	_colorSpace_Mode = mode;
}
int Raytracer::getColorSpaceMode() {
	return _colorSpace_Mode;
}
//_colorSpace_Mode

void Raytracer::setEnvMapMode(int mode) {
	switch ( mode ) {
	case ENV_MAP_CUBE_SKYBOX:
		printf("setEnvMapMode: ENV_MAP_CUBE_SKYBOX\n");
		break;

	case NONE:
		printf("setEnvMapMode: NONE\n");
		break;
	}

	_envMap_Mode = mode;
}
int Raytracer::getEnvMapMode() {
	return _envMap_Mode;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent,
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;

	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}

	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;

	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
			node->trans = node->trans*rotation;
			angle = -angle;
		}
		else {
			node->invtrans = rotation*node->invtrans;
		}
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;

	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation;
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans;
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;

	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale;
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans;
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view,
		Vector3D up ) {
	Matrix4x4 mat;
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat;
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	Raytracer::traverseScene(node, ray, false);
}

/**
 * Set onlyFirst=true for shadow checks
 */
void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray, bool onlyFirst ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel;
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
            //if texture mapping is on, find the corresponding pixel intesity from the Texture Map for the ray intersection point
            if (node->useTextureMapping && getShadingMode() == SCENE_MODE_PHONG_TEXTURE) {
                //find the colour of the intersection point from the texture map of the object
                Colour colourIntObj;
                colourIntObj = node->obj->computeColourAtImagePoint(ray.intersection.point, _worldToModel);
                //set abient, diffuse, and specular colour to the colour of the textureMap 
                node->mat->diffuse = colourIntObj;
                node->mat->specular = colourIntObj;
            }
            ray.intersection.mat = node->mat;            
		}
	}

	if ( !onlyFirst || ray.intersection.none ) {
		// Traverse the children.
		childPtr = node->child;
		while (childPtr != NULL) {
			traverseScene(childPtr, ray);
			childPtr = childPtr->next;
		}
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

/**
 * Calculate if we're within a shadow, and set the ray accordingly
 */
void Raytracer::setIfShadowed(LightListNode* curLight, Ray3D& ray)
{
//	if ( getShadows() == SHADOW_CAST ) {
//		/*
//		 * Basic idea is we create a ray from the intersection point to the
//		 * light, and see if we hit anything on the way.
//		 */
//		Ray3D lightray;
//		lightray.origin = ray.intersection.point;
//		lightray.dir = curLight->light->get_position() - ray.intersection.point;
//		lightray.dir.normalize();
//
//		traverseScene(_root, lightray, true);
//
//#if _DEBUG
//		printf("------- setIfShadowed[%d, %d]\n", ray.x, ray.y);
//
//		std::cout << "setIfShadowed: " << "Light Pos: " << curLight->light->get_position() << std::endl;
//		std::cout << "setIfShadowed: " << "Ray int point: " << ray.intersection.point << std::endl;
//
//		std::cout << "setIfShadowed: " << "LR Origin: " << lightray.origin << std::endl;
//		std::cout << "setIfShadowed: " << "LR DIR: " << lightray.dir << std::endl;
//		std::cout << "setIfShadowed: " << "LR int point: " << lightray.intersection.point << std::endl;
//#endif
//		// if we intersected with another object, the light is blocked
//		if (!lightray.intersection.none)
//		{
//#if _DEBUG
//			printf("SHADOW!\n");
//#endif
//			ray.shadowed = true;
//		}
//#if _DEBUG
//		printf("setIfShadowed: t_value: %f\n", lightray.intersection.t_value);
//#endif
//	}
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	while ( curLight != NULL ) {
		// do rest of shading
		switch ( getShadingMode() ) {
		case SCENE_MODE_SIGNATURE:
			curLight->light->shade_signature(ray);
			break;

		case SCENE_MODE_DIFFUSE:
			setIfShadowed(curLight, ray);
			curLight->light->shade_diffuse(ray);
			break;

		case SCENE_MODE_PHONG:
			setIfShadowed(curLight, ray);
			curLight->light->shade(ray);
			break;
		}

		curLight = curLight->next;
	}
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray, int reflectionDepth, int refractionDepth ) {
	Colour col(0.0, 0.0, 0.0);
	traverseScene(_root, ray);

	// Don't bother shading if the ray didn't hit
	// anything.
    if (REFRACTION_FLAG_DEBUG) {
        printf("\n\n            checking intersection       \n");
        printf("          incident ray: (%f, %f, %f) \n", ray.dir[0], ray.dir[1], ray.dir[2]);
    }

	Intersection inter = ray.intersection;
	if ( !inter.none ) {
		computeShading(ray);
        if (REFRACTION_FLAG_DEBUG) {
            printf("\n\n    ----------------------  INTERSECTION WITH OBJECT ------------------------------   \n\n");
        }

        //Spawn a ray for reflection
		if ( reflectionDepth < getReflDepth() ) {
			Vector3D n = inter.normal;
			Vector3D de = -ray.dir;
			Vector3D m = 2 * de.dot(n) * n - de; // reflect across normal

			// http://www.cdf.toronto.edu/~moore/csc418/Notes/Phong.pdf
			// reflection
			if ( inter.mat->refl_amt > 0.0f ) {
                if (REFRACTION_FLAG_DEBUG) {
                    printf("\n\n     ........  DOING REFLECTION ........   \n");                    
                }

				Ray3D refl_ray;
				refl_ray.origin = inter.point;// + (DBL_0_UPPER * m);
				refl_ray.dir = m;

				Colour rCol(0.0, 0.0, 0.0);
				rCol = shadeRay(refl_ray, reflectionDepth + 1, refractionDepth);

				ray.col = ray.col + (inter.mat->refl_amt * rCol);
			}

	        // Spawn a seperate ray for refraction
			if ( inter.mat->transparency > 0.0f ) {
				double th_1 = std::abs(de.dot(n));

				double c1, c2;
				if ( th_1 >= 0 ) {
					c1 = inter.mat->IOR;
					c2 = ray._cur_speed;

					ray._cur_speed = c1;
				} else {
					c2 = inter.mat->IOR;
					c1 = ray._cur_speed;

					ray._cur_speed = c2;
				}

				double IOR = c2 / c1;
	//				double IOR = 1; // will make it transparent, for testing

				double th_2 = sqrt(1.0 - (IOR * IOR) * (1.0 - (th_1 * th_1)));

				Ray3D refr_ray;
				refr_ray.dir = -de + IOR * (th_1 * n) + (-th_2 * n);
				refr_ray.origin = (inter.point) + (DBL_EPSILON * refr_ray.dir);

				Colour rCol(0.0, 0.0, 0.0);
				rCol = shadeRay(refr_ray, reflectionDepth + 1, refractionDepth + 1);

				ray.col = ray.col + (inter.mat->transparency * rCol);
			}

		}

		col = ray.col;
	} else {
		// env map
		if ( getEnvMapMode() == ENV_MAP_CUBE_SKYBOX ) {
			col = getEnvMap().getColorForRay(ray.dir);
		}
	}

	col.clamp();
	return col;
}

/**
 * Construct a ray from the eye through the image plane and get the shaded
 * resulting color.
 */
#if _DEBUG
Colour Raytracer::getColorFromRay(
		const Matrix4x4& viewToWorld,
		const Point3D& imagePlane,
		const Point3D& eye,
		int i, int j)
#else
Colour Raytracer::getColorFromRay(
		const Matrix4x4& viewToWorld,
		const Point3D& imagePlane,
		const Point3D& eye)
#endif
{
	// create ray from eye to place in image plane
	Ray3D ray;
	ray.origin = viewToWorld * imagePlane;
	ray.dir = ray.origin - eye;

	ray.dir.normalize();

#if _DEBUG
	// keep track of where the ray shot out from... to help with debugging later on
	ray.y = i;
	ray.x = j;
#endif

    Colour col;
    //If Depth of field is ON, also shoot few random rays around the campera pin hole pixel to approximate depth of field
    //http://ray-tracer-concept.blogspot.ca/2011/12/depth-of-field.html
    if (DEPTH_OF_FIELD_FLAG) {
        //generate NUM_RAYS random rays around the position
//        printf(" >>>>>>> ADDING DEPTH OF FIELD <<<<<<<<< \n");
        for (int di = 0; di < NUM_RAYS; di++) {
            //generate a random number within 5 pixel radius
            float du = rand()/float(RAND_MAX+1);//generating random number
            float dv = rand()/float(RAND_MAX+1);
            
            // creating new camera position(or ray start using jittering)
            float x = imagePlane[0] + du;
            float y = imagePlane[1] + dv;
            float z = imagePlane[2];
            //set the new pixel at the imagePlane
            Point3D blurr_img = Point3D(x, y, z);
            ray.origin = viewToWorld * blurr_img;
            ray.dir = ray.origin - eye;
            ray.dir.normalize();
            // create ray from eye to place in image plane
            col = col + shadeRay(ray, 0, 0);
        }
        //normalize for NUM_RAYS pixel colours
        double col_r = col[0]/ NUM_RAYS;
        double col_g = col[1]/ NUM_RAYS;
        double col_b = col[2]/ NUM_RAYS;
        col = Colour(col_r, col_g, col_b);
    }else{
        // shade pixel based on ray data
        col = shadeRay(ray, 0, 0);
    }

	return col;
}

/**
 * Gamma correct and return value
 */
int encodeTo_sRGB(double c) {
//	return int(powf(c, 0.4545454545454) * 255 + 0.5);
	int o = int(powf(c, 1/2.2) * 255);

	if ( o > 255 ) o = 255;
	return o;
}

void Raytracer::render_no_AA(int width, double factor, int height,
		const Matrix4x4& viewToWorld, const Point3D& eye)
{
	printf("\t -> No AA Mode\n");

	Point3D imagePlane;
	imagePlane[2] = -1;

	double w2 = -double(width)/2 + 0.5;
	double h2 = -double(height)/2 + 0.5;

#if _DEBUG
	// render out a single ray for debugging ----------------------------------

	int i = _scrHeight - 250;
	int j = 400;

	// Sets up ray origin and direction in view space,
	// image plane is at z = -1.
	imagePlane[0] = (w2 + j)/factor;
	imagePlane[1] = (h2 + i)/factor;

	// create ray from eye to place in image plane
	Colour col = getColorFromRay(viewToWorld, imagePlane, eye, i, j);

	_rbuffer[i*width+j] = int(col[0]*255);
	_gbuffer[i*width+j] = int(col[1]*255);
	_bbuffer[i*width+j] = int(col[2]*255);
#else
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {

			// Sets up ray origin and direction in view space,
			// image plane is at z = -1.
			imagePlane[0] = (w2 + j)/factor;
			imagePlane[1] = (h2 + i)/factor;

			// create ray from eye to place in image plane
			Colour col = getColorFromRay(viewToWorld, imagePlane, eye);
    

			switch (getColorSpaceMode()) {
				case COLOR_ENC_SRGB_GAMMA_CORRECT:
					_rbuffer[i*width+j] += int(encodeTo_sRGB(col[0]));
					_gbuffer[i*width+j] += int(encodeTo_sRGB(col[1]));
					_bbuffer[i*width+j] += int(encodeTo_sRGB(col[2]));
					break;

				case NONE:
					_rbuffer[i*width+j] += int(col[0]*255);
					_gbuffer[i*width+j] += int(col[1]*255);
					_bbuffer[i*width+j] += int(col[2]*255);
					break;
			}
		}

		if ( i % 50 == 0 ) printf("%d NO AA: Row Complete.\n", i);
	}
#endif
}

void Raytracer::render_SS(int width, double factor, int height,
		const Matrix4x4& viewToWorld, const Point3D& eye)
{
	printf("\t -> Super Sampling Mode\n");

	// based on super sampling from:
	// http://paulbourke.net/miscellaneous/aliasing/
	float ray_offset = 1.0f / _num_ss_samples * 2;
	float color_factor = 1.0f / _num_ss_samples;

	printf("\t\t Number SS Samples: %d\n", _num_ss_samples);
	printf("\t\t Ray Offset: %f\n", ray_offset);
	printf("\t\t Color Factor: %f\n", color_factor);

	Point3D imagePlane;
	imagePlane[2] = -1;

	double w2 = -double(width)/2 + 0.5;
	double h2 = -double(height)/2 + 0.5;

	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {

			// do extra super sampling
			for (float part_I = i; part_I < i + 1.0f; part_I += ray_offset) {
				for (float part_J = j; part_J < j + 1.0f; part_J += ray_offset) {
					imagePlane[0] = (w2 + part_J)/factor;
					imagePlane[1] = (h2 + part_I)/factor;

					// create ray from eye to place in image plane
#if _DEBUG
					Colour col = getColorFromRay(
							viewToWorld, imagePlane, eye, i, j);
#else
					Colour col = getColorFromRay(
							viewToWorld, imagePlane, eye);
#endif

					switch (getColorSpaceMode()) {
						case COLOR_ENC_SRGB_GAMMA_CORRECT:
							_rbuffer[i*width+j] += int(encodeTo_sRGB(col[0])*color_factor);
							_gbuffer[i*width+j] += int(encodeTo_sRGB(col[1])*color_factor);
							_bbuffer[i*width+j] += int(encodeTo_sRGB(col[2])*color_factor);
							break;

						case NONE:
							_rbuffer[i*width+j] += int(col[0]*255*color_factor);
							_gbuffer[i*width+j] += int(col[1]*255*color_factor);
							_bbuffer[i*width+j] += int(col[2]*255*color_factor);
							break;
					}
				}
			}
		} // for j

		if ( i % 50 == 0 ) printf("%d SS: Row Complete.\n", i);
	} // for i
}

void Raytracer::render( int width, int height, Point3D eye, Vector3D view,
		Vector3D up, double fov, char* fileName ) {
	printf("\n --------------------- \n RENDERING: %s\n", fileName);

	printf("\t view[%f, %f, %f]\n", view[0], view[1], view[2]);
	printf("\t eye[%f, %f, %f]\n", eye[0], eye[1], eye[2]);
	printf("\t up[%f, %f, %f]\n", up[0], up[1], up[2]);

	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);

	clock_t start = clock();

	// color each pixel
	switch ( getAAMode() ) {
		case Raytracer::AA_SUPER_SAMPLING:
			render_SS(width, factor, height, viewToWorld, eye);
			break;

		case Raytracer::NONE:
			render_no_AA(width, factor, height, viewToWorld, eye);
			break;
	}

	clock_t end = clock();

	double ms = ((end - start) * 1000)/CLOCKS_PER_SEC;

	printf("RENDER TIME :: %f ms\n", ms);

	flushPixelBuffer(fileName);
}

/**
 * Wooden Monkey Scene 1
 */
void wmonkey_scene_1()
{
	printf("WOODEN MONKEY SCENE : 1 ----------------------------------\n\n");
	Raytracer rt;
	int width = 16 * 20 * 2;
	int height = 12 * 20 * 2;

	// Camera parameters.
	Point3D eye1(0, 0, 1), eye2(4, 2, 1);
	Vector3D view1(0, 0, -1), view2(-4, -2, -6);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
			Colour(0.628281, 0.555802, 0.366065),
			51.2, 0.8 );
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63),
			Colour(0.316228, 0.316228, 0.316228),
			12.8);

	// Defines a point light source.
	double l0c = 0.5;
	PointLight * light0 = new PointLight(
			Point3D(-2, 2, 5),
			Colour(l0c, l0c, l0c),
			0.2);
	rt.addLightSource(light0);

	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = rt.addObject( new UnitSphere(), &gold );
	SceneDagNode* sphere2 = rt.addObject( new UnitSphere(), &gold );
	SceneDagNode* plane = rt.addObject( new UnitSquare(), &jade );

	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	rt.translate(sphere, Vector3D(0, 0, -5));
	rt.rotate(sphere, 'x', -45);
	rt.rotate(sphere, 'z', 45);
	rt.scale(sphere, Point3D(0, 0, 0), factor1);

	rt.translate(plane, Vector3D(0, 0, -7));
	rt.rotate(plane, 'z', 45);
	rt.scale(plane, Point3D(0, 0, 0), factor2);

	double f[3] = { 0.5, 0.5, 0.5 };
	rt.translate(sphere2, Vector3D(3, 0, -5));
	rt.scale(sphere2, Point3D(0, 0, 0), f);

	rt.setAAMode(Raytracer::AA_SUPER_SAMPLING);
	rt.setShadingMode(Raytracer::SCENE_MODE_PHONG);
	rt.setShadows(Raytracer::SHADOW_CAST);
	rt.setEnvMapMode(Raytracer::ENV_MAP_CUBE_SKYBOX);
	rt.setColorSpaceMode(Raytracer::COLOR_ENC_SRGB_GAMMA_CORRECT);
	rt.setReflDepth(4);

	if ( rt.getEnvMapMode() != Raytracer::NONE ) {
		// load images
		EnvMap env;
		if ( _DEBUG ) {
			env = EnvMap(
				"EnvMaps/DebugMaps/posx.bmp",
				"EnvMaps/DebugMaps/posy.bmp",
				"EnvMaps/DebugMaps/posz.bmp",
				"EnvMaps/DebugMaps/negx.bmp",
				"EnvMaps/DebugMaps/negy.bmp",
				"EnvMaps/DebugMaps/negz.bmp"
			);
		} else {
			env = EnvMap(
				"EnvMaps/SaintLazarusChurch/posx.bmp",
				"EnvMaps/SaintLazarusChurch/posy.bmp",
				"EnvMaps/SaintLazarusChurch/posz.bmp",
				"EnvMaps/SaintLazarusChurch/negx.bmp",
				"EnvMaps/SaintLazarusChurch/negy.bmp",
				"EnvMaps/SaintLazarusChurch/negz.bmp"
			);
		}

		rt.setEnvMap(env);
	}

	printf("WOODEN MONKEY SCENE : 1 :: Rendering...\n");
	rt.render(width, height, eye2, view2, up, fov, "wmonkey_1.bmp");
	printf("WOODEN MONKEY SCENE : 1 :: Done!\n");
}

int main(int argc, char* argv[])
{
	// Build your scene and setup your camera here, by calling
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the
	// assignment.
	Raytracer raytracer;
	int width = 16 * 20 * 2;
	int height = 12 * 20 * 2;

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye1(0, 0, 1), eye2(4, 2, 1);
	Vector3D view1(0, 0, -1), view2(-4, -2, -6);
//	Point3D eye1(0, 0, 1), eye2(4, 2, -6);
//	Vector3D view1(0, 0, -1), view2(-4, -2, 1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
			Colour(0.628281, 0.555802, 0.366065),
			51.2, LARGE_SPH_REFLECT, LARGE_SPH_REFRAC_INDX, LARGE_SPH_REFRACT);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63),
			Colour(0.316228, 0.316228, 0.316228),
			12.8);

    Material red( Colour(0, 0, 0), Colour(0.9, 0.05, 0.05),
                  Colour(0.4, 0.2, 0.2),
                  12.8);
    
	// Defines a point light source.
    Point3D light_pos;
    if (LIGHT_DEFAULT){
        light_pos = Point3D(0, 0, 5);    
    }else{
        light_pos = LIGHT_POS_TEST;
    }

	PointLight * light0 = new PointLight(
			light_pos,
			Colour(0.9, 0.9, 0.9),
			0.1);
	raytracer.addLightSource(light0);

	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
    
    //set the texture map for the objects of interest in the scene if texture map flag is ON
    if (TEXTURE_MAP_FLAG) {
        // load texture image
        TextureMap txtmp;
        txtmp = TextureMap(TEXTURE_IMG);
        raytracer.setTextureMap(txtmp);
        
        //for now, we are only using texture map for sphere
        sphere->useTextureMapping = true;
        sphere->obj->setTextureMap(txtmp);
    }


	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	raytracer.translate(sphere, Vector3D(0, 0, -5));
	raytracer.rotate(sphere, 'x', -45);
	raytracer.rotate(sphere, 'z', 45);
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane, Vector3D(0, 0, -7));
	raytracer.rotate(plane, 'z', 45);
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	double f[3] = { 0.5, 0.5, 0.5 };
	raytracer.translate(sphere2, Vector3D(0, 0, -8));
	raytracer.scale(sphere2, Point3D(0, 0, 0), f);
    

    bool DO_SIGNATURE 				= false;
    bool DO_SIGNATURE_SS 			= false;
    bool DO_DIFFUSE 				= false;
    bool DO_PHONG 					= false;
    bool DO_PHONG_SS 				= false;
    bool DO_FULL_FEATURED 			= false;
	bool DO_WOODEN_MONKEY_SCENES 	= true;
	bool DO_REFRACTION_SCENE 		= false;

	bool RENDER_FIRST_VIEW = true;
	bool RENDER_SECOND_VIEW = true;


	raytracer.setReflDepth(0);
	raytracer.setEnvMapMode(Raytracer::NONE);

	// render signature
	if ( DO_SIGNATURE ) {
		raytracer.setAAMode(Raytracer::NONE);
		raytracer.setShadingMode(Raytracer::SCENE_MODE_SIGNATURE);

		if ( RENDER_FIRST_VIEW )
			raytracer.render(width, height, eye1, view1, up, fov, "sig1.bmp");
		if ( RENDER_SECOND_VIEW )
			raytracer.render(width, height, eye2, view2, up, fov, "sig2.bmp");
	}

	// render signature with SS AA
	if ( DO_SIGNATURE_SS ) {
		raytracer.setAAMode(Raytracer::AA_SUPER_SAMPLING);
		raytracer.setShadingMode(Raytracer::SCENE_MODE_SIGNATURE);

		if ( RENDER_FIRST_VIEW )
			raytracer.render(width, height, eye1, view1, up, fov, "sigSS1.bmp");
		if ( RENDER_SECOND_VIEW )
			raytracer.render(width, height, eye2, view2, up, fov, "sigSS2.bmp");
	}

	// render diffuse
	if ( DO_DIFFUSE ) {
		raytracer.setAAMode(Raytracer::NONE);
		raytracer.setShadingMode(Raytracer::SCENE_MODE_DIFFUSE);

		if ( RENDER_FIRST_VIEW )
			raytracer.render(width, height, eye1, view1, up, fov, "diffuse1.bmp");
		if ( RENDER_SECOND_VIEW )
			raytracer.render(width, height, eye2, view2, up, fov, "diffuse2.bmp");
	}

	// render phong
	if ( DO_PHONG ) {
		raytracer.setAAMode(Raytracer::NONE);
		raytracer.setShadingMode(Raytracer::SCENE_MODE_PHONG);

		if ( RENDER_FIRST_VIEW )
			raytracer.render(width, height, eye1, view1, up, fov, "phong1.bmp");
		if ( RENDER_SECOND_VIEW )
			raytracer.render(width, height, eye2, view2, up, fov, "phong2.bmp");
	}

	// phong with super sampling AA
	if ( DO_PHONG_SS ) {
		raytracer.setAAMode(Raytracer::AA_SUPER_SAMPLING);
		raytracer.setShadingMode(Raytracer::SCENE_MODE_PHONG);

		if ( RENDER_FIRST_VIEW )
			raytracer.render(width, height, eye1, view1, up, fov, "phongSS1.bmp");
		if ( RENDER_SECOND_VIEW )
			raytracer.render(width, height, eye2, view2, up, fov, "phongSS2.bmp");
	}
    
    // refraction if it's turned on
    if (REFRACTION_FLAG) {
        raytracer.setRefractionMode(REFRACTION_FLAG);
    }

	// all features enabled or turned to max
	if ( DO_FULL_FEATURED ) {
		raytracer.setAAMode(Raytracer::NONE);
		raytracer.setAAMode(Raytracer::AA_SUPER_SAMPLING);
		raytracer.setShadingMode(Raytracer::SCENE_MODE_PHONG);
		raytracer.setShadows(Raytracer::SHADOW_CAST);
//		raytracer.setShadows(Raytracer::NONE);
		raytracer.setEnvMapMode(Raytracer::ENV_MAP_CUBE_SKYBOX);
//		raytracer.setEnvMapMode(Raytracer::NONE);
		raytracer.setReflDepth(4);

		if ( raytracer.getEnvMapMode() != Raytracer::NONE ) {
			// load images

			EnvMap env;
			if ( _DEBUG ) {
				env = EnvMap(
					"EnvMaps/DebugMaps/posx.bmp",
					"EnvMaps/DebugMaps/posy.bmp",
					"EnvMaps/DebugMaps/posz.bmp",
					"EnvMaps/DebugMaps/negx.bmp",
					"EnvMaps/DebugMaps/negy.bmp",
					"EnvMaps/DebugMaps/negz.bmp"
				);
			} else {
				env = EnvMap(
					"EnvMaps/SaintLazarusChurch/posx.bmp",
					"EnvMaps/SaintLazarusChurch/posy.bmp",
					"EnvMaps/SaintLazarusChurch/posz.bmp",
					"EnvMaps/SaintLazarusChurch/negx.bmp",
					"EnvMaps/SaintLazarusChurch/negy.bmp",
					"EnvMaps/SaintLazarusChurch/negz.bmp"
				);
			}

			raytracer.setEnvMap(env);
		}
        
		// adjust lighting?
		if ( raytracer.getReflDepth() > 0 ) {
			double l0i = 0.5;
			light0->setAmbient(Colour(l0i, l0i, l0i));
		}

		if ( RENDER_FIRST_VIEW )
			raytracer.render(width, height, eye1, view1, up, fov, "all1.bmp");
		if ( RENDER_SECOND_VIEW )
			raytracer.render(width, height, eye2, view2, up, fov, "all2.bmp");
	}

	// different scenes just for the wooden monkey thing
	if ( DO_WOODEN_MONKEY_SCENES ) {
//		wmonkey_scene_1();
        wmonkey_scene_2();
		// TODO add more scenes here as required...
	}
    
    //render the 2nd refraction scene
    if ( REFRACTION_FLAG && DO_REFRACTION_SCENE ) {
        refraction_scene_1();
    }

	printf("Press enter to terminate...\n");
	std::string s;
	std::getline(std::cin, s);

	return 0;
}


/**
 * Wooden Monkey Scene 1
 */
void refraction_scene_1()
{
	printf("REFRACTION SCENE : 1 ----------------------------------\n\n");
	Raytracer rt;
	int width = 16 * 20 * 2;
	int height = 12 * 20 * 2;
    
	// Camera parameters.
	Point3D eye1(0, 0, 1), eye2(4, 2, 1);
	Vector3D view1(0, 0, -1), view2(-4, -2, -6);
	Vector3D up(0, 1, 0);
	double fov = 60;
    
	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
                  Colour(0.628281, 0.555802, 0.366065),
                  51.2, LARGE_SPH_REFLECT, LARGE_SPH_REFRAC_INDX, LARGE_SPH_REFRACT);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63),
                  Colour(0.316228, 0.316228, 0.316228),
                  12.8);
    // Defines a material for shading.
	Material gold_nonRefract( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
                  Colour(0.628281, 0.555802, 0.366065),
                  51.2, 0.8 );
    
	// Defines a point light source.
	double l0c = 0.5;
	PointLight * light0 = new PointLight(
                                         Point3D(-2, 2, 5),
                                         Colour(l0c, l0c, l0c),
                                         0.2);
	rt.addLightSource(light0);
    
	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = rt.addObject( new UnitSphere(), &gold );
	SceneDagNode* sphere2 = rt.addObject( new UnitSphere(), &gold_nonRefract );
	SceneDagNode* plane = rt.addObject( new UnitSquare(), &jade );
	SceneDagNode* sphere3 = rt.addObject( new UnitSphere(), &RED);
	SceneDagNode* sphere4 = rt.addObject( new UnitSphere(), &GREEN_TRANSP);
    SceneDagNode* plane2 = rt.addObject( new UnitSquare(), &jade );
//    SceneDagNode* plane3 = rt.addObject( new UnitSquare(), &jade );
//    SceneDagNode* plane4 = rt.addObject( new UnitSquare(), &jade );
    
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	rt.translate(sphere, Vector3D(0, 0, -5));
	rt.rotate(sphere, 'x', -45);
	rt.rotate(sphere, 'z', 45);
	rt.scale(sphere, Point3D(0, 0, 0), factor1);
    
	rt.translate(plane, Vector3D(0, 0, -7));
	rt.rotate(plane, 'z', 45);
	rt.scale(plane, Point3D(0, 0, 0), factor2);
    
	double f[3] = { 0.5, 0.5, 0.5 };
	rt.translate(sphere2, Vector3D(3, 0, -5));
	rt.scale(sphere2, Point3D(0, 0, 0), f);

	rt.translate(sphere3, Vector3D(0, 2, -5));
	rt.scale(sphere3, Point3D(0, 0, 0), f);
    
	double f2[3] = { 0.6, 0.6, 0.6 };
    rt.translate(sphere4, Vector3D(-2, 1, -3));
	rt.scale(sphere4, Point3D(0, 0, 0), f2);
    
    double fp2[3] = { 3.0, 3.0, 3.0 };
    rt.translate(plane2,Vector3D(-4,1,-5));
    rt.rotate(plane2, 'z', 45);
    rt.rotate(plane2, 'y', 45);
	rt.scale(plane2, Point3D(0, 0, 0), fp2);

//    rt.translate(plane3,Vector3D(-2,0,-5));
//    rt.rotate(plane2, 'z', 45);
//    rt.rotate(plane3, 'x', 90);
//	rt.scale(plane3, Point3D(0, 0, 0), fp2);
//
//    rt.translate(plane4,Vector3D(-2,1,-5));
//    rt.rotate(plane2, 'z', 45);
//    rt.rotate(plane4, 'y', 90);
//	rt.scale(plane4, Point3D(0, 0, 0), fp2);
    
	rt.setAAMode(Raytracer::AA_SUPER_SAMPLING);
	rt.setShadingMode(Raytracer::SCENE_MODE_PHONG);
	rt.setShadows(Raytracer::SHADOW_CAST);
	rt.setEnvMapMode(Raytracer::ENV_MAP_CUBE_SKYBOX);
	rt.setColorSpaceMode(Raytracer::COLOR_ENC_SRGB_GAMMA_CORRECT);
	rt.setReflDepth(4);
    
    //set the texture map for the objects of interest in the scene if texture map flag is ON
    if (TEXTURE_MAP_FLAG) {
        // load texture image
        TextureMap txtmp;
        txtmp = TextureMap(TEXTURE_IMG);
        TextureMap txtmp2 = TextureMap(TEXTURE_IMG2);
        TextureMap txtmp3 = TextureMap(TEXTURE_IMG3);
        
        //for now, we are only using texture map for sphere
        sphere->useTextureMapping = true;
        sphere->obj->setTextureMap(txtmp);
        
        sphere2->useTextureMapping = false;
        
        sphere4->useTextureMapping = true;
        sphere4->setTextMapOfObject(txtmp2);
        
        plane2->useTextureMapping = true;
        plane2->setTextMapOfObject(txtmp3);
        
//        plane3->useTextureMapping = true;
//        plane3->setTextMapOfObject(txtmp3);
//        
//        plane4->useTextureMapping = true;
//        plane4->setTextMapOfObject(txtmp3);
    }
    // refraction if it's turned on
    if (REFRACTION_FLAG) {
        rt.setRefractionMode(REFRACTION_FLAG);
    }
    
	if ( rt.getEnvMapMode() != Raytracer::NONE ) {
		// load images
		EnvMap env;
		if ( _DEBUG ) {
			env = EnvMap(
                         "EnvMaps/DebugMaps/posx.bmp",
                         "EnvMaps/DebugMaps/posy.bmp",
                         "EnvMaps/DebugMaps/posz.bmp",
                         "EnvMaps/DebugMaps/negx.bmp",
                         "EnvMaps/DebugMaps/negy.bmp",
                         "EnvMaps/DebugMaps/negz.bmp"
                         );
		} else {
			env = EnvMap(
                         "EnvMaps/SaintLazarusChurch/posx.bmp",
                         "EnvMaps/SaintLazarusChurch/posy.bmp",
                         "EnvMaps/SaintLazarusChurch/posz.bmp",
                         "EnvMaps/SaintLazarusChurch/negx.bmp",
                         "EnvMaps/SaintLazarusChurch/negy.bmp",
                         "EnvMaps/SaintLazarusChurch/negz.bmp"
                         );
		}
        
		rt.setEnvMap(env);
	}
    
	printf("REFRACTION SCENE : 1 :: Rendering...\n");
	rt.render(width, height, eye2, view2, up, fov, "refraction_2.bmp");
    Point3D eye3(0, 0, 1);
    Vector3D view3(0, 0, -1);
    printf("REFRACTION SCENE : 2 :: Rendering...\n");
    rt.render(width, height, eye3, view3, up, fov, "refraction_1.bmp");
    
	printf("REFRACTION SCENE : 1 :: Done!\n");
}


/**
 * Wooden Monkey Scene 2
 */

void wmonkey_scene_2()
{
	printf("WOODEN MONKEY SCENE : 2 (path tracing) ----------------------------------\n\n");
	Raytracer rt;
	int width = 16 * 20 * 2;
	int height = 12 * 20 * 2;
    
	// room dimensions
	double wDiameter = 20;
	double wRadius = wDiameter / 2;
	double wallSize[] = {wDiameter + DBL_EPSILON, wDiameter + DBL_EPSILON, wDiameter + DBL_EPSILON};
    
	// Camera parameters.
	Point3D camera_pos(0, 0, 1);
	Vector3D camera_target(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 65;
    
	double l0c = .9;
	PointLight * light0 = new PointLight(
                                         //			Point3D(-wDiameter - roomRad + (roomRad / 2), 50 * 2, 50 * 2),
                                         Point3D(0, wRadius - 5, 0),
                                         Colour(l0c, l0c, l0c),
                                         0.2);
	rt.addLightSource(light0);
    
	// http://en.wikipedia.org/wiki/Cornell_Box
	// http://www.kevinbeason.com/smallpt/#moreinfo
	Material matLight(
                      Colour(1, 1, 1) 			// ambient
                      , Colour(1, 1, 1)			// diffuse
                      , Colour(1, 1, 1)			// spec
                      , 0, 0, 0, 0
                      );
    
	Material matMirror(
                       Colour(0.3, 0.3, 0.3) 			// ambient
                       , Colour(0.1, 0.1, 0.1)		// diffuse
                       , Colour(0.628281, 0.555802, 0.366065)		// spec
                       , 51.2, 1.0
                       );
    
	Material matGlass(
                      Colour(0, 0, 0) 			// ambient
                      , Colour(0.1, 0.1, 0.1)			// diffuse
                      , Colour(1.0, 1.0, 1.0)		// spec
                      , 100, 0.0, 1.01, 0.9
                      );
    
	Material matBeige(
                      Colour(0.607843, 0.549019, 0.372549) 			// ambient
                      , Colour(0.741176, 0.686274, 0.525490)		// diffuse
                      , Colour(0.933333, 0.901960, 0.807843)		// spec
                      , 12.8
                      );
    
	Material matReddishWall(
                            Colour(.25, .25, .25) 			// ambient
                            , Colour(.75, .25, .25)		// diffuse
                            , Colour(.3, .3, .3)		// spec
                            , 2, 0, 0, 0
                            );
    
	Material matBluishWall(
                           Colour(.25, .25, .25) 			// ambient
                           , Colour(.25, .25, .75)		// diffuse
                           , Colour(.3, .3, .3)		// spec
                           , 2, 0, 0, 0
                           );
    
	Material matGreenishWall(
                             Colour(.25, .25, .25) 			// ambient
                             , Colour(.25, .75, .25)		// diffuse
                             , Colour(.3, .3, .3)		// spec
                             , 2, 0, 0, 0
                             );
    
	Material matBaseWall(
                         Colour(.25, .25, .25) 			// ambient
                         , Colour(.25, .25, .25)		// diffuse
                         , Colour(.6, .6, .6)		// spec
                         , 2, 0, 0, 0
                         );
    
	// create and position the box
	SceneDagNode* wallLeft = rt.addObject( new UnitSquare(), &matReddishWall );
	SceneDagNode* wallRight = rt.addObject( new UnitSquare(), &matBluishWall );
	SceneDagNode* wallFront = rt.addObject( new UnitSquare(), &matBaseWall );
	SceneDagNode* wallTop = rt.addObject( new UnitSquare(), &matLight );
	SceneDagNode* wallBot = rt.addObject( new UnitSquare(), &matGreenishWall );
    
	rt.translate(wallFront, Vector3D(0, 0, -wDiameter - wRadius + DBL_EPSILON));
	rt.scale(wallFront, Point3D(0, 0, 0), wallSize);
    
	rt.translate(wallRight, Vector3D(wRadius, 0, -wDiameter + DBL_EPSILON));
	rt.rotate(wallRight, 'y', -90);
	rt.scale(wallRight, Point3D(0, 0, 0), wallSize);
    
	rt.translate(wallLeft, Vector3D(-wRadius, 0, -wDiameter + DBL_EPSILON));
	rt.rotate(wallLeft, 'y', 90);
	rt.scale(wallLeft, Point3D(0, 0, 0), wallSize);
    
	rt.translate(wallTop, Vector3D(0, wRadius, -wDiameter + DBL_EPSILON));
	rt.rotate(wallTop, 'x', 90);
	rt.scale(wallTop, Point3D(0, 0, 0), wallSize);
    
	rt.translate(wallBot, Vector3D(0, -wRadius, -wDiameter + DBL_EPSILON));
	rt.rotate(wallBot, 'x', -90);
	rt.scale(wallBot, Point3D(0, 0, 0), wallSize);
    
	// create some objects within the box...
	SceneDagNode* sphere_chrome = rt.addObject( new UnitSphere(), &matGlass );
	double _sChrome_size = 2;
	double _sChrome[] = {_sChrome_size, _sChrome_size, _sChrome_size};
	rt.translate(sphere_chrome,
                 Vector3D(wRadius / 3 * 2,
                          -wRadius + (_sChrome_size),
                          -wDiameter + (wRadius / 3)));
	rt.scale(sphere_chrome, Point3D(0, 0, 0), _sChrome);
    
	SceneDagNode* sphere_glass = rt.addObject( new UnitSphere(), &matMirror );
	double _sGlass_size = 2.5;
	double _sGlass[] = {_sGlass_size, _sGlass_size, _sGlass_size};
	rt.translate(sphere_glass,
                 Vector3D(-wRadius / 3 * 2,
                          -wRadius + (_sGlass_size),
                          -wDiameter - (wRadius / 3)));
	rt.scale(sphere_glass, Point3D(0, 0, 0), _sGlass);
    
	SceneDagNode* sphere_beige = rt.addObject( new UnitSphere(), &matBeige );
	double _sBeige_size = 3.5;
	double _sBeige[] = {_sBeige_size, _sBeige_size, _sBeige_size};
	rt.translate(sphere_beige,
                 Vector3D(wRadius / 3 * 2,
                          -wRadius + (_sGlass_size * 2),
                          -wDiameter - (wRadius / 3)));
	rt.scale(sphere_beige, Point3D(0, 0, 0), _sBeige);
    
	rt.setAAMode(Raytracer::AA_SUPER_SAMPLING);
	rt.setAAMode(Raytracer::NONE);
	rt.setShadingMode(Raytracer::SCENE_MODE_DIFFUSE);
	rt.setShadows(Raytracer::SHADOW_CAST);
	rt.setEnvMapMode(Raytracer::ENV_MAP_CUBE_SKYBOX);
    //	rt.setColorSpaceMode(Raytracer::COLOR_ENC_SRGB_GAMMA_CORRECT);
	rt.setReflDepth(4);
    
    // refraction if it's turned on
    if (REFRACTION_FLAG) {
        rt.setRefractionMode(REFRACTION_FLAG);
    }
    
	if ( rt.getEnvMapMode() != Raytracer::NONE ) {
		// load images
		EnvMap env;
		if ( _DEBUG ) {
			env = EnvMap(
                         "EnvMaps/DebugMaps/posx.bmp",
                         "EnvMaps/DebugMaps/posy.bmp",
                         "EnvMaps/DebugMaps/posz.bmp",
                         "EnvMaps/DebugMaps/negx.bmp",
                         "EnvMaps/DebugMaps/negy.bmp",
                         "EnvMaps/DebugMaps/negz.bmp"
                         );
		} else {
			env = EnvMap(
                         "EnvMaps/SaintLazarusChurch/posx.bmp",
                         "EnvMaps/SaintLazarusChurch/posy.bmp",
                         "EnvMaps/SaintLazarusChurch/posz.bmp",
                         "EnvMaps/SaintLazarusChurch/negx.bmp",
                         "EnvMaps/SaintLazarusChurch/negy.bmp",
                         "EnvMaps/SaintLazarusChurch/negz.bmp"
                         );
		}
        
		rt.setEnvMap(env);
	}
    
	printf("WOODEN MONKEY SCENE : 2 :: Rendering...\n");
	rt.render(width, height, camera_pos, camera_target, up, fov, "wmonkey_2.bmp");
	printf("WOODEN MONKEY SCENE : 2 :: Done!\n");
}
