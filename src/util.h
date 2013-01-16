/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		utility functions and structures
		(based on code from CGL, University of Waterloo),
		modify this file as you see fit.

***********************************************************/

#ifndef _UTIL_
#define _UTIL_

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <string.h>

#include "bmp_io.h"

#ifndef _DEBUG
#define _DEBUG 0
#endif

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#ifndef IOR_AIR
#define IOR_AIR	1.0
#endif

#ifndef DEFAULT_RELF
#define DEFAULT_REFL	0.0
#endif

#ifndef DEFAULT_TRANSP
#define DEFAULT_TRANSP	0.0
#endif



// small values for neg/pos floats for weird edge cases...bah
#ifndef DBL_EPSILON
#define DBL_EPSILON 0.000001
#endif

#ifndef DBL_0_UPPER
#define DBL_0_UPPER	0.000001
#endif
#ifndef DBL_0_LOWER
#define DBL_0_LOWER	-0.000001
#endif

#ifndef DBL_1_UPPER
#define DBL_1_UPPER	1.000001
#endif
#ifndef DBL_1_LOWER
#define DBL_1_LOWER	0.999999
#endif

class Point3D {
public:
	Point3D();
	Point3D(double x, double y, double z);
	Point3D(const Point3D& other);

	Point3D& operator =(const Point3D& other);
	double& operator[](int i);
	double operator[](int i) const;

private:
	double m_data[3];
};

class Vector3D {
public:
	Vector3D();
	Vector3D(double x, double y, double z);
	Vector3D(const Vector3D& other);

	Vector3D& operator =(const Vector3D& other);
	double& operator[](int i);
	double operator[](int i) const;

	double length() const;
	double normalize();
	double dot(const Vector3D& other) const;
	Vector3D cross(const Vector3D& other) const;

private:
	double m_data[3];
};

// standard operators on points and vectors
Vector3D operator *(double s, const Vector3D& v);
Vector3D operator +(const Vector3D& u, const Vector3D& v);
Point3D operator +(const Point3D& u, const Vector3D& v);
Vector3D operator -(const Point3D& u, const Point3D& v);
Vector3D operator -(const Vector3D& u, const Vector3D& v);
Vector3D operator -(const Vector3D& u);
Point3D operator -(const Point3D& u, const Vector3D& v);
Vector3D cross(const Vector3D& u, const Vector3D& v);
std::ostream& operator <<(std::ostream& o, const Point3D& p);
std::ostream& operator <<(std::ostream& o, const Vector3D& v);

class Vector4D {
public:
	Vector4D();
	Vector4D(double w, double x, double y, double z);
	Vector4D(const Vector4D& other);

	Vector4D& operator =(const Vector4D& other);
	double& operator[](int i);
	double operator[](int i) const;

private:
	double m_data[4];
};

class Matrix4x4 {
public:
  Matrix4x4();
  Matrix4x4(const Matrix4x4& other);
  Matrix4x4& operator=(const Matrix4x4& other);

  Vector4D getRow(int row) const;
  double *getRow(int row);
  Vector4D getColumn(int col) const;

  Vector4D operator[](int row) const;
  double *operator[](int row);

  Matrix4x4 transpose() const;

private:
  double m_data[16];
};

Matrix4x4 operator *(const Matrix4x4& M, const Matrix4x4& N);
Vector3D operator *(const Matrix4x4& M, const Vector3D& v);
Point3D operator *(const Matrix4x4& M, const Point3D& p);
// Multiply n by the transpose of M, useful for transforming normals.
// Recall that normals should be transformed by the inverse transpose
// of the matrix.
Vector3D transNorm(const Matrix4x4& M, const Vector3D& n);
std::ostream& operator <<(std::ostream& os, const Matrix4x4& M);

class Colour {
public:
	Colour();
	Colour(double r, double g, double b);
	Colour(const Colour& other);

	Colour& operator =(const Colour& other);
	Colour operator *(const Colour& other);
	double& operator[](int i);
	double operator[](int i) const;

	void clamp();

private:
	double m_data[3];
};

Colour operator *(double s, const Colour& c);
Colour operator +(const Colour& u, const Colour& v);
std::ostream& operator <<(std::ostream& o, const Colour& c);

struct Material;
void _init_material_IOR(Material * mat);

struct Material {
	// all settings
	Material(
			Colour ambient, Colour diffuse, Colour specular,
			double spec_exp,
			double amt_refl, double indexRefraction, double amt_transparent) :
		ambient(ambient), diffuse(diffuse), specular(specular),
		specular_exp(spec_exp),
		refl_amt(amt_refl), IOR(indexRefraction), _IOR_C(IOR_AIR / IOR),
		transparency(amt_transparent) {}

	// transparency related
	Material(
			Colour ambient, Colour diffuse, Colour specular,
			double spec_exp,
			double indexRefraction, double amt_transparent) :
		ambient(ambient), diffuse(diffuse), specular(specular),
		specular_exp(spec_exp),
		refl_amt(DEFAULT_REFL), IOR(indexRefraction), _IOR_C(IOR_AIR / IOR),
		transparency(amt_transparent) {}

	// reflection related
	Material(
			Colour ambient, Colour diffuse, Colour specular,
			double spec_exp,
			double amt_refl) :
				ambient(ambient), diffuse(diffuse), specular(specular),
				specular_exp(spec_exp),
				refl_amt(amt_refl), IOR(IOR_AIR), _IOR_C(IOR_AIR / IOR),
				transparency(DEFAULT_TRANSP) {}

	// basic phong
	Material(
			Colour ambient, Colour diffuse, Colour specular,
			double spec_exp) :
				ambient(ambient), diffuse(diffuse), specular(specular),
				specular_exp(spec_exp),
				refl_amt(DEFAULT_REFL), IOR(IOR_AIR), _IOR_C(IOR_AIR / IOR),
				transparency(DEFAULT_TRANSP)  {}
    
	// Ambient components for Phong shading.
	Colour ambient;

	// Diffuse components for Phong shading.
	Colour diffuse;

	// Specular components for Phong shading.
	Colour specular;

	// Specular expoent.
	double specular_exp;

	// reflection amount [0, 1] (0 = none)
	double refl_amt;

	// index of refraction
	double IOR;

	// value of C for this material, assuming c_air = IOR_AIR
	double _IOR_C;

	// amount of light to let through [0, 1] (0 = none)
	double transparency;
};

/**
 * A BMP image
 */
struct _BMP {
	_BMP(std::string file_in) :
		file_in_name(file_in)
	{
		bmp_read_status = false;
		width = 0;
		height = 0;
		rarray = NULL;
		garray = NULL;
		barray = NULL;

		bmp_read_status = bmp_read(
				(char *) file_in_name.c_str(),
				&width, &height,
				&rarray, &garray, &barray
				);

		if ( bmp_read_status ) {
			printf("An error occured while reading BMP file: %s.\n", file_in_name.c_str());
		} else {
			printf("Loaded file: %s.\n", file_in_name.c_str());
		}
	}

	_BMP() {
		bmp_read_status = false;
		width = 0;
		height = 0;
		rarray = NULL;
		garray = NULL;
		barray = NULL;
	}

	// true iff an error occured
	bool bmp_read_status;

	unsigned long int width;
	long int height;

	unsigned char * rarray;
	unsigned char * garray;
	unsigned char * barray;

	std::string file_in_name;
};

/**
 * An environment map (6 maps, one per cube face)
 */
struct EnvMap {
	EnvMap(
		std::string px,
		std::string py,
		std::string pz,
		std::string nx,
		std::string ny,
		std::string nz) {

		printf("Loading environment maps...\n");

		_px = new _BMP(px);
		_py = new _BMP(py);
		_pz = new _BMP(pz);
		_nx = new _BMP(nx);
		_ny = new _BMP(ny);
		_nz = new _BMP(nz);

		printf("...done.\n");
	}

	EnvMap() {
		_px = new _BMP();
		_py = new _BMP();
		_pz = new _BMP();
		_nx = new _BMP();
		_ny = new _BMP();
		_nz = new _BMP();
	}

	Colour getColorForRay(Vector3D rayDir) {
		double x,y,z;
		double u = 0, v = 0;

		double x1 = rayDir[0];
		double y1 = rayDir[1];
		double z1 = rayDir[2];

		x = fabs(x1);
		y = fabs(y1);
		z = fabs(z1);

		_BMP * side = _px;
		if ( x > y && x > z ) {
			u = (y1 + x1) / (2 * x1);
			v = 1 - ((z1 + x1) / (2 * x1));

			if ( x1 > 0 ) { // +x
				side = _px;
			} else { // -x
				u = 1 - u;
				side = _nx;
			}
		} else if ( y > x && y > z ) {
			u = 1 - ((z1 + y1) / (2 * y1));
			v = (x1 + y1) / (2 * y1);

			if ( y1 > 0 ) { // +y
				side = _py;
			} else { // -y
				v = 1 - v;
				side = _ny;
			}
		} else if ( z > x && z > y ) {
			u = (y1 + z1) / (2 * z1);
			v = (x1 + z1) / (2 * z1);

			if ( z1 > 0 ) { // +z
				side = _pz;
			} else { // -z
				u = 1 - u;
				side = _nz;
			}
		}

		if ( u > 1 ) u = 1;
		else if ( u < 0 ) u = 0;
		if ( v > 1 ) v = 1;
		else if ( v < 0 ) v = 0;

		return getColorAtUVonSide(int(u * side->width), int(v * side->height), side);
	}

	Colour getColorAtUVonSide(int u, int v, _BMP * side) {
		int r = side->rarray[u * side->width + v];
		int g = side->garray[u * side->width + v];
		int b = side->barray[u * side->width + v];

		double rd = r / 255.0f;
		double gd = g / 255.0f;
		double bd = b / 255.0f;

		return Colour(rd, gd, bd);
	}

	// each side of the cube
	struct _BMP * _px;
	struct _BMP * _py;
	struct _BMP * _pz;
	struct _BMP * _nx;
	struct _BMP * _ny;
	struct _BMP * _nz;
};

struct Intersection {
	// Location of intersection.
	Point3D point;
	// Normal at the intersection.
	Vector3D normal;
	// Material at the intersection.
	Material* mat;
	// Position of the intersection point on your ray.
	// (i.e. point = ray.origin + t_value * ray.dir)
	// This is used when you need to intersect multiply objects and
	// only want to keep the nearest intersection.
	double t_value;
	// Set to true when no intersection has occured.
	bool none;
};


/*
 Texture Map structure
 This data structure holds the texture map to be used for different scene objects
*/
struct TextureMap {
    TextureMap(std::string txtImg){
		printf("Loading Texture map...\n");
		_txtImg = new _BMP(txtImg);
        
		printf("...done.\n");
	}
    
    TextureMap(){
        _txtImg = new _BMP();
    } 
    
    struct _BMP *_txtImg;
    
    Colour getColourAtPxLocInImg( double u_px, double v_px){
        //make sure (u, v) passed in is
        int u, v;
        u = _txtImg->width*u_px;
        v = _txtImg->width*v_px;
        
        if(_DEBUG){
            printf(" ~~~~ computer pixel intensity at (u, v)  = (%d, %d \n", u, v);
        }

        int r = _txtImg->rarray[u *_txtImg->width + v];
		int g = _txtImg->garray[u *_txtImg->width + v];
		int b = _txtImg->barray[u *_txtImg->width + v];
        
		double rd = r/ 255.0f;
		double gd = g/ 255.0f;
		double bd = b/ 255.0f;
		return Colour(rd, gd, bd);
    }

    
};


// Ray structure.
struct Ray3D {
	Ray3D() : _cur_speed(IOR_AIR), x(0), y(0) {
		intersection.none = true;
		shadowed = false;
	}

	Ray3D( Point3D p, Vector3D v ) : origin(p), dir(v), _cur_speed(IOR_AIR), x(0), y(0) {
		intersection.none = true;
		shadowed = false;
	}

	// Origin and direction of the ray.
	Point3D origin;
	Vector3D dir;

	// Intersection status, should be computed by the intersection
	// function.
	Intersection intersection;

	// Current colour of the ray, should be computed by the shading
	// function.
	Colour col;

	// Set to true if we've hit a shadowed area
	bool shadowed;

	// Current light speed of the ray (IOR_AIR by default)
	double _cur_speed;

	int x;
	int y;
};
#endif





