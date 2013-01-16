/***********************************************************
	 Starter code for Assignment 3

	 This code was originally written by Jack Wang for
			CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

//#include <cmath>
#include "light_source.h"

/**
 * Phong Shading
 */
void PointLight::shade( Ray3D& ray ) {
	Intersection intPoint = ray.intersection;

	// construct vectors
	Vector3D n = intPoint.normal;
	Vector3D de = -ray.dir;
	Vector3D s = get_position() - intPoint.point; // light source
	Vector3D m = ((2 * n.dot(s)) * n) - s; // perfect mirror directions

	// normalize
	n.normalize();
	s.normalize();
	m.normalize();
	de.normalize();

	// do the diffuse shading
	do_diffuse(s, n, ray);

	// do the specular shading
	Colour specular(0, 0, 0);

	double mdde = m.dot(de);
	if ( mdde >= 0 ) {
		specular =
			pow(mdde, intPoint.mat->specular_exp)
			* intPoint.mat->specular * _col_specular;
	}

	if ( ray.shadowed ) {
		ray.col = ray.col + (_shadow_opacity * specular);
	} else {
		ray.col = ray.col + specular;
	}
}

/**
 * Signature Shading
 */
void PointLight::shade_signature( Ray3D& ray ) {
	// basically, if we hit an object, set the color to its diffuse
	ray.col = ray.intersection.mat->diffuse;
	ray.col.clamp();
}

/**
 * Diffuse Shading
 */
void PointLight::shade_diffuse( Ray3D& ray ) {
	Intersection intPoint = ray.intersection;

	// construct vectors
	Vector3D n = intPoint.normal;
	Vector3D s = get_position() - intPoint.point;

	// normalize
	n.normalize();
	s.normalize();

	do_diffuse(s, n, ray);
}

void PointLight::setAmbient(Colour col) {
	_col_ambient = col;

	printf("setAmbient[%f, %f, %f]\n",
			_col_ambient[0], _col_ambient[1], _col_ambient[2]);
}

/**
 * Calculate the diffuse value for a given ray, vector s and normal n
 */
void PointLight::do_diffuse(const Vector3D& s, const Vector3D& n, Ray3D& ray)
{
	// http://www.cdf.toronto.edu/~moore/csc418/Notes/Phong.pdf
	Colour ambient, diffuse;

	// ambient is simply mat * light
	ambient = ray.intersection.mat->ambient * _col_ambient;

	// for the diffuse light, depends on angle of light
	double sdn = s.dot(n);
	if (sdn >= 0)
		diffuse = sdn * ray.intersection.mat->diffuse * _col_diffuse;

	if ( ray.shadowed )
		ray.col = ray.col + (_shadow_opacity * diffuse) + ambient;
	else
		ray.col = ray.col + diffuse + ambient;
}
