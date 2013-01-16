/***********************************************************
	 Starter code for Assignment 3

	 This code was originally written by Jack Wang for
			CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

/**
 * Return true iff the given t_value is closer to the eye than the value that
 * already exists within the ray, if any.
 */
bool isCloser(Ray3D& ray, double t_value) {
	return
			(ray.intersection.none || t_value <= ray.intersection.t_value)
			&& t_value >= DBL_0_UPPER;
}


/**
 * A unit square (plane) centered at origin, on the x/y plane.
 */
bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// http://en.wikipedia.org/wiki/Line-plane_intersection#Algebraic_form
	Point3D p0(0, 0, 0); // a point on the plane
	Vector3D n(0, 0, 1); // normal
	double sMax = 0.5; // distance from one side to another
	double sMin = -0.5; // distance from one side to another

	Vector3D d = worldToModel * ray.dir;
	Point3D rayOrigin_model = worldToModel * ray.origin;

	double t = 0, numr, denom;

	// compute the intersection point
	numr = (p0 - rayOrigin_model).dot(n);
	denom = d.dot(n);

	if ( numr != 0 && denom != 0 ) {
		t = numr / denom;
		Point3D intPoint = rayOrigin_model + t * d;

		if ( isCloser(ray, t) ) {
			// check if we're within the bounds of the square
			if (
					intPoint[0] < sMax && intPoint[0] > sMin
					&& intPoint[1] < sMax && intPoint[1] > sMin
					) {

				ray.intersection.point = modelToWorld * intPoint;
				ray.intersection.normal = transNorm(worldToModel, n);
				ray.intersection.normal.normalize();
				ray.intersection.t_value = t;
				ray.intersection.none = false;

				return true;
			}
		}
	}

	return false;
}

/**
 * A unit sphere centered at origin.
 *
 * Goal here is to fill ray.intersection with correct values should an
 * intersection occur. This includes intersection.point, intersection.normal,
 * intersection.none, intersection.t_value.
 */
bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// http://www.cs.brown.edu/~ls/teaching08/LS14_RayTracingIntersections.pdf

	Point3D sphereOrigin(0, 0, 0);

	// set to obj space
	Vector3D d = worldToModel * ray.dir;
	Point3D rayOrigin_model = worldToModel * ray.origin;
	Vector3D a = rayOrigin_model - sphereOrigin;

	double A, B, C, D;
	A = d.dot(d);
	B = a.dot(d);
	C = a.dot(a) - 1;
	D = B * B - A * C;

	// D < 0: implicitly -> no hit
	if ( D >= 0 ) {
		double x, y, lambda_s_1, lambda_s_2, t = 0;
		x = (-B / A);
		y = (sqrt(D) / A);
		lambda_s_1 = x + y;

		if ( D == 0 || D <= DBL_0_UPPER ) {
			// 1 hit at r(lambda_s_1)
			t = lambda_s_1;
		} else {
			// 2 hits at r(lambda_s_1), r(lambda_s_2)
			lambda_s_2 = x - y;
			if ( abs(lambda_s_1) > 0 ) {
				if ( abs(lambda_s_2) > 0 ) {
					// sphere is in front of us, both hits valid
					t = lambda_s_2;
				} else {
					// we're inside the sphere, first hit valid
					t = lambda_s_1;
				}
			} else if ( lambda_s_2 < DBL_0_LOWER ) {
				// the sphere is behind us so we ignore it...
				return false;
			}
		}

		// we've got a hit, time to calculate the rest of the data
		if ( isCloser(ray, t) ) {
			Point3D intPoint = rayOrigin_model + t * d; // r(t)
			Vector3D normal = intPoint - sphereOrigin;

//			std::cout << "ip" << intPoint << " = " << rayOrigin_model << " + " << t << " * " << d << std::endl;
            //if ray intersects, change the material to that of the texture map
            
			ray.intersection.point = modelToWorld * intPoint;
			ray.intersection.normal = transNorm(worldToModel, normal);
			ray.intersection.normal.normalize();
			ray.intersection.t_value = t;
			ray.intersection.none = false;
            
		} else {
			return false;
		}

		return true;   
	}

	return false;
}

//Implement the functions that generates the colour from the texture map of the scene object

Colour UnitSphere::computeColourAtImagePoint(Point3D imgPoint, const Matrix4x4& worldToModel){
    
    //convert the image point to point in model space
    Point3D pointOnSphere = worldToModel*imgPoint;
    
    // compute the parametric alpha, beta coordinate of the unit sphere
    double x = pointOnSphere[0];
    double y = pointOnSphere[1];
    double z = pointOnSphere[2];
    
    double alpha = atan(y/x);
    double beta = acos(z);
    
    if (alpha < 0) {
        // we want alpha to be [0. 2PI)
        alpha = 2*M_PI + alpha;
    }
    
    if (beta < 0) {
        // we want beta to be [
        printf("WARNING!!! ---- beta should not be negative: beta = %f \n", beta);
        exit(0);
    }
    
    //the corresponding coordinate (u, v) on the texture map is
    double u = alpha/(2*M_PI);
    double v = beta/M_PI;
    
    Colour txtMapCl = _myTexture.getColourAtPxLocInImg(u, v);
  
    if (_DEBUG) {
          printf(" >>>    Texture Value at ( x, y, z) = ( %f, %f, %f), ( u, v) = ( %f, %f), Colour = ( %f, %f, %f) \n", x, y, z, u, v, txtMapCl[0], txtMapCl[1], txtMapCl[2]);  
    }

    txtMapCl.clamp(); //clamp colour to make sure it's within [0,1]
    
    return txtMapCl;
}

Colour UnitSquare::computeColourAtImagePoint(Point3D imgPoint, const Matrix4x4& worldToModel){
    
    Point3D pointOnSquare = worldToModel*imgPoint;
    
    //compute the parameteric (u, v) for the plane
    //the normal of the plane is (0, 0, 1), hence it's in x-y plane
    double u = pointOnSquare[0];
    double v = pointOnSquare[1];
    
    //ofset (u, v) by 0.5 so that (0,0) is at the corner of the image
    u = u + 0.5;
    v = v+ 0.5;
    
    Colour txMapCl = _myTexture.getColourAtPxLocInImg(u, v);
    
    return txMapCl;
}
