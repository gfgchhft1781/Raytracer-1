/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		light source classes

***********************************************************/

#include "util.h"

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
	virtual void shade( Ray3D& ) = 0;
	virtual void shade_signature( Ray3D& ray ) = 0;
	virtual void shade_diffuse( Ray3D& ray ) = 0;
	virtual Point3D get_position() const = 0; 
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col), _shadow_opacity(0) {}

	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular), _shadow_opacity(0) {}

	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular, double shadow_opacity )
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse),
	_col_specular(specular), _shadow_opacity(shadow_opacity) {}

	PointLight( Point3D pos, Colour col, double shadow_opacity ) : _pos(pos), _col_ambient(col),
	_col_diffuse(col), _col_specular(col), _shadow_opacity(shadow_opacity) {}

	void shade( Ray3D& ray );
	void shade_signature( Ray3D& ray );
	void shade_diffuse( Ray3D& ray );
	Point3D get_position() const { return _pos; }
	
	void setAmbient(Colour col);

private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 

	// how opaque a shadow is [0,1], with 1 = fully transparent
	double _shadow_opacity;

	void do_diffuse(const Vector3D& s, const Vector3D& n, Ray3D& ray);
};
