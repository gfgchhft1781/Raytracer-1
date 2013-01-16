/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		classes defining primitives in the scene

***********************************************************/

#include "util.h"

// All primitives should provide a intersection function.  
// To create more primitives, inherit from SceneObject.
// Namely, you can create, Sphere, Cylinder, etc... classes
// here.
class SceneObject {
public:
	// Returns true if an intersection occured, false otherwise.
	virtual bool intersect( Ray3D&, const Matrix4x4&, const Matrix4x4& ) = 0;
    
    //implement this function for objects that will have texture
    //this is used to set the texture of the object
    virtual void setTextureMap(TextureMap mp) = 0;
    
    //implement the following function if texture map is being used for the scene object
    virtual Colour computeColourAtImagePoint(Point3D imgPoint, const Matrix4x4& worldToModel) = 0;
};

// Example primitive you can create, this is a unit square on 
// the xy-plane.
class UnitSquare : public SceneObject {
public:
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
    
    Colour computeColourAtImagePoint(Point3D imgPoint, const Matrix4x4& worldToModel);
    
    void setTextureMap(TextureMap mp){
        _myTexture = mp;
        _useTextureMap = true;
    }
    
private:
    TextureMap _myTexture;
    bool _useTextureMap;
};

class UnitSphere : public SceneObject {
public:
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
    
    Colour computeColourAtImagePoint(Point3D imgPoint, const Matrix4x4& worldToModel);
    
    void setTextureMap(TextureMap mp){
        _myTexture = mp;
        _useTextureMap = true;
    }
private:
    TextureMap _myTexture;
    bool _useTextureMap;
};

