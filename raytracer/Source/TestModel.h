#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::distance;

// Used to describe a light
class PhongLightSource {
public:
  vec4 position;
  vec3 color;
	float ambientIntensity;
	float diffuseIntensity;
	float specularIntensity;

  PhongLightSource(const vec4 &p, const vec3 &c, const float &am, const float &di, const float &sp)
		: position(p), color(c), ambientIntensity(am), diffuseIntensity(di), specularIntensity(sp)
  {

	}
};

class PhongMaterial {
public:
	vec3 color;
	float ambientRef;
	float diffuseRef;
	float specularRef;
	float shininess;

	PhongMaterial(const vec3 &c, const float &am, const float &di, const float &sp, const float &sh)
		: color(c), ambientRef(am), diffuseRef(di), specularRef(sp), shininess(sh)
  {

	}
};

// Used to describe a triangular surface:
class PhongTriangle {
public:
	vec4 v0;
	vec4 v1;
	vec4 v2;
	vec4 normal;
	PhongMaterial material;

	PhongTriangle( vec4 v0, vec4 v1, vec4 v2, PhongMaterial material )
		: v0(v0), v1(v1), v2(v2), material(material)
	{
		ComputeNormal();
	}

	void ComputeNormal()
	{
	  vec3 e1 = glm::vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
	  vec3 e2 = glm::vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
	  vec3 normal3 = glm::normalize( glm::cross( e2, e1 ) );
	  normal.x = normal3.x;
	  normal.y = normal3.y;
	  normal.z = normal3.z;
	  normal.w = 1.0;
	}
};

// Used to describe a spherical surface:
class PhongSphere {
public:
	glm::vec4 centre;
	float radius;
	PhongMaterial material;

	PhongSphere( glm::vec4 centre, float radius, PhongMaterial material )
		: centre(centre), radius(radius), material(material)
	{

	}
};

class LightSource {
public:
	float watts;
	vec3 color;
	vec4 position;
	vec4 direction;
	float width;
	float length;
	float height;

	LightSource( float watts, vec3 color, vec4 position, vec4 direction, float width, float length, float height)
		: watts(watts), color(color), position(position), direction(direction), width(width), length(length), height(height)
	{

	}
};

class Material {
public:
	vec3 diffuseRef;
	vec3 specRef;

	Material( vec3 diffuseRef, vec3 specRef)
		: diffuseRef(diffuseRef), specRef(specRef)
	{

	}
};

// Used to describe a triangular surface:
class Triangle {
public:
	vec4 v0;
	vec4 v1;
	vec4 v2;
	vec4 normal;
	vec3 color;
	Material material;

	Triangle( vec4 v0, vec4 v1, vec4 v2, vec3 color, Material material )
		: v0(v0), v1(v1), v2(v2), color(color), material(material)
	{
		ComputeNormal();
	}

	void ComputeNormal()
	{
	  vec3 e1 = glm::vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
	  vec3 e2 = glm::vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
	  vec3 normal3 = glm::normalize( glm::cross( e2, e1 ) );
	  normal.x = normal3.x;
	  normal.y = normal3.y;
	  normal.z = normal3.z;
	  normal.w = 1.0;
	}
};

// Used to describe a spherical surface:
class Sphere {
public:
	vec4 centre;
	float radius;
	vec3 color;
	Material material;

	Sphere( vec4 centre, float radius, vec3 color, Material material )
		: centre(centre), radius(radius), color(color), material(material)
	{

	}
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles, std::vector<Sphere>& spheres, std::vector<LightSource>& lights)
{
	using glm::vec3;
	using glm::vec4;

	// Defines colors:
	vec3 red(    0.75f, 0.15f, 0.15f );
	vec3 yellow( 0.75f, 0.75f, 0.15f );
	vec3 green(  0.15f, 0.75f, 0.15f );
	vec3 cyan(   0.15f, 0.75f, 0.75f );
	vec3 blue(   0.15f, 0.15f, 0.75f );
	vec3 purple( 0.75f, 0.15f, 0.75f );
	vec3 white(  0.75f, 0.75f, 0.75f );

	float matteDiffuseRef = 0.2f;

	// Define materials
	Material matteWhite( white * matteDiffuseRef, vec3(0,0,0));
	Material matteRed( red * matteDiffuseRef, vec3(0,0,0));
	Material matteBlue( blue * matteDiffuseRef, vec3(0,0,0));
	Material matteGreen( green * matteDiffuseRef, vec3(0,0,0));
	Material matteYellow( yellow * matteDiffuseRef, vec3(0,0,0));
	Material mattePurple( purple * matteDiffuseRef, vec3(0,0,0));
	Material matteCyan( cyan * matteDiffuseRef, vec3(0,0,0));

	Material mirror( vec3(0,0,0), vec3(1,1,1));

	lights.clear();
	lights.reserve( 1 );

	triangles.clear();
	triangles.reserve( 5*2*3 );

	spheres.clear();
	spheres.reserve( 1 );

	// ---------------------------------------------------------------------------
	// Room

	float L = 555;			// Length of Cornell Box side.

	vec4 A(L,0,0,1);
	vec4 B(0,0,0,1);
	vec4 C(L,0,L,1);
	vec4 D(0,0,L,1);

	vec4 E(L,L,0,1);
	vec4 F(0,L,0,1);
	vec4 G(L,L,L,1);
	vec4 H(0,L,L,1);

	lights.push_back( LightSource( 30, vec3(1, 1, 1), vec4(0, -1.0, -0.5, 1), vec4(0, -1, 0, 1), 0.35, 0.35, 0.2) );

	spheres.push_back( Sphere(vec4(0.4,0,-0.2,1), 0.2, white, matteWhite) );

	// // Floor:
	// triangles.push_back( Triangle( C, B, A, yellow, matteYellow ) );
	// triangles.push_back( Triangle( C, D, B, yellow, matteYellow ) );
	//
	// // Left wall
	// triangles.push_back( Triangle( A, E, C, red, matteRed ) );
	// triangles.push_back( Triangle( C, E, G, red, matteRed ) );
	//
	// // Right wall
	// triangles.push_back( Triangle( F, B, D, blue, matteBlue ) );
	// triangles.push_back( Triangle( H, F, D, blue, matteBlue ) );
	//
	// // Ceiling
	// triangles.push_back( Triangle( E, F, G, white, matteWhite ) );
	// triangles.push_back( Triangle( F, H, G, white, matteWhite ) );
	//
	// // Back wall
	// triangles.push_back( Triangle( G, D, C, cyan, matteCyan ) );
	// triangles.push_back( Triangle( G, H, D, cyan, matteCyan ) );

	// Floor:
	triangles.push_back( Triangle( C, B, A, white, matteWhite ) );
	triangles.push_back( Triangle( C, D, B, white, matteWhite ) );

	// Left wall
	triangles.push_back( Triangle( A, E, C, red, matteRed ) );
	triangles.push_back( Triangle( C, E, G, red, matteRed ) );

	// Right wall
	triangles.push_back( Triangle( F, B, D, blue, matteBlue ) );
	triangles.push_back( Triangle( H, F, D, blue, matteBlue ) );

	// Ceiling
	triangles.push_back( Triangle( E, F, G, white, matteWhite ) );
	triangles.push_back( Triangle( F, H, G, white, matteWhite ) );

	// Back wall
	triangles.push_back( Triangle( G, D, C, white, matteWhite ) );
	triangles.push_back( Triangle( G, H, D, white, matteWhite ) );

	// ---------------------------------------------------------------------------
	// Short block

	A = vec4(290,0,114,1);
	B = vec4(130,0, 65,1);
	C = vec4(240,0,272,1);
	D = vec4( 82,0,225,1);

	E = vec4(290,165,114,1);
	F = vec4(130,165, 65,1);
	G = vec4(240,165,272,1);
	H = vec4( 82,165,225,1);

	// Front
	triangles.push_back( Triangle(E,B,A,white,matteWhite) );
	triangles.push_back( Triangle(E,F,B,white,matteWhite) );

	// Front
	triangles.push_back( Triangle(F,D,B,white,matteWhite) );
	triangles.push_back( Triangle(F,H,D,white,matteWhite) );

	// BACK
	triangles.push_back( Triangle(H,C,D,white,matteWhite) );
	triangles.push_back( Triangle(H,G,C,white,matteWhite) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,white,matteWhite) );
	triangles.push_back( Triangle(E,A,C,white,matteWhite) );

	// TOP
	triangles.push_back( Triangle(G,F,E,white,matteWhite) );
	triangles.push_back( Triangle(G,H,F,white,matteWhite) );

	// ---------------------------------------------------------------------------
	// Tall block

	A = vec4(423,0,247,1);
	B = vec4(265,0,296,1);
	C = vec4(472,0,406,1);
	D = vec4(314,0,456,1);

	E = vec4(423,330,247,1);
	F = vec4(265,330,296,1);
	G = vec4(472,330,406,1);
	H = vec4(314,330,456,1);

	// // Front
	// triangles.push_back( Triangle(E,B,A,purple,mattePurple) );
	// triangles.push_back( Triangle(E,F,B,purple,mattePurple) );
	//
	// // Front
	// triangles.push_back( Triangle(F,D,B,purple,mattePurple) );
	// triangles.push_back( Triangle(F,H,D,purple,mattePurple) );
	//
	// // BACK
	// triangles.push_back( Triangle(H,C,D,purple,mattePurple) );
	// triangles.push_back( Triangle(H,G,C,purple,mattePurple) );
	//
	// // LEFT
	// triangles.push_back( Triangle(G,E,C,purple,mattePurple) );
	// triangles.push_back( Triangle(E,A,C,purple,mattePurple) );
	//
	// // TOP
	// triangles.push_back( Triangle(G,F,E,purple,mattePurple) );
	// triangles.push_back( Triangle(G,H,F,purple,mattePurple) );

	// Front
	triangles.push_back( Triangle(E,B,A,white,matteWhite) );
	triangles.push_back( Triangle(E,F,B,white,matteWhite) );

	// Front
	triangles.push_back( Triangle(F,D,B,white,matteWhite) );
	triangles.push_back( Triangle(F,H,D,white,matteWhite) );

	// BACK
	triangles.push_back( Triangle(H,C,D,white,matteWhite) );
	triangles.push_back( Triangle(H,G,C,white,matteWhite) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,white,matteWhite) );
	triangles.push_back( Triangle(E,A,C,white,matteWhite) );

	// TOP
	triangles.push_back( Triangle(G,F,E,white,matteWhite) );
	triangles.push_back( Triangle(G,H,F,white,matteWhite) );

	// ----------------------------------------------
	// Scale to the volume [-1,1]^3
	// for( size_t i=0; i<lights.size(); ++i )
	// {
	// 	lights[i].position *= 2/L;
	//
	// 	lights[i].position -= vec4(1,1,1,1);
	//
	// 	lights[i].position.x *= -1;
	// 	lights[i].position.y *= -1;
	// 	lights[i].position.w = 1.0;
	// }

	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0 *= 2/L;
		triangles[i].v1 *= 2/L;
		triangles[i].v2 *= 2/L;

		triangles[i].v0 -= vec4(1,1,1,1);
		triangles[i].v1 -= vec4(1,1,1,1);
		triangles[i].v2 -= vec4(1,1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].v0.w = 1.0;
		triangles[i].v1.w = 1.0;
		triangles[i].v2.w = 1.0;

		triangles[i].ComputeNormal();
	}
}

void LoadTestModelPhong( std::vector<PhongTriangle>& triangles, std::vector<PhongSphere>& spheres, std::vector<PhongLightSource>& lights )
{
	using glm::vec3;
	using glm::vec4;

	// Defines colors:
	vec3 red(    0.75f, 0.15f, 0.15f );
	vec3 yellow( 0.75f, 0.75f, 0.15f );
	vec3 green(  0.15f, 0.75f, 0.15f );
	vec3 cyan(   0.15f, 0.75f, 0.75f );
	vec3 blue(   0.15f, 0.15f, 0.75f );
	vec3 purple( 0.75f, 0.15f, 0.75f );
	vec3 white(  0.75f, 0.75f, 0.75f );

  vec3 darkPurple(0.65f, 0.1f, 0.65f);

	PhongMaterial matteRed(red, 1, 2, 1, 1.1);
	PhongMaterial matteYellow(yellow, 1, 2, 1, 1.1);
	PhongMaterial matteGreen(green, 1, 2, 1, 1.1);
	PhongMaterial matteCyan(cyan, 1, 2, 1, 1.1);
	PhongMaterial matteBlue(blue, 1, 2, 1, 1.1);
	PhongMaterial mattePurple(purple, 1, 2, 1, 1.1);
	PhongMaterial matteWhite(white, 1, 2, 1, 1.1);

	PhongMaterial shinyPurple(darkPurple, 4, 4, 5, 50);
	PhongMaterial shinyRed(red, 4, 4, 5, 50);
	PhongMaterial shinyYellow(yellow, 4, 4, 5, 50);

	triangles.clear();
	triangles.reserve( 5*2*3 );

	spheres.clear();
	spheres.reserve(2);

	lights.clear();
	lights.reserve(1);

	// ---------------------------------------------------------------------------
	// Triangles
	// ---------------------------------------------------------------------------
	// Room

	float L = 555;			// Length of Cornell Box side.

	vec4 A(L,0,0,1);
	vec4 B(0,0,0,1);
	vec4 C(L,0,L,1);
	vec4 D(0,0,L,1);

	vec4 E(L,L,0,1);
	vec4 F(0,L,0,1);
	vec4 G(L,L,L,1);
	vec4 H(0,L,L,1);

	// Floor:
	triangles.push_back( PhongTriangle( C, B, A, matteGreen ) );
	triangles.push_back( PhongTriangle( C, D, B, matteGreen ) );

	// Left wall
	triangles.push_back( PhongTriangle( A, E, C, mattePurple ) );
	triangles.push_back( PhongTriangle( C, E, G, mattePurple ) );

	// Right wall
	triangles.push_back( PhongTriangle( F, B, D, matteYellow ) );
	triangles.push_back( PhongTriangle( H, F, D, matteYellow ) );

	// Ceiling
	triangles.push_back( PhongTriangle( E, F, G, matteCyan ) );
	triangles.push_back( PhongTriangle( F, H, G, matteCyan ) );

	// Back wall
	triangles.push_back( PhongTriangle( G, D, C, matteWhite ) );
	triangles.push_back( PhongTriangle( G, H, D, matteWhite ) );

	// ---------------------------------------------------------------------------
	// Short block

	A = vec4(290,0,114,1);
	B = vec4(130,0, 65,1);
	C = vec4(240,0,272,1);
	D = vec4( 82,0,225,1);

	E = vec4(290,165,114,1);
	F = vec4(130,165, 65,1);
	G = vec4(240,165,272,1);
	H = vec4( 82,165,225,1);

	// Front
	triangles.push_back( PhongTriangle(E,B,A,matteRed) );
	triangles.push_back( PhongTriangle(E,F,B,matteRed) );

	// Front
	triangles.push_back( PhongTriangle(F,D,B,matteRed) );
	triangles.push_back( PhongTriangle(F,H,D,matteRed) );

	// BACK
	triangles.push_back( PhongTriangle(H,C,D,matteRed) );
	triangles.push_back( PhongTriangle(H,G,C,matteRed) );

	// LEFT
	triangles.push_back( PhongTriangle(G,E,C,matteRed) );
	triangles.push_back( PhongTriangle(E,A,C,matteRed) );

	// TOP
	triangles.push_back( PhongTriangle(G,F,E,matteRed) );
	triangles.push_back( PhongTriangle(G,H,F,matteRed) );

	// ---------------------------------------------------------------------------
	// Tall block

	A = vec4(423,0,247,1);
	B = vec4(265,0,296,1);
	C = vec4(472,0,406,1);
	D = vec4(314,0,456,1);

	E = vec4(423,330,247,1);
	F = vec4(265,330,296,1);
	G = vec4(472,330,406,1);
	H = vec4(314,330,456,1);

	// Front
	triangles.push_back( PhongTriangle(E,B,A,matteBlue) );
	triangles.push_back( PhongTriangle(E,F,B,matteBlue) );

	// Front
	triangles.push_back( PhongTriangle(F,D,B,matteBlue) );
	triangles.push_back( PhongTriangle(F,H,D,matteBlue) );

	// BACK
	triangles.push_back( PhongTriangle(H,C,D,matteBlue) );
	triangles.push_back( PhongTriangle(H,G,C,matteBlue) );

	// LEFT
	triangles.push_back( PhongTriangle(G,E,C,matteBlue) );
	triangles.push_back( PhongTriangle(E,A,C,matteBlue) );

	// TOP
	triangles.push_back( PhongTriangle(G,F,E,matteBlue) );
	triangles.push_back( PhongTriangle(G,H,F,matteBlue) );

	// ---------------------------------------------------------------------------
	// Spheres
	// ---------------------------------------------------------------------------

	spheres.push_back( PhongSphere(vec4(0.4,0,-0.2,1), 0.3, shinyPurple) );
	spheres.push_back( PhongSphere(vec4(-0.4,0.7,-0.8,1), 0.2, shinyYellow) );

	// ---------------------------------------------------------------------------
	// Lights
	// ---------------------------------------------------------------------------

	//lights.push_back( Light(vec4(0.5,-0.5,-0.7,1.0), vec3(1,1,1), 0.02f, 0.3f, 1.0f ));
	//lights.push_back( Light(vec4(-0.5,-0.5,-0.7,1.0), vec3(1,1,1), 0.02f, 0.3f, 1.0f ));
	lights.push_back( PhongLightSource(vec4(0.0,-0.5,-0.7,1.0), vec3(1.0,1.0,1.0), 0.02f, 0.3f, 1.0f ));

	// ----------------------------------------------
	// Scale to the volume [-1,1]^3

	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0 *= 2/L;
		triangles[i].v1 *= 2/L;
		triangles[i].v2 *= 2/L;

		triangles[i].v0 -= vec4(1,1,1,1);
		triangles[i].v1 -= vec4(1,1,1,1);
		triangles[i].v2 -= vec4(1,1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].v0.w = 1.0;
		triangles[i].v1.w = 1.0;
		triangles[i].v2.w = 1.0;

		triangles[i].ComputeNormal();
	}
}


#endif
