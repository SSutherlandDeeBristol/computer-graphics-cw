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
  vec3 position;
  vec3 color;
	float ambientIntensity;
	float diffuseIntensity;
	float specularIntensity;

  PhongLightSource(const vec3 &p, const vec3 &c, const float &am, const float &di, const float &sp)
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
	vec3 v0;
	vec3 v1;
	vec3 v2;
	vec3 normal;
	PhongMaterial material;

	PhongTriangle( vec3 v0, vec3 v1, vec3 v2, PhongMaterial material )
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
	}
};

// Used to describe a spherical surface:
class PhongSphere {
public:
	vec3 centre;
	float radius;
	PhongMaterial material;

	PhongSphere( vec3 centre, float radius, PhongMaterial material )
		: centre(centre), radius(radius), material(material)
	{

	}
};

class LightSource {
public:
	float watts;
	vec3 color;
	vec3 position;
	vec3 direction;
	float width;
	float length;

	LightSource( float watts, vec3 color, vec3 position, vec3 direction, float width, float length)
		: watts(watts), color(color), position(position), direction(direction), width(width), length(length)
	{

	}
};

class Material {
public:
	vec3 diffuseRef;
	vec3 specRef;
  float refractiveIndex;

	Material( vec3 diffuseRef, vec3 specRef, float refractiveIndex)
		: diffuseRef(diffuseRef), specRef(specRef), refractiveIndex(refractiveIndex)
	{

	}
};

// Used to describe a triangular surface:
class Triangle {
public:
	vec3 v0;
	vec3 v1;
	vec3 v2;
	vec3 normal;
	vec3 color;
	Material material;

	Triangle( vec3 v0, vec3 v1, vec3 v2, vec3 color, Material material )
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
	}
};

// Used to describe a spherical surface:
class Sphere {
public:
	vec3 centre;
	float radius;
	vec3 color;
	Material material;

	Sphere( vec3 centre, float radius, vec3 color, Material material )
		: centre(centre), radius(radius), color(color), material(material)
	{

	}
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles, std::vector<Sphere>& spheres, std::vector<LightSource>& lights) {
	// Defines colors:
	vec3 red(    0.75f, 0.15f, 0.15f );
	vec3 yellow( 0.75f, 0.75f, 0.15f );
	vec3 green(  0.15f, 0.75f, 0.15f );
	vec3 cyan(   0.15f, 0.75f, 0.75f );
	vec3 blue(   0.15f, 0.15f, 0.75f );
	vec3 purple( 0.75f, 0.15f, 0.75f );
	vec3 white(  0.75f, 0.75f, 0.75f );

	float matteDiffuseRef = 0.05f;

	// Define materials
	Material matteWhite( white * matteDiffuseRef, vec3(0,0,0), 0.0f);
	Material matteRed( red * matteDiffuseRef, vec3(0,0,0), 0.0f);
	Material matteBlue( blue * matteDiffuseRef, vec3(0,0,0), 0.0f);
	Material matteGreen( green * matteDiffuseRef, vec3(0,0,0), 0.0f);
	Material matteYellow( yellow * matteDiffuseRef, vec3(0,0,0), 0.0f);
	Material mattePurple( purple * matteDiffuseRef, vec3(0,0,0), 0.0f);
	Material matteCyan( cyan * matteDiffuseRef, vec3(0,0,0), 0.0f);

	Material mirror( vec3(0,0,0), vec3(1,1,1), 0.0f);
  Material glass( vec3(0,0,0), vec3(0,0,0), 1.5f);

	lights.clear();
	lights.reserve( 2 );

	triangles.clear();
	triangles.reserve( 5*2*3 );

	spheres.clear();
	spheres.reserve( 1 );

	// ---------------------------------------------------------------------------
	// Room

	float L = 555;			// Length of Cornell Box side.

	vec3 A(L,0,0);
	vec3 B(0,0,0);
	vec3 C(L,0,L);
	vec3 D(0,0,L);

	vec3 E(L,L,0);
	vec3 F(0,L,0);
	vec3 G(L,L,L);
	vec3 H(0,L,L);

	//lights.push_back( LightSource( 40, vec3(1, 1, 1), vec3(0, -1.0, -0.5), vec3(0, 1, 0), 0.4, 0.4) );
  float middleX = 0.0;
  float middleZ = -0.5;
  float width = 0.25;
  float length = 0.25;
  float seperation = 0.03;
  vec3 direction(0, 1, 0);
  float lightPower = 60;

  float offsetX = (width / 2 + seperation);
  float offsetZ = (length / 2 + seperation);

  lights.push_back( LightSource( lightPower/4, vec3(1, 1, 1), vec3(middleX + offsetX, -1.0, middleZ + offsetZ), direction, width, length) );
  lights.push_back( LightSource( lightPower/4, vec3(1, 1, 1), vec3(middleX + offsetX, -1.0, middleZ - offsetZ), direction, width, length) );
  lights.push_back( LightSource( lightPower/4, vec3(1, 1, 1), vec3(middleX - offsetX, -1.0, middleZ + offsetZ), direction, width, length) );
  lights.push_back( LightSource( lightPower/4, vec3(1, 1, 1), vec3(middleX - offsetX, -1.0, middleZ - offsetZ), direction, width, length) );

	spheres.push_back( Sphere(vec3(0.5,0.75,-0.3), 0.2, white, glass) );

	// Floor:
	triangles.push_back( Triangle( C, B, A, white, matteWhite ) );
	triangles.push_back( Triangle( C, D, B, white, matteWhite ) );

	// Left wall
	triangles.push_back( Triangle( A, E, C, green, matteGreen ) );
	triangles.push_back( Triangle( C, E, G, green, matteGreen ) );

	// Right wall
	triangles.push_back( Triangle( F, B, D, red, matteRed ) );
	triangles.push_back( Triangle( H, F, D, red, matteRed ) );

	// Ceiling
	triangles.push_back( Triangle( E, F, G, white, matteWhite ) );
	triangles.push_back( Triangle( F, H, G, white, matteWhite ) );

	// Back wall
	triangles.push_back( Triangle( G, D, C, white, matteWhite ) );
	triangles.push_back( Triangle( G, H, D, white, matteWhite ) );

	// ---------------------------------------------------------------------------
	// Short block

	A = vec3(290,0,114);
	B = vec3(130,0, 65);
	C = vec3(240,0,272);
	D = vec3( 82,0,225);

	E = vec3(290,165,114);
	F = vec3(130,165, 65);
	G = vec3(240,165,272);
	H = vec3( 82,165,225);

	// Front
	// triangles.push_back( Triangle(E,B,A,white,matteWhite) );
	// triangles.push_back( Triangle(E,F,B,white,matteWhite) );
  //
	// // Front
	// triangles.push_back( Triangle(F,D,B,white,matteWhite) );
	// triangles.push_back( Triangle(F,H,D,white,matteWhite) );
  //
	// // BACK
	// triangles.push_back( Triangle(H,C,D,white,matteWhite) );
	// triangles.push_back( Triangle(H,G,C,white,matteWhite) );
  //
	// // LEFT
	// triangles.push_back( Triangle(G,E,C,white,matteWhite) );
	// triangles.push_back( Triangle(E,A,C,white,matteWhite) );
  //
	// // TOP
	// triangles.push_back( Triangle(G,F,E,white,matteWhite) );
	// triangles.push_back( Triangle(G,H,F,white,matteWhite) );

	// ---------------------------------------------------------------------------
	// Tall block

	A = vec3(423,0,247);
	B = vec3(265,0,296);
	C = vec3(472,0,406);
	D = vec3(314,0,456);

	E = vec3(423,330,247);
	F = vec3(265,330,296);
	G = vec3(472,330,406);
	H = vec3(314,330,456);

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

		triangles[i].v0 -= vec3(1,1,1);
		triangles[i].v1 -= vec3(1,1,1);
		triangles[i].v2 -= vec3(1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

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

  vec3 A(L,0,0);
	vec3 B(0,0,0);
	vec3 C(L,0,L);
	vec3 D(0,0,L);

	vec3 E(L,L,0);
	vec3 F(0,L,0);
	vec3 G(L,L,L);
	vec3 H(0,L,L);

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

  A = vec3(290,0,114);
	B = vec3(130,0, 65);
	C = vec3(240,0,272);
	D = vec3( 82,0,225);

	E = vec3(290,165,114);
	F = vec3(130,165, 65);
	G = vec3(240,165,272);
	H = vec3( 82,165,225);

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

  A = vec3(423,0,247);
	B = vec3(265,0,296);
	C = vec3(472,0,406);
	D = vec3(314,0,456);

	E = vec3(423,330,247);
	F = vec3(265,330,296);
	G = vec3(472,330,406);
	H = vec3(314,330,456);

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

	spheres.push_back( PhongSphere(vec3(0.4,0,-0.2), 0.3, shinyPurple) );
	spheres.push_back( PhongSphere(vec3(-0.4,0.7,-0.8), 0.2, shinyYellow) );

	// ---------------------------------------------------------------------------
	// Lights
	// ---------------------------------------------------------------------------

	//lights.push_back( Light(vec4(0.5,-0.5,-0.7,1.0), vec3(1,1,1), 0.02f, 0.3f, 1.0f ));
	//lights.push_back( Light(vec4(-0.5,-0.5,-0.7,1.0), vec3(1,1,1), 0.02f, 0.3f, 1.0f ));
	lights.push_back( PhongLightSource(vec3(0.0,-0.5,-0.7), vec3(1.0,1.0,1.0), 0.02f, 0.3f, 1.0f ));

	// ----------------------------------------------
	// Scale to the volume [-1,1]^3

	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0 *= 2/L;
		triangles[i].v1 *= 2/L;
		triangles[i].v2 *= 2/L;

		triangles[i].v0 -= vec3(1,1,1);
		triangles[i].v1 -= vec3(1,1,1);
		triangles[i].v2 -= vec3(1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].ComputeNormal();
	}
}


#endif
