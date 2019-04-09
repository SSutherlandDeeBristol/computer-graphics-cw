#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

class LightSource {
public:
	float watts;
	glm::vec3 color;
	glm::vec4 position;
	glm::vec4 direction;
	float width;
	float height;

	LightSource( float watts, glm::vec3 color, glm::vec4 position, glm::vec4 direction, float width, float height)
		: watts(watts), color(color), position(position), direction(direction), width(width), height(height)
	{

	}
};

class Material {
public:
	glm::vec3 diffuseRef;
	glm::vec3 specRef;

	Material( glm::vec3 diffuseRef, glm::vec3 specRef)
		: diffuseRef(diffuseRef), specRef(specRef)
	{

	}
};

// Used to describe a triangular surface:
class Triangle {
public:
	glm::vec4 v0;
	glm::vec4 v1;
	glm::vec4 v2;
	glm::vec4 normal;
	glm::vec3 color;
	Material material;

	Triangle( glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, glm::vec3 color, Material material )
		: v0(v0), v1(v1), v2(v2), color(color), material(material)
	{
		ComputeNormal();
	}

	void ComputeNormal()
	{
	  glm::vec3 e1 = glm::vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
	  glm::vec3 e2 = glm::vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
	  glm::vec3 normal3 = glm::normalize( glm::cross( e2, e1 ) );
	  normal.x = normal3.x;
	  normal.y = normal3.y;
	  normal.z = normal3.z;
	  normal.w = 1.0;
	}
};

// Used to describe a spherical surface:
class Sphere {
public:
	glm::vec4 centre;
	float radius;
	glm::vec3 color;

	Sphere( glm::vec4 centre, float radius, glm::vec3 color )
		: centre(centre), radius(radius), color(color)
	{

	}
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles, std::vector<LightSource>& lights)
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

	float matteDiffuseRef = 0.5f;

	// Define materials
	Material matteWhite( white * matteDiffuseRef, vec3(0,0,0));
	Material matteRed( red * matteDiffuseRef, vec3(0,0,0));
	Material matteBlue( blue * matteDiffuseRef, vec3(0,0,0));
	Material matteGreen( green * matteDiffuseRef, vec3(0,0,0));
	Material matteYellow( yellow * matteDiffuseRef, vec3(0,0,0));
	Material mattePurple( purple * matteDiffuseRef, vec3(0,0,0));
	Material matteCyan( cyan * matteDiffuseRef, vec3(0,0,0));

	lights.clear();
	lights.reserve( 1 );

	lights.push_back( LightSource( 2000, vec3(1, 1, 1), vec4( 0, -0.4, -0.5, 1.0 ), vec4(0, -1, 0, 1), 0.2, 0.2));

	triangles.clear();
	triangles.reserve( 5*2*3 );

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
