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
struct Light {
  glm::vec4 position;
  glm::vec3 color;
	float ambientIntensity;
	float diffuseIntensity;
	float specularIntensity;
  Light(const vec4 &p, const vec3 &c, const float &am, const float &di, const float &sp)
		: position(p), color(c), ambientIntensity(am), diffuseIntensity(di), specularIntensity(sp) {

	}
};

struct Material {
	vec3 color;
	float ambientRef;
	float diffuseRef;
	float specularRef;
	float shininess;
	Material(const vec3 &c, const float &am, const float &di, const float &sp, const float &sh)
		: color(c), ambientRef(am), diffuseRef(di), specularRef(sp), shininess(sh) {

	}
};

// Used to describe a triangular surface:
class Triangle {
public:
	glm::vec4 v0;
	glm::vec4 v1;
	glm::vec4 v2;
	glm::vec4 normal;
	Material material;

	Triangle( glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, Material material )
		: v0(v0), v1(v1), v2(v2), material(material)
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
	Material material;

	Sphere( glm::vec4 centre, float radius, Material material )
		: centre(centre), radius(radius), material(material)
	{

	}
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles, std::vector<Sphere>& spheres, std::vector<Light>& lights )
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

	Material matteRed(red, 1, 2, 1, 1);
	Material matteYellow(yellow, 1, 2, 1, 1);
	Material matteGreen(green, 1, 2, 1, 1);
	Material matteCyan(cyan, 1, 2, 1, 1);
	Material matteBlue(blue, 1, 2, 1, 1);
	Material mattePurple(purple, 1, 2, 1, 1);
	Material matteWhite(white, 1, 2, 2, 1);

	Material shinyPurple(darkPurple, 1.2, 3, 4, 7);

	triangles.clear();
	triangles.reserve( 5*2*3 );

	spheres.clear();
	spheres.reserve(1);

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
	triangles.push_back( Triangle( C, B, A, matteGreen ) );
	triangles.push_back( Triangle( C, D, B, matteGreen ) );

	// Left wall
	triangles.push_back( Triangle( A, E, C, mattePurple ) );
	triangles.push_back( Triangle( C, E, G, mattePurple ) );

	// Right wall
	triangles.push_back( Triangle( F, B, D, matteYellow ) );
	triangles.push_back( Triangle( H, F, D, matteYellow ) );

	// Ceiling
	triangles.push_back( Triangle( E, F, G, matteCyan ) );
	triangles.push_back( Triangle( F, H, G, matteCyan ) );

	// Back wall
	triangles.push_back( Triangle( G, D, C, matteWhite ) );
	triangles.push_back( Triangle( G, H, D, matteWhite ) );

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
	triangles.push_back( Triangle(E,B,A,matteRed) );
	triangles.push_back( Triangle(E,F,B,matteRed) );

	// Front
	triangles.push_back( Triangle(F,D,B,matteRed) );
	triangles.push_back( Triangle(F,H,D,matteRed) );

	// BACK
	triangles.push_back( Triangle(H,C,D,matteRed) );
	triangles.push_back( Triangle(H,G,C,matteRed) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,matteRed) );
	triangles.push_back( Triangle(E,A,C,matteRed) );

	// TOP
	triangles.push_back( Triangle(G,F,E,matteRed) );
	triangles.push_back( Triangle(G,H,F,matteRed) );

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
	triangles.push_back( Triangle(E,B,A,matteBlue) );
	triangles.push_back( Triangle(E,F,B,matteBlue) );

	// Front
	triangles.push_back( Triangle(F,D,B,matteBlue) );
	triangles.push_back( Triangle(F,H,D,matteBlue) );

	// BACK
	triangles.push_back( Triangle(H,C,D,matteBlue) );
	triangles.push_back( Triangle(H,G,C,matteBlue) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,matteBlue) );
	triangles.push_back( Triangle(E,A,C,matteBlue) );

	// TOP
	triangles.push_back( Triangle(G,F,E,matteBlue) );
	triangles.push_back( Triangle(G,H,F,matteBlue) );

	// ---------------------------------------------------------------------------
	// Spheres
	// ---------------------------------------------------------------------------

	spheres.push_back( Sphere(vec4(0.4,0,-0.2,1), 0.3, shinyPurple) );

	// ---------------------------------------------------------------------------
	// Lights
	// ---------------------------------------------------------------------------

	//lights.push_back( Light(vec4(0.5,-0.5,-0.7,1.0), vec3(1,1,1), 0.02f, 0.3f, 1.0f ));
	//lights.push_back( Light(vec4(-0.5,-0.5,-0.7,1.0), vec3(1,1,1), 0.02f, 0.3f, 1.0f ));
	lights.push_back( Light(vec4(0.0,-0.5,-0.9,1.0), vec3(1,1,1), 0.02f, 0.3f, 1.0f ));

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
