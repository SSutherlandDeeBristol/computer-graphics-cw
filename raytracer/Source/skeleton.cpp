#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>
#include "limits"
#include <math.h>
#include "mpi.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::distance;

SDL_Event event;

#define SCREEN_WIDTH 300
#define SCREEN_HEIGHT 300
#define FULLSCREEN_MODE false

struct Intersection {
  vec4 position;
  float distance;
  int triangleIndex;
};

float focalLength = SCREEN_HEIGHT;
vec4 cameraPos(0.0, 0.0, -3.0, 1.0);

vec4 lightPos( 0, -0.5, -0.7, 1.0 );
vec3 lightColor = 14.f * vec3( 1, 1, 1 );
vec3 indirectLight = 0.5f * vec3( 1, 1, 1 );

std::vector<Triangle> triangles;
mat4 R;

float yaw = 0;
float pitch = 0;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void updateRotation();
vec3 DirectLight( const Intersection& i );
void moveCameraRight(int direction);
void moveCameraUp(int direction);
void moveCameraForward(int direction);

int main(int argc, char* argv[]) {
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  LoadTestModel(triangles);

  while (Update()) {
    Draw(screen);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);

	return 0;
}

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
      vec4 d = normalize(R * vec4(x - SCREEN_WIDTH/2, y - SCREEN_HEIGHT/2, focalLength, 1));
      Intersection intersection;

      vec3 directLight(0.0, 0.0, 0.0); // The direct colour
			vec3 reflectedLight(0.0, 0.0, 0.0); // The visible colour
			vec3 colour(0.0, 0.0, 0.0); // The original colour of the triangle

      if (ClosestIntersection(cameraPos, d, triangles, intersection)) {
				directLight = DirectLight(intersection);
				colour = triangles[intersection.triangleIndex].color;

				reflectedLight = colour * (directLight + indirectLight);
      }

      PutPixelSDL(screen, x, y, reflectedLight);
    }
  }
}

bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection) {
  bool intersectionFound = false;
  closestIntersection = Intersection();
  closestIntersection.distance = std::numeric_limits<float>::max();

  for(std::vector<Triangle>::size_type i = 0; i < triangles.size(); i++) {
    vec4 v0 = triangles[i].v0;
    vec4 v1 = triangles[i].v1;
    vec4 v2 = triangles[i].v2;

    vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
    vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
    vec3 b = vec3(start.x-v0.x,start.y-v0.y,start.z-v0.z);
    vec3 d = vec3(dir.x, dir.y, dir.z);

    mat3 A = mat3(-d, e1, e2);

    mat3 A1(b, e1, e2);
    float detA = glm::determinant(A);

    float t = glm::determinant(A1) / detA;

    if (t > 0) {
			mat3 A2(-vec3(dir), b, e2);
			mat3 A3(-vec3(dir), e1, b);

			float u = glm::determinant(A2) / detA;
			float v = glm::determinant(A3) / detA;

			if (u >= 0 && v >= 0 && (u + v) <= 1) {
				// Intersection occured
				vec4 position = start + t * dir;
				float dist = distance(start, position);

				if (dist <= closestIntersection.distance) {
					intersectionFound = true;
					closestIntersection.distance = dist;
					closestIntersection.position = position;
					closestIntersection.triangleIndex = i;
				}
			}
    }
  }

  return intersectionFound;
}

void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R) {
	R[0][0] = cos(thetaY) * cos(thetaZ);
	R[0][1] = -cos(thetaX) * sin(thetaZ) + sin(thetaX) * sin(thetaY) * cos(thetaZ);
	R[0][2] = sin(thetaX) * sin(thetaZ) + cos(thetaX) * sin(thetaY) * cos(thetaZ);

	R[1][0] = cos(thetaY) * sin(thetaZ);
	R[1][1] = cos(thetaX) * cos(thetaZ) + sin(thetaX) * sin(thetaY) * sin(thetaZ);
	R[1][2] = -sin(thetaX) * cos(thetaZ) + cos(thetaX) * sin(thetaY) * sin(thetaZ);

	R[2][0] = -sin(thetaY);
	R[2][1] = sin(thetaX) * cos(thetaY);
	R[2][2] = cos(thetaX) * cos(thetaY);
}

void updateRotation() {
	mat3 RT;

	getRotationMatrix(pitch, yaw, 0, RT);

	R = mat4(RT);
	
	R = transpose(R);
}

void moveCameraRight(int direction, float distance) {
	vec4 right(R[0][0], R[0][1], R[0][2], 0);
	
	cameraPos += direction * distance * right;
}

void moveCameraUp(int direction, float distance) {
	vec4 up(R[1][0], R[1][1], R[1][2], 0);
	
	cameraPos += direction * distance * up;
}

void moveCameraForward(int direction, float distance) {
	vec4 forward(R[2][0], R[2][1], R[2][2], 0);
	
	cameraPos += direction * distance * forward;
}

vec3 DirectLight( const Intersection& i ) {
	Triangle triangle = triangles[i.triangleIndex];

	float r = distance(i.position, lightPos);

	float A = 4 * M_PI * pow(r, 2);

	vec4 lightDir = lightPos - i.position;

	vec4 rHat = normalize(lightDir);
	vec4 nHat = normalize(triangle.normal);

	vec3 B = lightColor / A;
	vec3 D = B * max(dot(rHat,nHat), 0.0f);
	vec3 C = D * triangle.color;

	Intersection intersection;
	vec3 black = vec3(0.0, 0.0, 0.0); // Initialise to black

	if (ClosestIntersection(i.position, rHat, triangles, intersection)) {
		if (intersection.triangleIndex != i.triangleIndex && intersection.distance < r) {
			C = black;
		}
	}

	return C;
}

/*Place updates of parameters here*/
bool Update() {
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  std::cout << "Render time: " << dt << " ms." << std::endl;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
	    return false;
	  } else if (e.type == SDL_KEYDOWN) {
	    int key_code = e.key.keysym.sym;
	    switch(key_code) {
	      case SDLK_UP:
					pitch -= M_PI / 18;
					updateRotation();
					break;
	      case SDLK_DOWN:
					pitch += M_PI / 18;
					updateRotation();
          break;
	      case SDLK_LEFT:
					yaw -= M_PI / 18;
					updateRotation();
          break;
	      case SDLK_RIGHT:
					yaw += M_PI / 18;
					updateRotation();
          break;
				case SDLK_w:
					moveCameraUp(-1, 0.25);
					break;
				case SDLK_s:
					moveCameraUp(1, 0.25);
					break;
				case SDLK_a:
					moveCameraRight(-1, 0.25);
					break;
				case SDLK_d:
					moveCameraRight(1, 0.25);
					break;
				case SDLK_EQUALS:
					moveCameraForward(1, 0.25);
					break;
				case SDLK_MINUS:
					moveCameraForward(-1, 0.25);
					break;
				case SDLK_i:
					lightPos.z += 0.5;
					break;
				case SDLK_k:
					lightPos.z -= 0.5;
					break;
				case SDLK_j:
					lightPos.x -= 0.5;
					break;
				case SDLK_l:
					lightPos.x += 0.5;
					break;
				case SDLK_o:
					lightPos.y -= 0.5;
					break;
				case SDLK_p:
					lightPos.y += 0.5;
					break;
	      case SDLK_ESCAPE:
          /* Move camera quit */
          return false;
      }
    }
  }
  return true;
}
