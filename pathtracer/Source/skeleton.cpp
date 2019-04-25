#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>
#include "limits"
#include <math.h>
#include <random>
#include "glm/gtx/string_cast.hpp" // std::cout<<glm::to_string(hello)<<std::endl;

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

const int PATH_TRACER_BOUNCES = 2;
const int PATH_TRACER_SAMPLES = 16;
std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);

const float focalLength = SCREEN_HEIGHT;
const float shadowBiasThreshold = 0.001f;
const vec4 defaultCameraPos(0.0, 0.0, -3.0, 1.0);
const vec4 defaultLightPos(0.0, -0.5, -0.7, 1.0);

vec4 cameraPos(0.0, 0.0, -3.0, 1.0);
vec4 lightPos(0.0, -0.5, -0.7, 1.0);
vec3 lightColor = 14.f * vec3( 1, 1, 1 );
// vec3 indirectLight = 0.5f * vec3( 1, 1, 1 );
vec3 indirectLight = 0.1f * vec3( 1, 1, 1 );

std::vector<Triangle> triangles;
mat4 R;

float yaw = 0;
float pitch = 0;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
vec3 CastRay(vec4 start, vec4 dir, const vector<Triangle>& triangles, int depth);
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void updateRotation();
vec3 DirectLight( const Intersection& i );
void moveCameraRight(int direction);
void moveCameraUp(int direction);
void moveCameraForward(int direction);
void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb);
vec3 uniformSampleHemisphere(const float &r1, const float &r2);
void lookAt(vec3 toPos);
void resetView();

int main(int argc, char* argv[]) {
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  resetView();
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
      // Intersection intersection;

      // vec3 directLight(0.0, 0.0, 0.0); // The direct colour
			// vec3 reflectedLight(0.0, 0.0, 0.0); // The visible colour
			// vec3 colour(0.0, 0.0, 0.0); // The original colour of the triangle

      // if (ClosestIntersection(cameraPos, d, triangles, intersection)) {
			// 	directLight = DirectLight(intersection);
			// 	colour = triangles[intersection.triangleIndex].color;

			// 	reflectedLight = colour * (directLight + indirectLight);
      // }

      vec3 colour = CastRay(cameraPos, d, triangles, 0);

      // PutPixelSDL(screen, x, y, reflectedLight);
      PutPixelSDL(screen, x, y, colour);
    }
  }
}

vec3 CastRay(vec4 start, vec4 dir, const vector<Triangle>& triangles, int depth) {

  if (depth > PATH_TRACER_BOUNCES) return vec3(0, 0, 0);

  Intersection intersection;

  if (ClosestIntersection(start, dir, triangles, intersection)) {
    vec3 directLight(0.0, 0.0, 0.0); // The direct colour
    vec3 reflectedLight(0.0, 0.0, 0.0); // The visible colour
    vec3 colour(0.0, 0.0, 0.0); // The original colour of the triangle
    directLight = DirectLight(intersection);

    vec3 normal = vec3(triangles[intersection.triangleIndex].normal); // * glm::vec4(-1,-1,-1,1));

    vec3 Nt, Nb;

    createCoordinateSystem(normal, Nt, Nb);
    float pdf = 1 / (2 * M_PI);

    // number of samples N
    uint32_t N = PATH_TRACER_SAMPLES;
    vec3 indirectDiffuse = vec3(0, 0, 0);
    for (uint32_t i = 0; i < N; i++) {
      // step 2: create sample in world space
      float r1 = distribution(generator), r2 = distribution(generator);
      vec3 sample = uniformSampleHemisphere(r1, r2);

      // step 3: transform sample from world space to shaded point local coordinate system
      vec4 sampleWorld = vec4(
        sample.x * Nb.x + sample.y * normal.x + sample.z * Nt.x,
        sample.x * Nb.y + sample.y * normal.y + sample.z * Nt.y,
        sample.x * Nb.z + sample.y * normal.z + sample.z * Nt.z, 1);
      // step 4 & 5: cast a ray in this direction
      indirectDiffuse += r1 * CastRay(intersection.position,
        normalize(sampleWorld), triangles, depth + 1) / pdf;
    }
    // step 7: divide the sum by the total number of samples N
    indirectDiffuse /= (float) N;

    // std::cout<<glm::to_string(indirectDiffuse) << ", " << depth <<std::endl;

    colour = triangles[intersection.triangleIndex].color;
    reflectedLight = (directLight / (float) M_PI + 2.0f * indirectDiffuse + indirectLight) * colour;
    return reflectedLight;
  } else {
    return vec3(0, 0, 0);
  }
}

vec3 uniformSampleHemisphere(const float &r1, const float &r2) {
  // cos(theta) = u1 = y
  // cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
  float sinTheta = sqrtf(1 - r1 * r1);
  float phi = 2 * M_PI * r2;
  float x = sinTheta * cosf(phi);
  float z = sinTheta * sinf(phi);
  return vec3(x, r1, z);
}

/* http://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing/global-illumination-path-tracing-practical-implementation */
void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb) {
  if (std::fabs(N.x) > std::fabs(N.y)) {
    Nt = vec3(N.z, 0, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
  } else {
    Nt = vec3(0, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z);
  }
  Nb = cross(N, Nt);
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

				if (dist <= closestIntersection.distance && dist > shadowBiasThreshold) {
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
	R = transpose(mat4(RT));
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
	vec3 D = B * max(dot(rHat, nHat), 0.0f);
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

void lookAt(vec3 toPos) {
	vec3 fromPos = vec3(cameraPos);
	vec3 forward = -normalize(fromPos - toPos);

  pitch = asin(-forward.y);
  yaw = atan2(forward.x, forward.z);

  updateRotation();
}

void resetView() {
  cameraPos = defaultCameraPos;
  lightPos = defaultLightPos;
  pitch = 0;
  yaw = 0;
  updateRotation();
}

/*Place updates of parameters here*/
bool Update() {
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  std::cout << "Render time: " << dt << " ms." << std::endl;
	// std::cout << "cx: " << cameraPos.x << ", cy:" << cameraPos.y << ", cz:"<< cameraPos.z << std::endl;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
	    return false;
	  } else if (e.type == SDL_KEYDOWN) {
	    int key_code = e.key.keysym.sym;
	    switch(key_code) {
	      case SDLK_UP:
					pitch += M_PI / 18;
					updateRotation();
					break;
	      case SDLK_DOWN:
					pitch -= M_PI / 18;
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
				case SDLK_r:
					/* Look-At function, points camera to 0,0,0 */
					lookAt(vec3(0, 0, 0));
					break;
        case SDLK_t:
          // Reset camera position
          resetView();
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