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

int PATH_TRACER_BOUNCES = 1;
int PATH_TRACER_SAMPLES = 16;
std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);

const float focalLength = SCREEN_HEIGHT * 3/2;
const float shadowBiasThreshold = 0.001f;
const float bounceBiasThreshold = 0.01f;
const vec4 defaultCameraPos(0.0, 0.0, -4.0, 1.0);
const vec4 defaultLightPos(0.0, -0.5, -0.25, 1.0);

vec4 cameraPos(0.0, 0.0, -3.0, 1.0);
vec4 lightPos(0.0, -0.5, -0.25, 1.0);
vec3 lightColor = 60.f * vec3(1, 1, 1);
vec3 distantEnvironmentLight = 0.3f * vec3(1, 1, 1);

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
void createCoordSystem(vec3 N, vec3 &nX, vec3 &nY);
void uniformSampleHemisphere(float r1, float r2, vec3 &sample);
void getSampleTransform(vec3 normal, vec3 nX, vec3 nY, mat3 &transform);
void lookAt(vec3 toPos);
void reflect(vec3 dir, vec3 normal, vec3& reflectionDir);
void resetView();

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " PATH_TRACER_BOUNCES PATH_TRACER_SAMPLES" << endl;
    return 1;
  } else {
    try {
      PATH_TRACER_BOUNCES = stoi(argv[1]);
      PATH_TRACER_SAMPLES = stoi(argv[2]);
    } catch (exception const &e) {
      cerr << "Could not parse arguments." << endl;
    }
  }

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  resetView();
  LoadTestModel(triangles);

  while (Update()) {
    Draw(screen);
    SDL_Renderframe(screen);
    SDL_SaveImage(screen, "mainout.bmp");
  }

  KillSDL(screen);

	return 0;
}

void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  float lastpercentage = 0.0f;
  float screenSize = (float) (SCREEN_WIDTH * SCREEN_HEIGHT);

  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
      /* Calculate direction of ray */
      vec4 d = normalize(R * vec4(x - SCREEN_WIDTH/2, y - SCREEN_HEIGHT/2, focalLength, 1));

      /* Cast a bouncing ray in this direction, retrieving the colour*/
      vec3 colour = CastRay(cameraPos, d, triangles, 0);

      PutPixelSDL(screen, x, y, colour);

      // print out the progress of the render
      float pixelNum = x * SCREEN_HEIGHT + y;
      float percentage = floor(100.0 * pixelNum / screenSize);

      if (percentage > lastpercentage) {
        cout << "\33[2k";
        cout << percentage << "% complete\r" << flush;
      }

      lastpercentage = percentage;
    }
  }
}

vec3 CastRay(vec4 start, vec4 dir, const vector<Triangle>& triangles, int depth) {
  vec3 black = vec3(0, 0, 0);
  /* Escape from recursion if depth exceeds limit */
  if (depth > PATH_TRACER_BOUNCES) return black;

  Intersection intersection;
  float pdf = 1 / (2 * M_PI);
  vec3 reflectedLight = black;

  /* Find the closest intersection of the ray */
  if (ClosestIntersection(start, dir, triangles, intersection)) {

    /* Calculate the direct light at this point */
    Triangle intersectedTriangle = triangles[intersection.triangleIndex];
    int firedRays = 0;
    vec3 directLight = DirectLight(intersection);
    vec3 normal = vec3(intersectedTriangle.normal);
    vec3 nX, nY;
    vec3 indirectLight = black;

    /* Create a coordinate system at this point, with the y-axis aligning with the normal */
    createCoordSystem(normal, nX, nY);

    for (int i = 0; i < PATH_TRACER_SAMPLES; i++) {
      if (intersectedTriangle.mirror) {
        vec3 reflectedDir; reflect(vec3(dir), normal, reflectedDir);
        indirectLight += CastRay(intersection.position - (dir * bounceBiasThreshold), normalize(vec4(reflectedDir, 1)), triangles, depth + 1) / pdf;
        firedRays++;
      }
      /* Generate two random values between 0 and 1 */
      float r1 = distribution(generator), r2 = distribution(generator);

      /* Create a sample vector from these random values */
      vec3 sample; uniformSampleHemisphere(r1, r2, sample);

      /* Transform the sample to world space */
      mat3 sampleTransform; getSampleTransform(normal, nX, nY, sampleTransform);
      vec4 sampleWorldSpace = vec4(sample * sampleTransform, 1);

      /* Fire a recursive ray in this sample direction, weighting the result */
      indirectLight += r1 * CastRay(intersection.position - (dir * bounceBiasThreshold), normalize(sampleWorldSpace), triangles, depth + 1) / pdf;
      firedRays++;
    }

    /* Average the sum of the samples */
    indirectLight /= (float) firedRays;

    /* Compute and return reflected light */
    vec3 colour = intersectedTriangle.color;
    reflectedLight = (directLight / (float) M_PI + 2.0f * indirectLight) * colour;
  } else reflectedLight = distantEnvironmentLight;

  return reflectedLight / (float) (pow((depth + 1), 2));
}

void reflect(vec3 dir, vec3 normal, vec3& reflectionDir) {
  reflectionDir = normalize(dir - 2 * dot(dir, normal) * normal);
}

void getSampleTransform(vec3 normal, vec3 nX, vec3 nY, mat3 &transform) {
  transform[0][0] = nY.x; transform[0][1] = normal.x; transform[0][2] = nX.x;
  transform[1][0] = nY.y; transform[1][1] = normal.y; transform[1][2] = nX.y;
  transform[2][0] = nY.z; transform[2][1] = normal.z; transform[2][2] = nX.z;
}

void uniformSampleHemisphere(float r1, float r2, vec3 &sample) {
  float phi = 2 * M_PI * r2;
  float sineTheta = sqrtf(1 - pow(r1, 2));
  sample = vec3(sineTheta * cosf(phi), r1, sineTheta * sinf(phi));
}

void createCoordSystem(vec3 N, vec3 &nX, vec3 &nY) {
  bool nXB = fabs(N.x) > fabs(N.y);
  nX = vec3(nXB ? N.z : 0, nXB ? 0 : -N.z, nXB ? -N.x : N.y) / sqrtf((nXB ? N.x * N.x : N.y * N.y) + N.z * N.z);
  nY = cross(N, nX);
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
	mat3 RT; getRotationMatrix(pitch, yaw, 0, RT);
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
