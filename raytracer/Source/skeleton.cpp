#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include "PhongModel.h"
#include <stdint.h>
#include "limits"
#include <math.h>
#include "glm/gtx/string_cast.hpp"
#include <tuple>

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

#define PHOTON_MAPPER false

#define NUM_PHOTONS 10000

struct Intersection {
  vec4 position;
  float distance;
  int triangleIndex;
};

struct Photon {
  vec4 position;
  vec3 power;
  vec4 direction;
};

const float focalLength = SCREEN_HEIGHT;
const vec4 defaultCameraPos(0.0, 0.0, -3.0, 1.0);
vec4 cameraPos(0.0, 0.0, -3.0, 1.0);

const vec4 defaultLightPos( 0, -0.5, -0.7, 1.0 );
vec4 lightPos( 0, -0.5, -0.7, 1.0 );
vec3 lightColor = 14.f * vec3( 1, 1, 1 );
vec3 indirectLight = 0.5f * vec3( 1, 1, 1 );
const float shadowBiasThreshold = 0.001f;

std::vector<Triangle> triangles;
std::vector<LightSource> lights;

std::vector<TrianglePhong> phongTriangles;
std::vector<LightSourcePhong> phongLights;

std::vector<Photon> photonMap;

mat4 R;
float yaw = 0;
float pitch = 0;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, Intersection& closestIntersection);
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void updateRotation();
vec3 DirectLight( const Intersection& i );
void moveCameraRight(int direction);
void moveCameraUp(int direction);
void moveCameraForward(int direction);
void lookAt(mat4& ctw);
void emitPhotons();
void emitPhotonsFromLight(LightSource &light, int numPhotons);
bool tracePhoton(vec3 power, vec4 start, vec4 direction);
vec3 calculateRadiance(Intersection& intersection);
void drawPhotons(screen* screen);
vec3 getClosestPhotonPower(Intersection& intersection);
vec3 getNearestPhotonsPower(Intersection& intersection, int numNearest, float maxRadius);
void getNearestPhotonsIndex(Intersection& intersection, int numNearest, vector<int>& indices);
float getDist(vec4 a, vec4 b);
vec3 computeLight(const Intersection &i, const LightSourcePhong &l);

int main(int argc, char* argv[]) {
  screen *mainscreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  screen *photonscreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  if (PHOTON_MAPPER) LoadTestModel(triangles, lights);
  else LoadTestModelPhong(phongTriangles, phongLights);

  // Populate the photon map
  if (PHOTON_MAPPER) emitPhotons();

  while (Update()) {
    Draw(mainscreen);
    SDL_Renderframe(mainscreen);
    if (PHOTON_MAPPER) {
      drawPhotons(photonscreen);
      SDL_Renderframe(photonscreen);
    }
  }

  SDL_SaveImage(mainscreen, "mainout.bmp");
  if (PHOTON_MAPPER) SDL_SaveImage(photonscreen, "photonmapout.bmp");

  KillSDL(mainscreen);
  if (PHOTON_MAPPER) KillSDL(photonscreen);

	return 0;
}

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  float maxPixelVal = 0;
  vector<vec3> pixelVals;

  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
      vec4 d = normalize(R * vec4(x - SCREEN_WIDTH/2, y - SCREEN_HEIGHT/2, focalLength, 1));
      Intersection intersection;

      vec3 directLight(0.0, 0.0, 0.0); // The direct colour
			vec3 reflectedLight(0.0, 0.0, 0.0); // The visible colour
			vec3 colour(0.0, 0.0, 0.0); // The original colour of the triangle

      if (PHOTON_MAPPER) {
        if (ClosestIntersection(cameraPos, d, intersection)) {
          reflectedLight = calculateRadiance(intersection);
        }

        if (reflectedLight.x > maxPixelVal) maxPixelVal = reflectedLight.x;
        if (reflectedLight.y > maxPixelVal) maxPixelVal = reflectedLight.y;
        if (reflectedLight.z > maxPixelVal) maxPixelVal = reflectedLight.z;

        pixelVals.push_back(reflectedLight);
      } else {
        if (ClosestIntersection(cameraPos, d, intersection)) {
          colour = phongTriangles[intersection.triangleIndex].material.color;

          for (std::vector<LightSourcePhong>::size_type i = 0; i < phongLights.size(); i++) {
            directLight += computeLight(intersection, phongLights[i]);
          }

  				directLight = colour * (directLight + indirectLight);
        }
          PutPixelSDL(screen, x, y, directLight);
      }
    }
  }

  if (PHOTON_MAPPER) {
    for (int x = 0; x < SCREEN_WIDTH; x++) {
      for (int y = 0; y < SCREEN_HEIGHT; y++) {
        PutPixelSDL(screen, x, y, pixelVals[x*SCREEN_HEIGHT + y]/maxPixelVal);
      }
    }
  }
}

vec3 computeLight( const Intersection &i, const LightSourcePhong &l ) {
  // Using equation from https://en.wikipedia.org/wiki/Phong_reflection_model
  vec3 light(0.0,0.0,0.0);

  vec4 Lm = normalize(l.position - i.position);
  vec4 V = normalize(cameraPos - i.position);

  float dist = distance(l.position, i.position);

  float id = l.diffuseIntensity;
  float is = l.specularIntensity;
  float ia = l.ambientIntensity;

  float kd = 1;
  float ks = 1;
  float ka = 1;
  float alpha;

  vec4 normal;
  vec4 Rm;

  float LmNormalDot;
  float RmVDot;

  TrianglePhong triangle = phongTriangles[i.triangleIndex];

  ka = triangle.material.ambientRef;

  light += ka * (ia * l.color);

  kd = triangle.material.diffuseRef;

  normal = normalize(triangle.normal);

  LmNormalDot = dot(Lm, normal);

  ks = triangle.material.specularRef;

  Rm = normalize(2 * max(0.0f, LmNormalDot) * normal - Lm);

  alpha = triangle.material.shininess;

  RmVDot = dot(Rm, V);

  if (LmNormalDot > 0) {
    light += (kd * LmNormalDot * (id * l.color)) / dist;

    if (RmVDot > 0) {
      light += (ks * pow(RmVDot, alpha) * (is * l.color)) / dist;
    }
  }

  Intersection intersection;

  vec4 lightDir = l.position - i.position;
  vec4 rHat = normalize(lightDir);
  float r = distance(i.position, l.position);

  if (ClosestIntersection(i.position, rHat, intersection)) {
    if (intersection.triangleIndex != i.triangleIndex && intersection.distance < r) {
      light = ka * (ia * l.color);
    }
  }

  return light;
}

void drawPhotons(screen* screen) {
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  for (int i = 0; i < photonMap.size(); i++) {
    Photon p = photonMap[i];

    vec4 pos = R * p.position - cameraPos;

    int x = (focalLength * pos.x) / pos.z + SCREEN_WIDTH/2;
    int y = (focalLength * pos.y) / pos.z + SCREEN_WIDTH/2;
    if (x > 0 && x < SCREEN_WIDTH && y > 0 && y < SCREEN_HEIGHT) PutPixelSDL(screen, x, y, p.power);
  }

}

vec3 getClosestPhotonPower(Intersection& intersection) {
  float dist = numeric_limits<float>::max();
  vec3 accumPower;

  for (int i = 0; i < photonMap.size(); i++) {
    vec4 it = intersection.position;
    vec4 phop = photonMap[i].position;

    float distVec = getDist(it, phop);

    if (distVec < dist) {
      accumPower = photonMap[i].power;
      dist = distVec;
    }
  }

  return accumPower;
}

float getDist(vec4 a, vec4 b) {
  return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

void getNearestPhotonsIndex(Intersection& intersection, int numNearest, vector<int>& indices) {
  vector<float> distances;

  indices.reserve(numNearest);
  indices.clear();

  for (int i = 0; i < photonMap.size(); i++) {
    vec4 ipos = intersection.position;
    vec4 ppos = photonMap[i].position;

    float dist = getDist(ipos, ppos);

    if (indices.size() < numNearest) {
      indices.push_back(i);
      distances.push_back(dist);
    } else {
      int furthestIndex = 0;
      float maxDist = -numeric_limits<float>::max();

      for (int j = 0; j < numNearest; j++) {
        if (distances[j] > maxDist) {
          maxDist = distances[j];
          furthestIndex = j;
        }
      }

      if (dist < maxDist) {
        indices[furthestIndex] = i;
        distances[furthestIndex] = dist;
      }
    }
  }
}

vec3 getNearestPhotonsPower(Intersection& intersection, int numNearest, float maxRadius) {
  vector<int> nearestPhotonsIndex;
  vec3 accumPower = vec3(0,0,0);
  float radius = 0;

  getNearestPhotonsIndex(intersection, numNearest, nearestPhotonsIndex);

  for (int i = 0; i < numNearest; i++) {
    accumPower += photonMap[nearestPhotonsIndex[i]].power;
    float dist = getDist(intersection.position, photonMap[nearestPhotonsIndex[i]].position);
    if (dist > radius) {
      radius = dist;
    }
  }

  vec3 unitPower = accumPower / (float) ((4/3) * M_PI * pow(radius, 3));

  float maxPower = NUM_PHOTONS / 3;

  if (unitPower.x > maxPower) unitPower = unitPower * (maxPower / unitPower.x);
  if (unitPower.y > maxPower) unitPower = unitPower * (maxPower / unitPower.y);
  if (unitPower.z > maxPower) unitPower = unitPower * (maxPower / unitPower.z);

  return unitPower;
}

vec3 calculateRadiance(Intersection& intersection) {
  return getNearestPhotonsPower(intersection, 100, 0.05);
}

void emitPhotons() {
  photonMap.reserve( NUM_PHOTONS );

  int numLights = lights.size();

  cout << "-------------------" << endl;
  cout << "emitting photons from lights" << endl;

  for (int i = 0; i < numLights; i++) {
    emitPhotonsFromLight(lights[i], NUM_PHOTONS / numLights);
  }

  cout << "finished emitting photons from lights" << endl;
  cout << "-------------------" << endl;
}

void emitPhotonsFromLight(LightSource &light, int numPhotons) {
  int numIntersections = 0;

  for(int i = 0; i < numPhotons; i++) {
    vec4 direction = vec4(1,1,1,1);

    while ((pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2)) > 1) {
      direction.x = ((float) rand() / (RAND_MAX)) * 2 - 1;
      direction.y = ((float) rand() / (RAND_MAX)) * 2 - 1;
      direction.z = ((float) rand() / (RAND_MAX)) * 2 - 1;
    }

    direction.w = 1;

    if (tracePhoton((light.watts / numPhotons) * light.color, light.position, normalize(direction))) numIntersections++;
  }

  cout << "number of photon intersections: " << numIntersections << endl;
}

bool tracePhoton(vec3 power, vec4 start, vec4 direction) {
  Intersection intersection;

  if (ClosestIntersection(start, direction, intersection)) {
    Triangle triangle = triangles[intersection.triangleIndex];

    Material material = triangle.material;

    vec3 diffuseRef = material.diffuseRef;
    vec3 specRef = material.specRef;

    // Probability of reflection
    float Pr = max(diffuseRef.x + specRef.x, max(diffuseRef.y + specRef.y, diffuseRef.z + specRef.z));

    // Probability of diffuse reflection
    float Pd = Pr * ((diffuseRef.x + diffuseRef.y + diffuseRef.z) /
                      (diffuseRef.x + diffuseRef.y + diffuseRef.z + specRef.x + specRef.y + specRef.z));

    // Probability of specular reflection
    float Ps = Pr - Pd;

    float rnd = ((float) rand() / (RAND_MAX));

    if (rnd < Pd) {
      //Diffuse reflection
      vec4 reflectionDir = normalize(2 * dot(triangle.normal, direction) * triangle.normal - direction);
      tracePhoton(power, intersection.position, reflectionDir);

    } else if (rnd < Ps + Pd) {
      // Specular reflection
      vec4 reflectionDir = normalize(2 * dot(triangle.normal, direction) * triangle.normal - direction);

      vec3 specPower = vec3(power.x * specRef.x / Ps, power.y * specRef.y / Ps, power.z * specRef.z / Ps);
      tracePhoton(specPower, intersection.position, reflectionDir);

    } else {
      // Absorbtion
      Photon p;
      p.position = intersection.position;
      p.direction = direction;
      p.power = triangle.color;
      photonMap.push_back(p);
    }

    return true;
  }

  return false;
}

bool ClosestIntersection(vec4 start, vec4 dir, Intersection& closestIntersection) {
  bool intersectionFound = false;
  closestIntersection = Intersection();
  closestIntersection.distance = std::numeric_limits<float>::max();

  float size = (PHOTON_MAPPER) ? triangles.size() : phongTriangles.size();

  for(std::vector<Triangle>::size_type i = 0; i < size; i++) {
    vec4 v0 = (PHOTON_MAPPER) ? triangles[i].v0 : phongTriangles[i].v0;
    vec4 v1 = (PHOTON_MAPPER) ? triangles[i].v1 : phongTriangles[i].v1;
    vec4 v2 = (PHOTON_MAPPER) ? triangles[i].v2 : phongTriangles[i].v2;

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

	if (ClosestIntersection(i.position, rHat, intersection)) {
		if (intersection.triangleIndex != i.triangleIndex && intersection.distance < r) {
			C = black;
		}
	}

	return C;
}

void lookAt(mat4& ctw) { /* TODO! */
	vec3 toPos = vec3(0, 0, 0);
	// vec3 fromPos = vec3(cameraPos);
	vec3 fromPos = vec3(cameraPos.x, cameraPos.y, cameraPos.z);

	vec3 forward = -normalize(fromPos - toPos);

	vec3 tmp = vec3(0, 1, 0);
	vec3 right = cross(normalize(tmp), forward);

	vec3 up = cross(forward, right);

	ctw[0][0] = right.x;
	ctw[0][1] = right.y;
	ctw[0][2] = right.z;
	ctw[1][0] = up.x;
	ctw[1][1] = up.y;
	ctw[1][2] = up.z;
	ctw[2][0] = forward.x;
	ctw[2][1] = forward.y;
	ctw[2][2] = forward.z;

	ctw[3][0] = fromPos.x;
	ctw[3][1] = fromPos.y;
	ctw[3][2] = fromPos.z;

	ctw[3][3] = 1;

  //ctw = transpose(ctw);
	R = ctw;

  pitch = asin(-forward.y);
  yaw = atan2(forward.x, forward.z);

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

	mat4 ctw;

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
					lookAt(ctw);
					break;
        case SDLK_t:
          // Reset camera position
          cameraPos = defaultCameraPos;
          lightPos = defaultLightPos;
          pitch = 0;
          yaw = 0;
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
