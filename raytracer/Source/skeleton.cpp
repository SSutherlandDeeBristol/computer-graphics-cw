#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
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

// Maximum number of times a photon can bounce
#define MAX_PHOTON_DEPTH 50

// Maximum number of times a ray from the camera can reflect/refract
#define MAX_CAMERA_RAY_DEPTH 50

// Filter constant when estimating radiance from the photon maps
#define FILTER_CONSTANT 1.005

// Number of rays to fire at each light source when
// calculating direct illumination
// Must be square number
#define NUM_SHADOW_RAYS 16

// Angle of the spotlight cutoff for the lightsources
#define SPOTLIGHT_CUTOFF M_PI/2

// Global refractive index
#define GLOBAL_REF_INDEX 1

enum geometry {triangle, sphere};
enum bounce {diffuse, specular, none};

struct Intersection {
  vec3 position;
  vec3 normal;
  vec3 direction;
  float distance;
  int index;
  geometry intersectionType;
};

struct Photon {
  vec3 position;
  vec3 power;
  vec3 direction;
};

// Variables reassigned from command line
bool PHOTON_MAPPER = false;
int NUM_PHOTONS = 5000;
int NUM_NEAREST_PHOTONS = 200;

// Camera values
const float focalLength = SCREEN_HEIGHT;
const vec4 defaultCameraPos(0.0, 0.0, -3.0, 1.0);
vec4 cameraPos(0.0, 0.0, -4.0, 1.0);

// Indirect lighting for Phong
vec3 indirectLight = 0.5f * vec3( 1, 1, 1 );

// Shadow Bias
const float shadowBiasThreshold = 0.001f;

// Standard geometry
vector<Triangle> triangles;
vector<Sphere> spheres;
vector<LightSource> lights;

// Phong geometry
vector<PhongTriangle> phongTriangles;
vector<PhongSphere> phongSpheres;
vector<PhongLightSource> phongLights;

// Photon maps
vector<Photon> photonMap;
vector<Photon> causticMap;

// Camera rotation
mat4 R;
float yaw = 0;
float pitch = 0;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);

void emitPhotons();
void emitPhotonsFromLight(LightSource &l, int NUM_PHOTONS);
bool tracePhoton(vec3 power, vec3 start, vec3 direction, int depth, bounce bounce);
void drawPhotons(screen* screen);

vec3 getPhongPixelValue(vec3 direction);
vec3 getPixelValue(vec3 start, vec3 direction, int depth);

vec3 phongComputeLight(const Intersection &i, const PhongLightSource &l);
vec3 getDirectLight(const Intersection& i, const LightSource& l);

vec3 sampleLightSource(const LightSource& l);
void sampleSquareLightSource(const LightSource& l, vector<vec3>& points);

void getNNearestPhotons(Intersection& intersection, vector<int>& indices, vector<Photon>& map, float& maxDistance);
vec3 getNearestPhotonsPower(Intersection& intersection, vector<Photon>& map);

vec3 reflect(vec3 dir, vec3 normal);
vec3 refract(vec3 dir, vec3 normal, Material material);
float fresnel(vec3 dir, vec3 normal, Material material);

bool closestIntersection(vec3 start, vec3 dir, Intersection& closestIntersection);
bool intersectTriangle(Intersection& closestIntersection, vec3 start, vec3 dir, vec3 v0, vec3 v1, vec3 v2, int index, vec3 normal);
bool intersectSphere(Intersection& intersection, vec3 start, vec3 dir, vec3 centre, float radius, int index);
bool intersectSquare(Intersection& intersection, vec3 start, vec3 dir, vec3 position, vec3 normal, float width, float length);

Material getMaterial(Intersection& intersection);

void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void updateRotation();
void moveCameraRight(int direction);
void moveCameraUp(int direction);
void moveCameraForward(int direction);
void lookAt(mat4& ctw);

int main(int argc, char* argv[]) {
  if (argc > 1) {
    try {
      NUM_PHOTONS = std::stoi(argv[1]);
      PHOTON_MAPPER = true;
      if (argc > 2) {
        try {
          NUM_NEAREST_PHOTONS = std::stoi(argv[2]);
        } catch (std::exception const &e) {
          cout << "Could not parse number of nearest photons." << endl;
        }
      }
    } catch (std::exception const &e) {
      cout << "Could not parse number of photons, using phong." << endl;
      PHOTON_MAPPER = false;
    }
  } else {
    PHOTON_MAPPER = false;
  }

  screen *mainscreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  screen *photonscreen = NULL;

  if (PHOTON_MAPPER) {
    cout << "Number of photons: " << NUM_PHOTONS << endl;
    cout << "Number of nearest photons in radiance estimate: " << NUM_NEAREST_PHOTONS << endl;
    photonscreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
    LoadTestModel(triangles, spheres, lights);
    //LoadBunny(triangles);
    emitPhotons();
  } else {
    LoadTestModelPhong(phongTriangles, phongSpheres, phongLights);
  }

  cameraPos = defaultCameraPos;

  while (Update()) {
    Draw(mainscreen);
    SDL_Renderframe(mainscreen);
    if (PHOTON_MAPPER) {
      drawPhotons(photonscreen);
      SDL_Renderframe(photonscreen);
    }
  }

  SDL_SaveImage(mainscreen, "mainout.bmp");
  KillSDL(mainscreen);

  if (PHOTON_MAPPER) {
    SDL_SaveImage(photonscreen, "photonmapout.bmp");
    KillSDL(photonscreen);
  }

	return 0;
}

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {

      vec3 d = normalize(vec3(R * vec4(x - SCREEN_WIDTH/2, y - SCREEN_HEIGHT/2, focalLength, 1)));

      if (PHOTON_MAPPER) {
        PutPixelSDL(screen, x, y, getPixelValue(vec3(cameraPos), d, 0));
      } else {
        PutPixelSDL(screen, x, y, getPhongPixelValue(d));
      }
    }
  }
}

vec3 getPhongPixelValue(vec3 direction) {
  vec3 colour(0.0,0.0,0.0);
  vec3 light(0.0,0.0,0.0);

  Intersection intersection;

  if (closestIntersection(vec3(cameraPos), direction, intersection)) {
    if (intersection.intersectionType == triangle) {
      colour = phongTriangles[intersection.index].material.color;
    } else if (intersection.intersectionType == sphere) {
      colour = phongSpheres[intersection.index].material.color;
    }

    for (size_t i = 0; i < phongLights.size(); i++) {
      light += phongComputeLight(intersection, phongLights[i]);
    }

    light = colour * (light + indirectLight);
  }

  return light;
}

vec3 getPixelValue(vec3 start, vec3 direction, int depth) {
  if (depth > MAX_CAMERA_RAY_DEPTH) return vec3(0.0,0.0,0.0);

  vec3 colour(0.0,0.0,0.0);
  vec3 emmittedLight(0.0,0.0,0.0);
  vec3 diffuseLight(0.0,0.0,0.0);
  vec3 causticLight(0.0,0.0,0.0);
  vec3 directLight(0.0,0.0,0.0);

  Intersection lightIntersection;
  Intersection intersection;

  for (size_t i = 0; i < lights.size(); i++) {
    LightSource l = lights[i];

    if (intersectSquare(lightIntersection, start, direction, l.position, l.direction, l.width, l.length)) {
      emmittedLight += l.watts * l.color;
    }
  }

  if (closestIntersection(start, direction, intersection)) {
    Material material = getMaterial(intersection);

    // If the surface is a mirror
    if (glm::equal(material.specRef, vec3(1.0f,1.0f,1.0f))[0]) {
      vec3 reflectDir = reflect(direction, intersection.normal);

      return getPixelValue(intersection.position, reflectDir, depth + 1);
    }

    if (material.refractiveIndex > 0.0f) {
      vec3 refractDir = refract(direction, intersection.normal, material);
      vec3 reflectDir = reflect(direction, intersection.normal);

      vec3 refractedColour(0.0,0.0,0.0);
      vec3 reflectedColour(0.0,0.0,0.0);

      float fresnelCoeff = fresnel(direction, intersection.normal, material);

      if (fresnelCoeff < 1.0f) {
        refractedColour = getPixelValue(intersection.position, refractDir, depth + 1);
      }

      reflectedColour = getPixelValue(intersection.position, reflectDir, depth + 1);

      return fresnelCoeff * reflectedColour + (1 - fresnelCoeff) * refractedColour;
    }

    for (size_t i = 0; i < lights.size(); i++) {
      directLight += getDirectLight(intersection, lights[i]);
    }

    causticLight += getNearestPhotonsPower(intersection, causticMap);
    diffuseLight += getNearestPhotonsPower(intersection, photonMap);
  } else {
    return vec3(0.0f,0.0f,0.0f);
  }

  colour = causticLight + directLight + diffuseLight + emmittedLight;

  return colour * (1.0f - max(0.0f, dot(normalize(intersection.normal), normalize(intersection.direction))));
}

vec3 reflect(vec3 dir, vec3 normal) {
  dir = normalize(dir);
  normal = normalize(normal);

  return normalize(dir - 2.0f * dot(dir, normal) * normal);
}

vec3 refract(vec3 dir, vec3 normal, Material material) {
  dir = normalize(dir);
  normal = normalize(normal);

  float n1 = (float) GLOBAL_REF_INDEX;
  float n2 = material.refractiveIndex;

  float c1 = dot(normal, dir);

  if (c1 > 0.0f) {
    // inside the material
    n1 = material.refractiveIndex;
    n2 = (float) GLOBAL_REF_INDEX;
    normal = -normal;
  } else {
    // outside the material
    c1 = -c1;
  }

  float n = (float) n1 / n2;
  float s = 1.0f - (float) (n * n * (1.0f - c1 * c1));

  if (s < 0.0f) {
    // total internal reflection
    return reflect(dir, normal);
  } else {
    return normalize(n * dir + (n * c1 - (float) sqrtf(s)) * normal);
  }
}

// Calculate the fresnel coefficient
float fresnel(vec3 dir, vec3 normal, Material material) {
  dir = normalize(dir);
  normal = normalize(normal);

  float n1 = (float) GLOBAL_REF_INDEX;
  float n2 = material.refractiveIndex;

  float c1 = dot(normal, dir);

  if (c1 > 0.0f) {
    // inside the material
    n1 = material.refractiveIndex;
    n2 = (float) GLOBAL_REF_INDEX;
    normal = -normal;
  } else {
    // outside the material
    c1 = -c1;
  }

  float n = (float) n1 / n2;

  float sint = n * sqrtf(max(0.0f, 1.0f - c1 * c1));

  if (sint >= 1) {
    // total internal reflection
    return 1.0f;
  } else {
    float cost = sqrtf(max(0.0f, 1.0f - sint * sint));
    c1 = fabsf(c1);
    float Rs = ((n2 * c1) - (n1 * cost)) / ((n2 * c1) + (n1 * cost));
    float Rp = ((n1 * c1) - (n2 * cost)) / ((n1 * c1) + (n2 * cost));
    return (float) (Rs * Rs + Rp * Rp) / 2.0f;
  }
}

vec3 sampleLightSource(const LightSource& l) {
  vec3 position(1.0f,1.0f,1.0f);

  position.x = ((float) rand() / (RAND_MAX)) * l.width + (l.position.x - l.width/2);
  position.y = l.position.y + shadowBiasThreshold;
  position.z = ((float) rand() / (RAND_MAX)) * l.length + (l.position.z - l.length/2);

  return position;
}

void sampleSquareLightSource(const LightSource& l, vector<vec3>& points) {
  points.clear();

  int gridWidth = sqrt(NUM_SHADOW_RAYS);

  // Create a randomised grid
  for (int i = 0; i < gridWidth; i++) {
    for (int j = 0; j < gridWidth; j++) {
      float x = (l.position.x - l.width/2 + (l.width/gridWidth)/2) + (i * (l.width/gridWidth))
                + (((float) rand() / (RAND_MAX)) * l.width/(gridWidth * 4) - l.width/(gridWidth * 8));
      float y = l.position.y + 0.01f;
      float z = (l.position.z - l.length/2 + (l.length/gridWidth)/2) + (j * (l.length/gridWidth))
                + (((float) rand() / (RAND_MAX)) * l.length/(gridWidth * 4) - l.length/(gridWidth * 8));

      points.push_back(vec3(x,y,z));
    }
  }
}

void getNNearestPhotons(Intersection& intersection, vector<int>& indices, vector<Photon>& map, float& maxDistance) {
  vector<float> distances;

  indices.clear();

  // Iterate over all photons in the map
  for (size_t i = 0; i < map.size(); i++) {
    vec3 ipos = intersection.position;
      vec3 ppos = map[i].position;

      float dist = distance(ipos, ppos);

      if ((int) indices.size() < NUM_NEAREST_PHOTONS) {
        // There are fewer than needed photons in the list so add
        indices.push_back(i);
        distances.push_back(dist);
      } else {
        // Find the furthest photon currently in the list and replace it with this photon
        // if it is closer
        int furthestIndex = 0;
        float maxDist = -numeric_limits<float>::max();

        for (int j = 0; j < NUM_NEAREST_PHOTONS; j++) {
          if (distances[j] > maxDist) {
            maxDist = distances[j];
            furthestIndex = j;
          }
        }

        maxDistance = maxDist;

        if (dist < maxDist) {
          indices[furthestIndex] = i;
          distances[furthestIndex] = dist;
        }
      }
  }
}

vec3 getNearestPhotonsPower(Intersection& intersection, vector<Photon>& map) {
  vector<int> nearestPhotonsIndex;
  vec3 accumPower = vec3(0.0f,0.0f,0.0f);
  float radius = 0.0f;
  float dist = 0.0f;
  vec3 photonPos(1.0,1.0,1.0);

  // Retrieve the nearest photons
  getNNearestPhotons(intersection, nearestPhotonsIndex, map, radius);

  for (size_t i = 0; i < nearestPhotonsIndex.size(); i++) {
    photonPos = map[nearestPhotonsIndex[i]].position;
    dist = distance(intersection.position, photonPos);

    // Weight each photons power on its distance from the intersection
    float wpc = 1.0f - (dist / (FILTER_CONSTANT * radius));

    accumPower += (map[nearestPhotonsIndex[i]].power * wpc);
  }

  // Calculate the unit power by dividing by the area of the projected sphere
  vec3 unitPower = accumPower / (float) ((1.0f - FILTER_CONSTANT * (float) 2/3) * (float) M_PI * pow(radius, 2.0f));

  return unitPower;
}

// Emit photons from all light sources
void emitPhotons() {
  int numLights = lights.size();

  cout << "-------------------" << endl;
  cout << "emitting photons from lights" << endl;

  for (int i = 0; i < numLights; i++) {
    emitPhotonsFromLight(lights[i], NUM_PHOTONS / numLights);
  }

  cout << "finished emitting photons from lights" << endl;
  cout << "-------------------" << endl;
  cout << "global photon map size: " << photonMap.size() << endl;
  cout << "caustics photon map size: " << causticMap.size() << endl;
  cout << "-------------------" << endl;
}

// Emit numPhotons from the light source
void emitPhotonsFromLight(LightSource &l, int numPhotons) {

  for(int i = 0; i < numPhotons; i++) {
    vec3 direction = -l.direction;
    vec3 position = sampleLightSource(l);

    vec3 direction3d(direction.x, direction.y, direction.z);
    vec3 lightDirection3d(l.direction.x, l.direction.y, l.direction.z);

    while (dot(normalize(direction3d), normalize(lightDirection3d)) < cos((float) SPOTLIGHT_CUTOFF)) {
      direction3d.x = ((float) rand() / (RAND_MAX)) * 2 - 1;
      direction3d.y = ((float) rand() / (RAND_MAX)) * 2 - 1;
      direction3d.z = ((float) rand() / (RAND_MAX)) * 2 - 1;
    }

    tracePhoton((l.watts / numPhotons) * l.color, position, direction3d, 0, none);
  }
}

// Recursive function to trace a photon in the scene
bool tracePhoton(vec3 power, vec3 start, vec3 direction, int depth, bounce bounce) {
  if (depth > MAX_PHOTON_DEPTH) return false;

  Intersection intersection;

  if (closestIntersection(start, direction, intersection)) {
    vec3 color(0.0,0.0,0.0);
    vec3 normal = intersection.normal;
    Material material(color,color, 0.0);

    // Get the material parameters
    if (intersection.intersectionType == triangle) {
      Triangle triangle = triangles[intersection.index];
      material = triangle.material;
      color = triangle.color;
    } else if(intersection.intersectionType == sphere) {
      Sphere sphere = spheres[intersection.index];
      material = sphere.material;
      color = sphere.color;
    }

    vec3 diffuseRef = material.diffuseRef;
    vec3 specRef = material.specRef;

    // Probability of reflection
    float Pr = max(diffuseRef.x + specRef.x, max(diffuseRef.y + specRef.y, diffuseRef.z + specRef.z));

    // Probability of diffuse reflection
    float Pd = Pr * ((diffuseRef.x + diffuseRef.y + diffuseRef.z) /
                      (diffuseRef.x + diffuseRef.y + diffuseRef.z + specRef.x + specRef.y + specRef.z));

    // Probability of specular reflection
    float Ps = Pr - Pd;

    vec3 reflectionDir = reflect(direction, normal);
    vec3 specPower = vec3(power.x * specRef.x / Ps, power.y * specRef.y / Ps, power.z * specRef.z / Ps);

    float rnd = ((float) rand() / RAND_MAX);

    if (material.refractiveIndex > 0.0f) {
      // Transmit photon
      vec3 refractDir = refract(direction, normal, material);

      tracePhoton(power, intersection.position, refractDir, depth + 1, specular);
    } else {
      if (rnd < Pd) {
        //Diffuse reflection
        tracePhoton(power * color, intersection.position, reflectionDir, depth + 1, diffuse);
      } else if (rnd < Ps + Pd) {
        // Specular reflection
        tracePhoton(specPower * color, intersection.position, reflectionDir, depth + 1, specular);
      } else {
        // Absorbtion
        if (depth > 0) {
          Photon p;
          p.position = intersection.position;
          p.direction = direction;
          p.power = power;

          if (bounce == specular) causticMap.push_back(p);
          else photonMap.push_back(p);
        }
      }
    }
  }

  return false;
}

// Get the direct light component for a given intersection and lightsource
vec3 getDirectLight(const Intersection& i, const LightSource& l) {
  vec3 nHat;
  vec3 color;
  vec3 directLight;

  if (i.intersectionType == triangle) {
    Triangle triangle = triangles[i.index];

    nHat = normalize(triangle.normal);
    color = triangle.color;
  } else if (i.intersectionType == sphere) {
    Sphere sphere = spheres[i.index];

    nHat = normalize(i.position - sphere.centre);
    color = sphere.color;
  }

  vec3 lightDir = l.position - i.position;

  vec3 rHat = normalize(lightDir);

  float r = distance(i.position, l.position);

  float A = 4 * (float) M_PI * pow(r, 2);

  vec3 B = l.color * l.watts / A;
  vec3 D = B * ((dot(rHat, normalize(-l.direction)) > cos((float) SPOTLIGHT_CUTOFF))
                ? (min(0.5f, max(dot(rHat,normalize(-l.direction)), 0.0f))) : 0.0f);

  vec3 C = D * color;

  directLight += C;

  // Soft Shadows
  vector<vec3> lightPoints;

  sampleSquareLightSource(l, lightPoints);
  int hitsLight = lightPoints.size();

  for (size_t h = 0; h < lightPoints.size(); h++) {
    vec3 lightDir = lightPoints[h] - i.position;

    vec3 rHat = normalize(lightDir);

    float r = distance(i.position, lightPoints[h]);

    Intersection intersection;

    if (closestIntersection(i.position, rHat, intersection)) {
      if (intersection.distance < r) {
        hitsLight--;
      } else {
        A = 4 * (float) M_PI * pow(r, 2);

        B = l.color * l.watts / A;
        D = B * ((dot(rHat, normalize(-l.direction)) > cos((float) SPOTLIGHT_CUTOFF))
                      ? (min(0.5f, max(dot(rHat,normalize(-l.direction)), 0.0f))) : 0.0f);
        C = D * color;

        directLight += C;
      }
    }
  }

	return ((float) hitsLight / lightPoints.size()) * (directLight / (float) (hitsLight + 1));
}

// Calculate the phong lighting for a given intersection and light source
// Using equation from https://en.wikipedia.org/wiki/Phong_reflection_model
vec3 phongComputeLight( const Intersection &i, const PhongLightSource &l ) {
  vec3 light(0.0,0.0,0.0);

  vec3 Lm = normalize(l.position - i.position);
  vec3 V = normalize(vec3(cameraPos) - i.position);

  float dist = distance(l.position, i.position);

  float id = l.diffuseIntensity;
  float is = l.specularIntensity;
  float ia = l.ambientIntensity;

  float kd = 1;
  float ks = 1;
  float ka = 1;
  float alpha;

  vec3 normal;
  vec3 Rm;

  float LmNormalDot;
  float RmVDot;

  if (i.intersectionType == triangle) {
    PhongTriangle triangle = phongTriangles[i.index];

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
  } else if (i.intersectionType == sphere) {
    PhongSphere sphere = phongSpheres[i.index];

    ka = sphere.material.ambientRef;

    light += ka * (ia * l.color);

    kd = sphere.material.diffuseRef;

    normal = normalize(i.position - sphere.centre);

    LmNormalDot = dot(Lm, normal);

    ks = sphere.material.specularRef;

    Rm = normalize(2 * max(0.0f, LmNormalDot) * normal - Lm);

    alpha = sphere.material.shininess;

    RmVDot = dot(Rm, V);

    if (LmNormalDot > 0) {
      light += (kd * max(0.0f, LmNormalDot) * (id * l.color)) / dist;

      if (RmVDot > 0) {
        light += (ks * pow(RmVDot, alpha) * (is * l.color)) / dist;
      }
    }
  }

  Intersection intersection;

  vec3 lightDir = l.position - i.position;
  vec3 rHat = normalize(lightDir);
  float r = distance(i.position, l.position);

  if (closestIntersection(i.position, rHat, intersection)) {
    if (intersection.index != i.index && intersection.distance < r) {
      light = ka * (ia * l.color);
    }
  }

  return light;
}

// Intersect triangle
bool intersectTriangle(Intersection& closestIntersection, vec3 start, vec3 dir, vec3 v0, vec3 v1, vec3 v2, int index, vec3 normal) {
  bool intersectionFound = false;

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
      vec3 position = start + t * dir;
      float dist = distance(start, position);

      if (dist <= closestIntersection.distance && dist > shadowBiasThreshold) {
        intersectionFound = true;
        closestIntersection.distance = dist;
        closestIntersection.position = position;
        closestIntersection.normal = normal;
        closestIntersection.index = index;
        closestIntersection.intersectionType = triangle;
        closestIntersection.direction = dir;
      }
    }
  }

  return intersectionFound;
}

// Intersect sphere
bool intersectSphere(Intersection& closestIntersection, vec3 start, vec3 dir, vec3 ce, float ra, int index) {
  bool intersectionFound = false;

  vec3 oc = start - ce;
  float a = dot(dir, dir);
  float b = 2.0f * dot(oc, dir);
  float c = dot(oc,oc) - ra*ra;
  float discriminant = b*b - 4.0f*a*c;

  if (discriminant < 0.0f) {
    return false;
  } else {
    float t = (-b - sqrt(discriminant)) / (2.0f*a);
    if (t >= 0) {
      vec3 position = start + t * dir;
      float dist = distance(start, position);

      if (dist <= closestIntersection.distance && dist > shadowBiasThreshold) {
        intersectionFound = true;
        closestIntersection.distance = dist;
        closestIntersection.position = position;
        closestIntersection.normal = normalize(position - ce);
        closestIntersection.index = index;
        closestIntersection.intersectionType = sphere;
        closestIntersection.direction = dir;
      }
    }
  }

  return intersectionFound;
}

// Intersect a square
// Requires axis aligned square at the moment
bool intersectSquare(Intersection& intersection, vec3 start, vec3 dir, vec3 position, vec3 normal, float width, float length) {
  bool intersectionFound = false;

  vec3 corner(position.x - width/2, position.y, position.z - length/2);
  float t = dot((normalize(corner) - normalize(start)), normalize(normal))
            / dot(normalize(dir), normalize(normal));

  if (t >= 0) {
    vec3 pos = start + t * dir;

    float dist = distance(start, pos);

    vec3 v = pos - corner;

    if (v.x >= 0 && v.x <= width && v.z >= 0 && v.z <= length) {
      intersectionFound = true;
      intersection.distance = dist;
      intersection.position = pos;
      intersection.direction = dir;
    }
  }

  return intersectionFound;
}

bool closestIntersection(vec3 start, vec3 dir, Intersection& closestIntersection) {
  bool intersectionFound = false;
  closestIntersection = Intersection();
  closestIntersection.distance = std::numeric_limits<float>::max();

  int trianglesSize = (PHOTON_MAPPER) ? triangles.size() : phongTriangles.size();

  for(int i = 0; i < trianglesSize; i++) {
    vec3 v0 = (PHOTON_MAPPER) ? triangles[i].v0 : phongTriangles[i].v0;
    vec3 v1 = (PHOTON_MAPPER) ? triangles[i].v1 : phongTriangles[i].v1;
    vec3 v2 = (PHOTON_MAPPER) ? triangles[i].v2 : phongTriangles[i].v2;
    vec3 normal = (PHOTON_MAPPER) ? triangles[i].normal : phongTriangles[i].normal;

    intersectionFound = intersectionFound | intersectTriangle(closestIntersection, start, dir, v0, v1, v2, i, normal);
  }

  int spheresSize = (PHOTON_MAPPER) ? spheres.size() : phongSpheres.size();

  for(int i = 0; i < spheresSize; i++) {
    vec3 centre = (PHOTON_MAPPER) ? spheres[i].centre : phongSpheres[i].centre;
    float radius = (PHOTON_MAPPER) ? spheres[i].radius : phongSpheres[i].radius;

    intersectionFound = intersectionFound | intersectSphere(closestIntersection, start, dir, centre, radius, i);
  }

  return intersectionFound;
}

// Return the material at the intersection
Material getMaterial(Intersection& intersection) {
  if (intersection.intersectionType == triangle) {
    return triangles[intersection.index].material;
  } else {
    return spheres[intersection.index].material;
  }
}

// Draw the photons on a new screen
void drawPhotons(screen* screen) {
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // Display the global photons
  for (size_t i = 0; i < photonMap.size(); i++) {
    Photon p = photonMap[i];

    vec4 pos = R * vec4(p.position,1) - cameraPos;

    int x = (focalLength * pos.x) / pos.z + SCREEN_WIDTH/2;
    int y = (focalLength * pos.y) / pos.z + SCREEN_WIDTH/2;
    if (x > 0 && x < SCREEN_WIDTH && y > 0 && y < SCREEN_HEIGHT) PutPixelSDL(screen, x, y, normalize(p.power));
  }

  // Display the caustic photons
  for (size_t i = 0; i < causticMap.size(); i++) {
    Photon p = causticMap[i];

    vec4 pos = R * vec4(p.position,1) - cameraPos;

    int x = (focalLength * pos.x) / pos.z + SCREEN_WIDTH/2;
    int y = (focalLength * pos.y) / pos.z + SCREEN_WIDTH/2;
    if (x > 0 && x < SCREEN_WIDTH && y > 0 && y < SCREEN_HEIGHT) PutPixelSDL(screen, x, y, normalize(p.power));
  }

  // Show the light sources
  /* for (int i = 0; i < lights.size(); i++) {
    float width = lights[i].width;
    float length = lights[i].length;

    for (float j = -width/2; j <= width/2; j += width/100) {
      for (float k = -length/2; k <= length/2; k += length/100) {
        vec4 pixelPos(lights[i].position.x + j, lights[i].position.y, lights[i].position.z + k, 1);
        vec4 pos = R * pixelPos - cameraPos;

        int x = (focalLength * pos.x) / pos.z + SCREEN_WIDTH/2;
        int y = (focalLength * pos.y) / pos.z + SCREEN_WIDTH/2;
        if (x > 0 && x < SCREEN_WIDTH && y > 0 && y < SCREEN_HEIGHT) PutPixelSDL(screen, x, y, lights[i].color);
      }
    }
  }*/

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
					pitch += (float) M_PI / 18;
					updateRotation();
					break;
	      case SDLK_DOWN:
					pitch -= (float) M_PI / 18;
					updateRotation();
          break;
	      case SDLK_LEFT:
					yaw -= (float) M_PI / 18;
					updateRotation();
          break;
	      case SDLK_RIGHT:
					yaw += (float) M_PI / 18;
					updateRotation();
          break;
				case SDLK_r:
					/* Look-At function, points camera to 0,0,0 */
					lookAt(ctw);
					break;
        case SDLK_t:
          // Reset camera position
          cameraPos = defaultCameraPos;
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
	      case SDLK_ESCAPE:
          /* Move camera quit */
          return false;
      }
    }
  }
  return true;
}
