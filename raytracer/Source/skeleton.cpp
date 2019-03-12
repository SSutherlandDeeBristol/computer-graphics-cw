#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>
#include "limits"
#include <math.h>

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

enum type {triangle, sphere};

struct Intersection {
  vec4 position;
  float distance;
  int index;
  type type;
};

const float focalLength = SCREEN_HEIGHT;
const vec4 defaultCameraPos(0.0, 0.0, -3.0, 1.0);
vec4 cameraPos(0.0, 0.0, -3.0, 1.0);

// const vec4 defaultLightPos( 0, -0.5, -0.7, 1.0 );
// vec4 lightPos( 0, -0.5, -0.7, 1.0 );
// vec3 lightColor = 14.f * vec3( 1, 1, 1 );
vec3 indirectLight = 0.3f * vec3( 1, 1, 1 );
const float shadowBiasThreshold = 0.001f;

std::vector<Triangle> triangles;
std::vector<Sphere> spheres;
std::vector<Light> lights;

mat4 R;
float yaw = 0;
float pitch = 0;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, const vector<Sphere>& spheres, Intersection& closestIntersection);
bool intersectTriangle(const Triangle &triangle, vec4 start, vec4 dir, float &dist, vec4 &position);
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void updateRotation();
vec3 computeLight( const Intersection &i, const Light &l );
vec3 DirectLight( const Intersection &i, const Light &l );
void moveCameraRight(int direction);
void moveCameraUp(int direction);
void moveCameraForward(int direction);
void lookAt(mat4& ctw);

int main(int argc, char* argv[]) {
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  LoadTestModel(triangles, spheres, lights);

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

      vec3 light(0.0, 0.0, 0.0); // The direct colour
			vec3 colour(0.0, 0.0, 0.0); // The original colour of the triangle

      if (ClosestIntersection(cameraPos, d, triangles, spheres, intersection)) {
				switch(intersection.type) {
          case triangle :
            colour = triangles[intersection.index].material.color;

            //light = triangles[intersection.index].material.ambientRef * indirectLight;
            break;
          case sphere :
            colour = spheres[intersection.index].material.color;

            //light = spheres[intersection.index].material.ambientRef * indirectLight;
            break;
          default :
            break;
        }

        for (int i = 0; i < lights.size(); i++) {
          //light += computeLight(intersection, lights[i]);
          light += DirectLight(intersection, lights[i]);
        }

				light = colour * light;
      }

      PutPixelSDL(screen, x, y, light);
    }
  }
}

bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, const vector<Sphere>& spheres, Intersection& closestIntersection) {
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
					closestIntersection.index = i;
          closestIntersection.type = triangle;
				}
			}
    }
  }

  for(std::vector<Sphere>::size_type i = 0; i < spheres.size(); i++) {
    vec3 ro = vec3(start);
    vec3 rd = vec3(dir);

    vec3 ce = vec3(spheres[i].centre);
    float ra = spheres[i].radius;

    vec3 oc = ro - ce;
    float b = dot( oc, rd );
    float c = dot( oc, oc ) - ra*ra;
    float h = b*b - c;

    if( h < 0.0 ) {
      // no intersection
      break;
    }

    h = sqrt( h );

    // -b - h or -b + h
    for(int sign = -1; sign <= 1; sign += 2) {
      // To prevent doing the same sum twice
      if (!(sign == 1 && h == 0)) {
        float sol = -b + sign * h;
        vec4 position = start + sol * dir;
        float dist = distance(start, position);

        if (dist <= closestIntersection.distance && dist > shadowBiasThreshold) {
          intersectionFound = true;
          closestIntersection.distance = dist;
          closestIntersection.position = position;
          closestIntersection.index = i;
          closestIntersection.type = sphere;
        }
      }
    }
  }

  return intersectionFound;
}

vec3 computeLight( const Intersection &i, const Light &l ) {
  // Using equation from https://en.wikipedia.org/wiki/Phong_reflection_model
  vec3 light(0.0,0.0,0.0);

  vec4 Lm = normalize(l.position - i.position);
  vec4 V = normalize(cameraPos - i.position);

  float id = l.diffuseIntensity;
  float is = l.specularIntensity;

  float kd;
  float ks;
  float alpha;

  vec4 normal;
  vec4 Rm;

  if (i.type == triangle) {
    Triangle triangle = triangles[i.index];

    kd = triangle.material.diffuseRef;

    normal = normalize(triangle.normal);

    light += kd * dot(Lm, normal) * id;

    ks = triangle.material.specularRef;

    Rm = normalize(2 * dot(Lm, normal) * normal - Lm);

    alpha = triangle.material.shininess;

    light += ks * pow(dot(Rm, V), alpha) * is;

  } else if (i.type == sphere) {
    Sphere sphere = spheres[i.index];

    kd = sphere.material.diffuseRef;

    normal = normalize(i.position - sphere.centre);

    light += kd * dot(Lm, normal) * id;

    ks = sphere.material.specularRef;

    Rm = normalize(2 * dot(Lm, normal) * normal - Lm);

    alpha = sphere.material.shininess;

    light += ks * pow(dot(Rm, V), alpha) * is;
  }

  return light;
}

vec3 DirectLight( const Intersection &i, const Light &l ) {
  vec3 color = vec3(1,1,1);
  vec3 black = vec3(0.0,0.0,0.0);

  switch (i.type) {
    case triangle: {
      Triangle triangle = triangles[i.index];

      float r = distance(i.position, l.position);

      float A = 4 * M_PI * pow(r, 2);

      vec4 lightDir = l.position - i.position;

      vec4 rHat = normalize(lightDir);
      vec4 nHat = normalize(triangle.normal);

      vec3 B = l.color / A;
      vec3 D = B * max(dot(rHat,nHat), 0.0f);

      color = l.diffuseIntensity * D;

      Intersection intersection;

      if (ClosestIntersection(i.position, rHat, triangles, spheres, intersection)) {
        if (intersection.index != i.index && intersection.distance < r) {
          color = black;
        }
      }

      break;
    }
    case sphere: {
      Sphere sphere = spheres[i.index];

      float r = distance(i.position, l.position);

      float A = 4 * M_PI * pow(r, 2);

      vec4 lightDir = l.position - i.position;

      vec4 rHat = normalize(lightDir);
      vec4 nHat = normalize(i.position - sphere.centre);

      vec3 B = l.color / A;
      vec3 D = B * max(dot(rHat,nHat), 0.0f);

      color = l.diffuseIntensity * D;

      Intersection intersection;

      if (ClosestIntersection(i.position, rHat, triangles, spheres, intersection)) {
        if (intersection.index != i.index && intersection.distance < r) {
          color = black;
        }
      }

      break;
    }
    default:
      break;
  }

  return color;
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
          //lightPos = defaultLightPos;
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
					//lightPos.z += 0.5;
					break;
				case SDLK_k:
					//lightPos.z -= 0.5;
					break;
				case SDLK_j:
					//lightPos.x -= 0.5;
					break;
				case SDLK_l:
					//lightPos.x += 0.5;
					break;
				case SDLK_o:
					//lightPos.y -= 0.5;
					break;
				case SDLK_p:
					//lightPos.y += 0.5;
					break;
	      case SDLK_ESCAPE:
          /* Move camera quit */
          return false;
      }
    }
  }
  return true;
}
