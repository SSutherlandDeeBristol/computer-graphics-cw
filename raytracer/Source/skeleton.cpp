#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>
#include "limits"
#include <math.h>
#include <unistd.h>
#include "glm/gtx/string_cast.hpp" // std::cout<<glm::to_string(hello)<<std::endl;

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::distance;

SDL_Event event;

#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 500
#define FULLSCREEN_MODE false

struct Intersection {
  vec4 position;
  float distance;
  int triangleIndex;
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
std::vector<Sphere> spheres;

mat4 R;
float yaw = 0;
float pitch = 0;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, const vector<Sphere>& spheres, Intersection& closestIntersection);
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void updateRotation();
vec3 DirectLight( const Intersection& i );
void moveCameraRight(int direction);
void moveCameraUp(int direction);
void moveCameraForward(int direction);
void lookAt(mat4& ctw);
void readBunny(std::vector<Triangle>& triangles);

int main(int argc, char* argv[]) {
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  LoadTestModel(triangles, spheres);
	readBunny(triangles);

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

      if (ClosestIntersection(cameraPos, d, triangles, spheres, intersection)) {
				directLight = DirectLight(intersection);
				colour = triangles[intersection.triangleIndex].color;

				reflectedLight = colour * (directLight + indirectLight);
      }

      PutPixelSDL(screen, x, y, reflectedLight);
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
					closestIntersection.triangleIndex = i;
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

	if (ClosestIntersection(i.position, rHat, triangles, spheres, intersection)) {
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

	// std::cout << "cx: " << cameraPos.x << ", cy:" << cameraPos.y << ", cz:"<< cameraPos.z << std::endl;

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

void readBunny(std::vector<Triangle>& triangles) {
	std::ifstream infile("Resources/bunny.obj");
	std::string line;

	vec3 white(  0.75f, 0.75f, 0.75f );

	/* Vector to store all vertices */
	vector<vec4> vs;

	while (std::getline(infile, line)) {

		/* Try parse line */
		std::istringstream iss(line);
		string code;
		float a, b, c;
		if (!(iss >> code >> a >> b >> c)) {
			cout << "Could not read file" << endl;
			break;
		}

		cout << code << ", " + to_string(a) << ", " + to_string(b) << ", " + to_string(c) << endl;

		if (code.compare("v") == 0) {
			/* Then this is a vertex definition */
			vec4 vertex = vec4(a, b, c, 1.0);
			vs.push_back(vertex);
		} else {
			/* This is a face definition */
			Triangle triangle = Triangle(vs[a-1], vs[b-1], vs[c-1], white);
			float L = 0.5;
			triangle.v0 *= 2/L;
			triangle.v1 *= 2/L;
			triangle.v2 *= 2/L;

			triangle.v0 -= vec4(0,0.5,0.5,1);
			triangle.v1 -= vec4(0,0.5,0.5,1);
			triangle.v2 -= vec4(0,0.5,0.5,1);

			triangle.v0.x *= -1;
			triangle.v1.x *= -1;
			triangle.v2.x *= -1;

			triangle.v0.y *= -1;
			triangle.v1.y *= -1;
			triangle.v2.y *= -1;

			triangle.v0.w = 1.0;
			triangle.v1.w = 1.0;
			triangle.v2.w = 1.0;

			triangle.ComputeNormal();
			triangle.ReverseNormal();
			triangles.push_back(triangle);
		}

	}

}