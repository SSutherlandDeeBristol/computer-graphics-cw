#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>
#include "limits"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::distance;

SDL_Event event;

#define SCREEN_WIDTH 160
#define SCREEN_HEIGHT 128
#define FULLSCREEN_MODE false

struct Intersection {
  vec4 position;
  float distance;
  int triangleIndex;
};

float focalLength = SCREEN_WIDTH/2;
vec4 cameraPos(0.0, 0.0, -2, 1.0);
std::vector<Triangle> triangles;
mat4 R;
float yaw;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);

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
    //vec3 x = glm::inverse(A) * b;
      
    mat3 A1(b, e1, e2);
    float detA = glm::determinant(A);
      
    float t = glm::determinant(A1) / detA;
      
    if (t >= 0) {
        mat3 A2(-vec3(dir), b, e2);
        mat3 A3(-vec3(dir), e1, b);
          
        float u = glm::determinant(A2) / detA;
        float v = glm::determinant(A3) / detA;
          
        if (u > 0 && v > 0 && (u + v) < 1) {
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

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
      vec4 d = vec4(x - SCREEN_WIDTH/2, y - SCREEN_HEIGHT/2, focalLength, 1);
      Intersection intersection;
  
      vec3 colour(0.0, 0.0, 0.0); // Initialise to black

      if (ClosestIntersection(cameraPos, d, triangles, intersection)) {
        colour = triangles[intersection.triangleIndex].color;
      }

      PutPixelSDL(screen, x, y, colour);
    }
  }
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
                cameraPos.y -= 0.05;
          /* Move camera forward */
          break;
	      case SDLK_DOWN:
                cameraPos.y += 0.05;
          /* Move camera backwards */
          break;
	      case SDLK_LEFT:
                cameraPos.x -= 0.05;
          /* Move camera left */
          break;
	      case SDLK_RIGHT:
                cameraPos.x += 0.05;
          /* Move camera right */
          break;
	      case SDLK_ESCAPE:
          /* Move camera quit */
          return false;
      }
    }
  }
  return true;
}

void updateRotation() {
    
}
