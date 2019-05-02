#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include <stdint.h>
#include "limits"
#include <math.h>
#include "glm/gtx/string_cast.hpp"
#include <tuple>
#include <vector>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::distance;
using std::vector;

SDL_Event event;

#define SCREEN_WIDTH 600
#define SCREEN_HEIGHT 600
#define FULLSCREEN_MODE false

struct Intersection {
  vec3 position;
  vec3 direction;
  float distance;
  int index;
};

struct Line {
  vec3 start;
  vec3 end;
};

vector<Line> sceneGeometry;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void DrawScene(screen* screen);

bool closestIntersection(vec3 start, vec3 dir, Intersection& closestIntersection);

int main(int argc, char* argv[]) {
  screen *mainscreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  screen *scenescreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  cameraPos = defaultCameraPos;

  while (Update()) {
    Draw(mainscreen);
    SDL_Renderframe(mainscreen);
    SDL_SaveImage(mainscreen, "mainout.bmp");

    DrawScene(scenescreen);
    SDL_Renderframe(scenescreen);
  }

  KillSDL(mainscreen);
  KillSDL(scenescreen);

	return 0;
}

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  float lastpercentage = 0.0f;
  float screenSize = (float) (SCREEN_WIDTH * SCREEN_HEIGHT);

  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {

      vec3 d = normalize(vec3(R * vec4(x - SCREEN_WIDTH/2, y - SCREEN_HEIGHT/2, focalLength, 1)));

    }
  }
}

void DrawScene(screen* screen) {

}

bool closestIntersection(vec3 start, vec3 dir, Intersection& closestIntersection) {
  bool intersectionFound = false;
  closestIntersection = Intersection();
  closestIntersection.distance = std::numeric_limits<float>::max();

  return intersectionFound;
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
