#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
glm::mat4x4 TransformationMatrix(glm::vec4 camPos, glm::mat3x3 rot);
void Interpolate( glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result );

int main( int argc, char* argv[] ) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  while ( Update() ) {
    Draw(screen);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
}

glm::mat4x4 TransformationMatrix(glm::vec4 camPos, glm::mat3x3 rot) {
  vec4 zeroVec(0,0,0,0);
  glm::mat4x4 m1(zeroVec, zeroVec, zeroVec, camPos);

  vec4 rot0(rot[0], 0);
  vec4 rot1(rot[1], 0);
  vec4 rot2(rot[2], 0);

  glm::mat4x4 m2(rot0, rot1, rot2, vec4(0,0,0,1));

  vec4 minusCamPos(-camPos[0], -camPos[1], -camPos[2], 1);
  glm::mat4x4 m3(zeroVec, zeroVec, zeroVec, minusCamPos);

  return ((m1 * m2) * m3);
}

void Interpolate( glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result ) {
  int N = result.size();
  glm::vec2 step = glm::vec2(b-a) / float(max(N-1,1));
  glm::vec2 current( a );

  for( int i=0; i<N; ++i ) {
    result[i] = current;
    current += step;
  }
}

/*Place updates of parameters here*/
bool Update() {
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
      return false;
    } else if (e.type == SDL_KEYDOWN) {
      int key_code = e.key.keysym.sym;
      switch(key_code) {
	      case SDLK_UP:
          /* Move camera forward */
          break;
	      case SDLK_DOWN:
		      /* Move camera backwards */
		        break;
	      case SDLK_LEFT:
		      /* Move camera left */
		        break;
	      case SDLK_RIGHT:
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
