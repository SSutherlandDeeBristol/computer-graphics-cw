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
using glm::ivec2;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

const float focalLength = SCREEN_HEIGHT;

std::vector<Triangle> triangles;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, ivec2& p );
bool isWithinBounds(ivec2 v);
glm::mat4x4 TransformationMatrix(glm::vec4 camPos, glm::mat3x3 rot);

int main( int argc, char* argv[] ) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  LoadTestModel(triangles);

  while ( Update() ) {
    Draw(screen);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

void VertexShader( const vec4& v, ivec2& p ) {
  int x = 0, y = 0;

  x = focalLength * (v.x / v.z) + SCREEN_WIDTH / 2;
  y = focalLength * (v.y / v.z) + SCREEN_HEIGHT / 2;

  p.x = x;
  p.y = y;
}

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  for( uint32_t i=0; i<triangles.size(); ++i ) {
    vector<vec4> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;
    for(int v=0; v<3; ++v) {
      ivec2 projPos;
      VertexShader( vertices[v], projPos );
      vec3 color(1,1,1);
      if (isWithinBounds(projPos)) PutPixelSDL( screen, projPos.x, projPos.y, color );
    }
  }
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

bool isWithinBounds(ivec2 v) {
  return v.x > 0 && v.x < SCREEN_WIDTH && v.y > 0 && v.y < SCREEN_HEIGHT;
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
