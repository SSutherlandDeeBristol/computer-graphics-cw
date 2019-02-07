#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;
const float velocity = 0.001f;
vector<vec3> stars(1000);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen);
void DrawRainbow(screen* screen);
void Interpolate(float a, float b, vector<float>& result);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);

int main(int argc, char* argv[]) {
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  for (size_t i = 0; i < stars.size(); i++) {
    vec3 star(0.0,0.0,0.0);
    star.x = float(rand()) / float(RAND_MAX) * (float(rand()) / float(RAND_MAX) > 0.5 ? 1 : -1);
    star.y = float(rand()) / float(RAND_MAX) * (float(rand()) / float(RAND_MAX) > 0.5 ? 1 : -1);
    star.z = float(rand()) / float(RAND_MAX);
    stars[i] = star;
  }

  while (NoQuitMessageSDL()){
    Draw(screen);
    Update();
    SDL_Renderframe(screen);
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));


  int projectedX, projectedY;
  float focalDistance = SCREEN_HEIGHT / 2;

  for (size_t s = 0; s < stars.size(); s++) {
    vec3 star = stars[s];
    vec3 colour = 0.2f * vec3(1,1,1) / (star.z * star.z);

    projectedX = (focalDistance * (star.x / star.z)) + (SCREEN_WIDTH / 2);
    projectedY = (focalDistance * (star.y / star.z)) + (SCREEN_HEIGHT / 2);

    if (projectedX > 0 && projectedX < SCREEN_WIDTH && projectedY > 0 && projectedY < SCREEN_HEIGHT)
      PutPixelSDL(screen, projectedX, projectedY, colour);
  }

}

void DrawRainbow(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 topLeft(1,0,0); // red
  vec3 topRight(0,0,1); // blue
  vec3 bottomRight(0,1,0); // green
  vec3 bottomLeft(1,1,0); // yellow

  vector<vec3> leftSide(SCREEN_HEIGHT);
  vector<vec3> rightSide(SCREEN_HEIGHT);
  Interpolate (topLeft, bottomLeft, leftSide);
  Interpolate (topRight, bottomRight, rightSide);

  for (int i = 0; i < SCREEN_HEIGHT; i++) {
    vector<vec3> row( SCREEN_WIDTH );
    Interpolate(leftSide[i], rightSide[i], row);

    for (std::vector<vec3>::size_type j = 0; j < row.size(); j++) {
      PutPixelSDL(screen, j, i, row[j]);
    }
  }
}

/*Place updates of parameters here*/
void Update() {
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  // std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/

  for (size_t s = 0; s < stars.size(); ++s) {
    // Add code for update of stars
    stars[s].z -= velocity * dt;
    if (stars[s].z <= 0) stars[s].z += 1;
    if (stars[s].z > 1) stars[s].z -= 1;
  }

}

void Interpolate(float a, float b, vector<float>& result) {
  int size = result.size();
  float dist = b - a, incr = dist / (float) (size - 1);

  if (size <= 0) {
  } else if (size == 1) {
    result[0] = a;
  } else if (size == 2) {
    result[0] = a; result[1] = b;
  } else {
    result[0] = a; result[size - 1] = b;
    for (int i = 1; i < size - 1; i++) result[i] = a + (incr * i);
  }
}

void Interpolate(vec3 a, vec3 b, vector<vec3>& result) {
  int size = result.size();
  vector<float> x(size), y(size), z(size);

  Interpolate(a.x, b.x, x); Interpolate(a.y, b.y, y); Interpolate(a.z, b.z, z);

  for (int i = 0; i < size; i++) {
    result[i].x = x[i]; result[i].y = y[i]; result[i].z = z[i];
  }
}