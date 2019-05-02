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
using glm::mat2;
using glm::vec3;
using glm::mat3;
using glm::distance;
using std::vector;
using glm::ivec2;
using glm::vec2;

SDL_Event event;

#define SCREEN_WIDTH 600
#define SCREEN_HEIGHT 600
#define FULLSCREEN_MODE false

struct Intersection {
  ivec2 position;
  ivec2 direction;
  float distance;
  int index;
};

struct Line {
  ivec2 start;
  ivec2 end;
  vec3 colour;
};

struct Scene {
  vector<Line> geometry;
};

struct Camera {
  ivec2 position;
  vec2 direction;
  float FOV;
};

Scene scene;

Camera camera;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void DrawScene(screen* screen, Scene scene);
void Draw3DView(screen* screen, vector<Intersection> intersections);

void DrawRect(screen* screen, ivec2 start, ivec2 end, vec3 colour);
void DrawLine(screen* screen, Line line);

void getTestIntersections(vector<Intersection> &intersections);

void InterpolatePixels(ivec2 a, ivec2 b, vector<ivec2>& result);
void Interpolate(float a, float b, vector<float>& result);

bool closestIntersection(vec2 start, vec2 dir, Intersection& closestIntersection);

float crossProd2D(vec2 v, vec2 w);

void GetRaysToCast(float FOV, vec2 dir, vector<vec2>& rays);
void LoadTestScene(Scene& scene);
void FillBlack(screen* screen);
void GetRotationMatrix(float theta, mat2& rotation);

int main(int argc, char* argv[]) {
  screen *mainscreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  screen *scenescreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  camera.position = ivec2(SCREEN_WIDTH/2, SCREEN_HEIGHT/2);
  camera.direction = vec2(1.0,1.0);
  camera.FOV = M_PI / 4;

  LoadTestScene(scene);

  while (Update()) {
    Draw(mainscreen);
    SDL_Renderframe(mainscreen);
    SDL_SaveImage(mainscreen, "mainout.bmp");

    DrawScene(scenescreen, scene);
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

  vec2 position = camera.position;
  vec2 direction = camera.direction;
  vector<vec2> rays;
  GetRaysToCast(camera.FOV, direction, rays);

  vector

  for (size_t i = 0; i < rays.size(); i++) {

  }

}

void DrawScene(screen* screen, Scene scene) {

  FillBlack(screen);

  for (size_t i = 0; i < scene.geometry.size(); i++) {
    DrawLine(screen, scene.geometry[i]);
  }

  /* Cast Ray Example */

  vec2 position = camera.position;
  vec2 direction = camera.direction;
  vector<vec2> rays;
  GetRaysToCast(camera.FOV, direction, rays);

  for (size_t i = 0; i < rays.size(); i++) {
    Line line;
    line.start = position;
    line.end = position + (50.0f * rays[i]);
    line.colour = vec3(1.0, 0.0, 0.0);
    DrawLine(screen, line);
  }

  /* Cast ray example end */

}


void Draw3DView(screen* screen, vector<Intersection> intersections) {
  if (intersections.size() < 1) return;

  int length = intersections.size();
  int width = floor(SCREEN_WIDTH / length);

  for (int i = 0; i < length; i++) {
    Intersection intersection = intersections[i];

    int height = floor(1 / intersection.distance);
    vec3 colour = scene.geometry[intersection.index].colour / (intersection.distance);

    int vGap = (SCREEN_HEIGHT - height) / 2;

    ivec2 topLeft(i * width, vGap);
    ivec2 botRight((i + 1) * width - 1, SCREEN_HEIGHT - vGap);

    DrawRect(screen, topLeft, botRight, colour);
  }
}

void DrawRect(screen* screen, ivec2 start, ivec2 end, vec3 colour) {
  int xDiff = end.x - start.x;
  int yDiff = end.y - start.y;

  for (int x = start.x; x <= xDiff; x++) {
    for (int y = start.y; y <= yDiff; y++) {
      PutPixelSDL(screen, x, y, colour);
    }
  }
}

void DrawLine(screen* screen, Line line) {
  int deltaX = abs(line.end.x - line.start.x);
  int deltaY = abs(line.end.y - line.start.y);

  int lineLength = (deltaX > deltaY) ? deltaX : deltaY;

  vector<ivec2> pixels(lineLength + 1);

  InterpolatePixels(line.start, line.end, pixels);

  for (int i = 0; i < pixels.size(); i++) {
    PutPixelSDL(screen, pixels[i].x, pixels[i].y, line.colour);
  }
}

void FillBlack(screen* screen) {
  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
      vec3 black = vec3(0.0f, 0.0f, 0.0f);
      PutPixelSDL(screen, x, y, black);
    }
  }
}

void GetRotationMatrix(float theta, mat2& rotation) {
  rotation[0][0] = cos(theta);
  rotation[0][1] = sin(theta);
  rotation[1][0] = -sin(theta);
  rotation[1][1] = cos(theta);
}

void GetRaysToCast(float FOV, vec2 dir, vector<vec2>& rays) {
  float step = FOV / (float) SCREEN_WIDTH;
  mat2 rotationMatrix;

  for (int i = 0; i < SCREEN_WIDTH; i++) {
    float theta = step * i;
    GetRotationMatrix(theta, rotationMatrix);
    vec2 ray = dir * rotationMatrix;
    rays.push_back(ray);
  }
}

void InterpolatePixels(ivec2 a, ivec2 b, vector<ivec2>& result) {
  int size = result.size();
  vector<float> x(size), y(size);

  Interpolate(a.x, b.x, x); Interpolate(a.y, b.y, y);

  /* Populate results vector */
  for (int i = 0; i < size; i++) {
    result[i].x = round(x[i]); result[i].y = round(y[i]);
  }
}

void Interpolate(float a, float b, vector<float>& result) {
  int size = result.size();
  float dist = b - a, incr = dist / float(max(1, (size - 1)));

  for (int i = 0; i < size; i++) result[i] = a + (incr * i);
}

bool closestIntersection(vec2 start, vec2 dir, Intersection& closestIntersection) {
  bool intersectionFound = false;
  closestIntersection = Intersection();
  closestIntersection.distance = std::numeric_limits<float>::max();

  for (int i = 0; i < scene.geometry.size(); i++) {
    Line line = scene.geometry[i];

    vec2 s = line.end - line.start;

    vec2 q = line.start;

    vec2 p = start;
    vec2 r = dir;

    float csr = crossProd2D(r,s);

    float t = crossProd2D((q - p), s) / csr;
    float u = crossProd2D(q - p, r) / csr;

    if (csr != 0 && t >= 0 && u >= 0 && u <= 1) {
      //the two line segments meet at the point p + t r = q + u s
      intersectionFound = true;

      vec2 intersectionPoint = p + t * r;
      float dist = distance(start, intersectionPoint);

      if (dist < closestIntersection.distance) {
        closestIntersection.position = intersectionPoint;
        closestIntersection.distance = dist;
        closestIntersection.index = i;
      }
    }

    if (csr == 0 && crossProd2D(q - p, r) != 0) {
      // parallel and non-intersecting
      break;
    }

  }

  return intersectionFound;
}

float crossProd2D(vec2 v, vec2 w) {
  return v.x * w.y - v.y * w.x;
}

void rotateCamera(Camera& camera, int direction, float delta) {
  mat2 rotation;

  GetRotationMatrix((float) direction * delta, rotation);

  camera.direction = camera.direction * rotation;

}

/*Place updates of parameters here*/
bool Update() {
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  //std::cout << "Render time: " << dt << " ms." << std::endl;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
	    return false;
	  } else if (e.type == SDL_KEYDOWN) {
	    int key_code = e.key.keysym.sym;
	    switch(key_code) {
        case SDLK_w:      camera.position.y -= SCREEN_HEIGHT / 100; break;
        case SDLK_s:      camera.position.y += SCREEN_HEIGHT / 100; break;
        case SDLK_a:      camera.position.x -= SCREEN_HEIGHT / 100; break;
        case SDLK_d:      camera.position.x += SCREEN_HEIGHT / 100; break;
        case SDLK_LEFT:   rotateCamera(camera, 1, (float) M_PI / 32); break;
        case SDLK_RIGHT:  rotateCamera(camera, -1, (float) M_PI / 32); break;
      }
    }
  }
  return true;
}

void LoadTestScene(Scene& scene) {

  struct FLine {
    vec2 start;
    vec2 end;
    vec3 colour;
  };

  vector<FLine> lines;

  vec3 black  = vec3(0.00f, 0.00f, 0.00f);
  vec3 red    = vec3(0.75f, 0.15f, 0.15f);
    vec3 yellow = vec3(0.75f, 0.75f, 0.15f);
    vec3 green  = vec3(0.15f, 0.75f, 0.15f);
    vec3 cyan   = vec3(0.15f, 0.75f, 0.75f);
    vec3 blue   = vec3(0.15f, 0.15f, 0.75f);
    vec3 purple = vec3(0.75f, 0.15f, 0.75f);
    vec3 white  = vec3(0.75f, 0.75f, 0.75f);

  /* ---------- WALLS ---------- */

  FLine leftWall, rightWall, topWall, bottomWall;

  leftWall.start   = vec2(0.0f, 0.0f); leftWall.end   = vec2(0.0f, 1.0f);
  rightWall.start  = vec2(1.0f, 0.0f); rightWall.end  = vec2(1.0f, 1.0f);
  topWall.start    = vec2(0.0f, 0.0f); topWall.end    = vec2(1.0f, 0.0f);
  bottomWall.start = vec2(0.0f, 1.0f); bottomWall.end = vec2(1.0f, 1.0f);

  leftWall.colour = red;
  rightWall.colour = yellow;
  topWall.colour = green;
  bottomWall.colour = cyan;

  lines.push_back(leftWall);
  lines.push_back(rightWall);
  lines.push_back(topWall);
  lines.push_back(bottomWall);

  /* ----------  BOX  ---------- */

  FLine A, B, C, D, E, F, G, H;

  A.start = vec2(5.0f/6.0f, 1.0f/6.0f); A.end = vec2(5.0f/6.0f, 2.0f/6.0f);
  B.start = vec2(5.0f/6.0f, 2.0f/6.0f); B.end = vec2(2.0f/6.0f, 2.0f/6.0f);
  C.start = vec2(2.0f/6.0f, 2.0f/6.0f); C.end = vec2(2.0f/6.0f, 4.0f/6.0f);
  D.start = vec2(2.0f/6.0f, 4.0f/6.0f); D.end = vec2(5.0f/6.0f, 4.0f/6.0f);
  E.start = vec2(5.0f/6.0f, 4.0f/6.0f); E.end = vec2(5.0f/6.0f, 5.0f/6.0f);
  F.start = vec2(5.0f/6.0f, 5.0f/6.0f); F.end = vec2(1.0f/6.0f, 5.0f/6.0f);
  G.start = vec2(1.0f/6.0f, 5.0f/6.0f); G.end = vec2(1.0f/6.0f, 1.0f/6.0f);
  H.start = vec2(1.0f/6.0f, 1.0f/6.0f); H.end = vec2(5.0f/6.0f, 1.0f/6.0f);

  A.colour = white;
  B.colour = white;
  C.colour = white;
  D.colour = white;
  E.colour = white;
  F.colour = white;
  G.colour = white;
  H.colour = white;

  lines.push_back(A);
  lines.push_back(B);
  lines.push_back(C);
  lines.push_back(D);
  lines.push_back(E);
  lines.push_back(F);
  lines.push_back(G);
  lines.push_back(H);

  /* Scale up scene */
  for (size_t i = 0; i < lines.size(); i++) {
    Line line;
    line.start = ivec2(lines[i].start.x * (SCREEN_WIDTH - 1), lines[i].start.y * (SCREEN_HEIGHT - 1));
    line.end = ivec2(lines[i].end.x * (SCREEN_WIDTH - 1), lines[i].end.y * (SCREEN_HEIGHT - 1));
    line.colour = lines[i].colour;
    scene.geometry.push_back(line);
  }
}
