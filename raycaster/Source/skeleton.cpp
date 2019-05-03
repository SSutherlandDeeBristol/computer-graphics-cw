#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include <stdint.h>
#include "limits"
#include <math.h>
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

#define FULLSCREEN_MODE false

int SCREEN_WIDTH = 600;
int SCREEN_HEIGHT = 600;

const vec3 SKY_COLOUR = vec3(0.78f, 0.88f, 0.91f);
const vec3 FLOOR_COLOUR = vec3(0.10f, 0.10f, 0.10f);
const vec3 CONE_COLOUR = vec3(1.00f, 0.60f, 0.10f);
const vec3 CAMERA_COLOUR = vec3(1.00f, 1.00f, 1.00f);

const vec2 DEFAULT_DIRECTION = vec2(1.0f, 1.0f);

struct Intersection {
  vec2 position;
  vec2 direction;
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
void Draw3DView(screen* screen, Scene scene, vector<vec2> rays, vector<Intersection>& intersections);
void DrawMinimap(screen* screen, Scene scene, vector<vec2> rays, vector<Intersection> intersections, bool floodCone);

void GetRaysToCast(float FOV, vec2 dir, vector<vec2>& rays);
void LoadTestScene(Scene& scene);
void DrawFloorAndSky(screen* screen);
void TranslatePositionToMinimap(ivec2 topLeft, ivec2 botRight, ivec2 position, ivec2& result);
void TranslateLineToMinimap(ivec2 topLeft, ivec2 botRight, Line line, Line& result);

void FillBlack(screen* screen);
void DrawLine(screen* screen, Line line);
void DrawRect(screen* screen, ivec2 start, ivec2 end, vec3 colour);

bool ClosestIntersection(Scene scene, vec2 start, vec2 dir, Intersection& closestIntersection);
float CrossProduct(vec2 v, vec2 w);
void InterpolatePixels(ivec2 a, ivec2 b, vector<ivec2>& result);
void Interpolate(float a, float b, vector<float>& result);
bool IsWithinScreenBounds(ivec2 p);
bool IsWithinScreenBounds(int x, int y);

float GetFarPlane();
void GetRotationMatrix(float theta, mat2& rotation);
void RotateCamera(Camera& camera, int direction, float delta);
void MoveCameraForward(Camera& camera, int distance);
void MoveCameraBackward(Camera& camera, int distance);
void MoveCameraLeft(Camera& camera, int distance);
void MoveCameraRight(Camera& camera, int distance);
void ResetCamera(Camera& camera);
void LookAt(Camera& camera, ivec2 position);

int main(int argc, char* argv[]) {

  /* Parse arguments */
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " SCREEN_WIDTH SCREEN_HEIGHT" << endl;
    return 1;
  } else {
    try {
      SCREEN_WIDTH = stoi(argv[1]);
      SCREEN_HEIGHT = stoi(argv[2]);
    } catch (exception const &e) {
      cerr << "Could not parse arguments." << endl;
    }
  }

  screen *mainscreen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);

  ResetCamera(camera);
  LoadTestScene(scene);

  while (Update()) {
    Draw(mainscreen);
    SDL_Renderframe(mainscreen);
    SDL_SaveImage(mainscreen, "mainout.bmp");
  }

  KillSDL(mainscreen);

	return 0;
}

void Draw(screen* screen) {

  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  /* Draw the background */
  DrawFloorAndSky(screen);

  /* Gets and draw the ray intersections */
  vector<vec2> rays;
  vector<Intersection> intersections;
  GetRaysToCast(camera.FOV, camera.direction, rays);
  Draw3DView(screen, scene, rays, intersections);
  DrawMinimap(screen, scene, rays, intersections, true);
}

void DrawMinimap(screen* screen, Scene scene, vector<vec2> rays, vector<Intersection> intersections, bool floodCone) {
  ivec2 topLeft = ivec2(SCREEN_WIDTH * 0.68f, SCREEN_HEIGHT * 0.68f);
  ivec2 botRight = ivec2(SCREEN_WIDTH * 0.98f, SCREEN_HEIGHT * 0.98f);

  /* Draw minimap background */
  DrawRect(screen, topLeft, botRight, vec3(0.0f, 0.0f, 0.0f));

  /* Draw scene geometry */
  for (size_t i = 0; i < scene.geometry.size(); i++) {
    Line scaledLine;
    TranslateLineToMinimap(topLeft, botRight, scene.geometry[i], scaledLine);
    DrawLine(screen, scaledLine);
  }

  /* Draw viewing cone from rays */
  if (rays.size() == intersections.size()) {
    for (size_t i = 0; i < rays.size(); i++) {
      if (intersections[i].distance < GetFarPlane()) {
        Line line, scaledLine;
        line.start = camera.position;
        line.end = floodCone ? intersections[i].position : (vec2) camera.position + (20.0f * rays[i]);
        line.colour = CONE_COLOUR;

        TranslateLineToMinimap(topLeft, botRight, line, scaledLine);
        DrawLine(screen, scaledLine);
      }
    }
  }

  /* Draw player position */
  ivec2 cameraTopLeft, cameraBotRight;
  TranslatePositionToMinimap(topLeft, botRight, camera.position - ivec2(2, 2), cameraTopLeft);
  TranslatePositionToMinimap(topLeft, botRight, camera.position + ivec2(2, 2), cameraBotRight);

  DrawRect(screen, cameraTopLeft, cameraBotRight, CAMERA_COLOUR);
}

void TranslateLineToMinimap(ivec2 topLeft, ivec2 botRight, Line line, Line& result) {
  TranslatePositionToMinimap(topLeft, botRight, line.start, result.start);
  TranslatePositionToMinimap(topLeft, botRight, line.end, result.end);
  result.colour = line.colour;
}

void TranslatePositionToMinimap(ivec2 topLeft, ivec2 botRight, ivec2 position, ivec2& result) {
  int width = botRight.x - topLeft.x, height = botRight.y - topLeft.y;

  float sfx = (float) width / (float) SCREEN_WIDTH;
  float sfy = (float) height / (float) SCREEN_HEIGHT;

  result.x = (int) (position.x * sfx);
  result.y = (int) (position.y * sfy);
  result += topLeft;
}

void DrawRect(screen* screen, ivec2 start, ivec2 end, vec3 colour) {
  if (start.x > SCREEN_WIDTH) start.x = SCREEN_WIDTH;
  if (start.x < 0) start.x = 0;
  if (start.y > SCREEN_HEIGHT) start.y = SCREEN_HEIGHT;
  if (start.y < 0) start.y = 0;

  for (int x = start.x; x <= end.x; x++) {
    for (int y = start.y; y <= end.y; y++) {
      if (IsWithinScreenBounds(x, y)) {
        PutPixelSDL(screen, x, y, colour);
      }
    }
  }
}

void Draw3DView(screen* screen, Scene scene, vector<vec2> rays, vector<Intersection>& intersections) {
  if ((int) rays.size() != SCREEN_WIDTH) return;
  intersections.clear();

  for (size_t i = 0; i < rays.size(); i++) {
    Intersection intersection; ClosestIntersection(scene, camera.position, rays[i], intersection);

    /* Multiply by cosine of angle between ray and direction to counteract fisheye effect */
    intersection.distance *= (dot(rays[i], camera.direction) / (length(rays[i]) * length(camera.direction)));
    intersections.push_back(intersection);
  }

  for (int i = 0; i < SCREEN_WIDTH; i++) {
    Intersection intersection = intersections[i];

    /* Clip to far plane */
    if (intersection.distance < GetFarPlane()) {
      int height = min(SCREEN_HEIGHT, (int) floor(SCREEN_HEIGHT * 50 / intersection.distance));
      vec3 colour = 10.0f * scene.geometry[intersection.index].colour / sqrt(intersection.distance);

      int vGap = max(0, (SCREEN_HEIGHT - height) / 2);

      Line line;
      line.start = ivec2(i, vGap);
      line.end = ivec2(i, vGap + height);
      line.colour = colour;

      DrawLine(screen, line);
    }
  }
}

float GetFarPlane() {
  return SCREEN_HEIGHT * 10;
}

/* Determines whether a pixel is within the screen bounds */
bool IsWithinScreenBounds(ivec2 p) {
  return IsWithinScreenBounds(p.x, p.y);
}

/* Determines whether a pixel is within the screen bounds */
bool IsWithinScreenBounds(int x, int y) {
  return x > 0 && x < SCREEN_WIDTH && y > 0 && y < SCREEN_HEIGHT;
}

void DrawLine(screen* screen, Line line) {
  int deltaX = abs(line.end.x - line.start.x);
  int deltaY = abs(line.end.y - line.start.y);

  int lineLength = (deltaX > deltaY) ? deltaX : deltaY;

  vector<ivec2> pixels(lineLength + 1);

  InterpolatePixels(line.start, line.end, pixels);

  for (size_t i = 0; i < pixels.size(); i++) {
    if (IsWithinScreenBounds(pixels[i])) PutPixelSDL(screen, pixels[i].x, pixels[i].y, line.colour);
  }
}

void FillBlack(screen* screen) {
  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
      PutPixelSDL(screen, x, y, vec3(0.0f, 0.0f, 0.0f));
    }
  }
}

void DrawFloorAndSky(screen* screen) {
  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
      PutPixelSDL(screen, x, y, y < SCREEN_HEIGHT / 2 ? SKY_COLOUR : FLOOR_COLOUR);
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

  float theta = step * SCREEN_WIDTH / 2;
  GetRotationMatrix(theta, rotationMatrix);
  vec2 startDir = dir * rotationMatrix;

  /* Populate rays vector */
  for (int i = 0; i < SCREEN_WIDTH; i++) {
    float theta = -step * i;
    GetRotationMatrix(theta, rotationMatrix);
    vec2 ray = startDir * rotationMatrix;
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

bool ClosestIntersection(Scene scene, vec2 start, vec2 dir, Intersection& closestIntersection) {
  bool intersectionFound = false;
  closestIntersection = Intersection();
  closestIntersection.distance = std::numeric_limits<float>::max();

  for (size_t i = 0; i < scene.geometry.size(); i++) {
    Line line = scene.geometry[i];

    vec2 s = line.end - line.start;
    vec2 q = line.start;
    vec2 p = start, r = dir;

    float csr = CrossProduct(r, s);
    float t = CrossProduct((q - p), s) / csr;
    float u = CrossProduct(q - p, r) / csr;

    if (csr != 0 && t >= 0 && u >= 0 && u <= 1) {
      /* The two line segments meet at the point p + t r = q + u s */
      intersectionFound = true;

      vec2 intersectionPoint = p + t * r;
      float dist = distance(start, intersectionPoint);

      if (dist < closestIntersection.distance) {
        closestIntersection.position = intersectionPoint;
        closestIntersection.distance = dist;
        closestIntersection.index = i;
      }
    }
  }

  return intersectionFound;
}

float CrossProduct(vec2 v, vec2 w) {
  return v.x * w.y - v.y * w.x;
}

void RotateCamera(Camera& camera, int direction, float delta) {
  mat2 rotation; GetRotationMatrix((float) direction * delta, rotation);
  camera.direction = normalize(camera.direction * rotation);
}

void MoveCameraForward(Camera& camera, int distance) {
  camera.position += (float) distance * camera.direction;
}

void MoveCameraBackward(Camera& camera, int distance) {
  camera.position -= (float) distance * camera.direction;
}

void MoveCameraLeft(Camera& camera, int distance) {
  camera.position += (float) distance * vec2(camera.direction.y, -camera.direction.x);
}

void MoveCameraRight(Camera& camera, int distance) {
  camera.position -= (float) distance * vec2(camera.direction.y, -camera.direction.x);
}

void ResetCamera(Camera& camera) {
  camera.position = ivec2(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);
  camera.direction = DEFAULT_DIRECTION;
  camera.FOV = M_PI / 4;
}

void LookAt(Camera& camera, ivec2 position) {
  vec2 newDirection = normalize((vec2) (ivec2(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2) - camera.position));
  camera.direction = length(newDirection) > 0.0f ? newDirection : DEFAULT_DIRECTION;
}

bool Update() {
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  std::cout << "Render time: " << dt << " ms." << std::endl;

  SDL_Event e;
  while (SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
	    return false;
	  } else if (e.type == SDL_KEYDOWN) {
	    int key_code = e.key.keysym.sym;
	    switch(key_code) {
        case SDLK_w:      MoveCameraForward(camera, SCREEN_HEIGHT / 100); break;
        case SDLK_s:      MoveCameraBackward(camera, SCREEN_HEIGHT / 100); break;
        case SDLK_a:      MoveCameraLeft(camera, SCREEN_HEIGHT / 100); break;
        case SDLK_d:      MoveCameraRight(camera, SCREEN_HEIGHT / 100); break;
        case SDLK_LEFT:   RotateCamera(camera, 1, (float) M_PI / 64); break;
        case SDLK_RIGHT:  RotateCamera(camera, -1, (float) M_PI / 64); break;
        case SDLK_t:      ResetCamera(camera); break;
        case SDLK_r:      LookAt(camera, ivec2(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2)); break;
        case SDLK_ESCAPE: return false;
      }
    }
  }
  return true;
}

void LoadTestScene(Scene& scene) {

  /* Float line struct */
  struct FLine {
    vec2 start;
    vec2 end;
    vec3 colour;
  };

  vector<FLine> lines;

  /* Define colours */
  vec3 red    = vec3(0.75f, 0.15f, 0.15f);
  vec3 yellow = vec3(0.75f, 0.75f, 0.15f);
  vec3 green  = vec3(0.15f, 0.75f, 0.15f);
  vec3 cyan   = vec3(0.15f, 0.75f, 0.75f);
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
    line.start = ivec2(lines[i].start.x * SCREEN_WIDTH, lines[i].start.y * SCREEN_HEIGHT);
    line.end = ivec2(lines[i].end.x * SCREEN_WIDTH, lines[i].end.y * SCREEN_HEIGHT);
    line.colour = lines[i].colour;
    scene.geometry.push_back(line);
  }
}
