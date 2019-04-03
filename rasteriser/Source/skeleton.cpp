#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include "glm/gtx/string_cast.hpp" // std::cout<<glm::to_string(hello)<<std::endl;

using namespace std;
using glm::ivec2;
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::mat3;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 600
#define SCREEN_HEIGHT 600
#define FULLSCREEN_MODE false

std::vector<Triangle> triangles;

const float focalLength = SCREEN_HEIGHT;
const vec4 defaultCameraPos(0.0, 0.0, -3.001, 1.0);

vec4 cameraPos(0, 0, -3.001, 1.0);
mat3 rotation;
mat4 transformMat;
mat4 projection;
float yaw = 0;
float pitch = 0;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, ivec2& p );
bool isWithinBounds(ivec2 v);
void TransformationMatrix(vec4 camPos, mat3 rot, mat4 &T);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void DrawLineSDL(screen* surface, ivec2 a, ivec2 b, vec3 colour);
void DrawPolygonEdges(screen* screen, const vector<vec4>& vertices, vec3 colour);
void getProjectionMatrix(mat4 &mat);
void clipTriangle(Triangle triangle, vector<vec4>& clippedTriangles);
bool isInside(float i, float w, float maxVal);

int main( int argc, char* argv[] ) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  // LoadTestModel(triangles);
  LoadTestTriangle(triangles);

  while ( Update() ) {
    Draw(screen);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

// void VertexShader(const vec4& v, ivec2& p) {
//   p.x = v.x + SCREEN_WIDTH / 2;
//   p.y = v.y + SCREEN_HEIGHT / 2;
// }

void VertexShader(const vec4& v, ivec2& p) {
  p.x = focalLength * (v.x / v.z) + SCREEN_WIDTH / 2;
  p.y = focalLength * (v.y / v.z) + SCREEN_HEIGHT / 2;
}

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  getRotationMatrix(pitch, yaw, 0, rotation);
  TransformationMatrix(cameraPos, rotation, transformMat);
  getProjectionMatrix(projection);

  for( uint32_t i=0; i<triangles.size(); ++i ) {
    vector<vec4> vertices(4);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

    vec3 colour = vec3(1, 1, 1);

    vector<vec4> clippedVertices;

    clipTriangle(triangles[i], clippedVertices);

    // TODO: Convert polygon to triangles and draw

    DrawPolygonEdges(screen, clippedVertices, colour);
  }
}

void clipTriangle(Triangle triangle, vector<vec4>& inVertices) {

  inVertices.clear();

  vector<vec4> vertices(3);
  vertices[0] = triangle.v0;
  vertices[1] = triangle.v1;
  vertices[2] = triangle.v2;

  for (int i = 0; i < 3; i++) {
    /* Perform transform on co-ordinate */
    vec4 transformedCoord = transformMat * vertices[i];
    /* Project co-ordinate to Homogenous space */
    vec4 homogenousCoord = projection * transformedCoord;
    std::cout << glm::to_string(homogenousCoord) << std::endl;
    /* Perform homogenous divide (projects to plane at w = 1) */
    vec4 homogenousDivide = (1/homogenousCoord.w) * homogenousCoord; // [u, v, f, 1]

    bool isInX = isInside(homogenousCoord.x, homogenousCoord.w, SCREEN_WIDTH/2);
    bool isInY = isInside(homogenousCoord.y, homogenousCoord.w, SCREEN_WIDTH/2);
    std::cout << "isInX: " << isInX << ", isInY: " << isInY << std::endl;

    inVertices.push_back(vertices[i]);
  }
}

bool isInside(float i, float w, float maxVal) {
  return (i <= (w * maxVal)) && (i >= (w * -maxVal));
}

bool isWithinBounds(ivec2 v) {
  return v.x > 0 && v.x < SCREEN_WIDTH && v.y > 0 && v.y < SCREEN_HEIGHT;
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result ) {
  int N = result.size();
  vec2 step = glm::vec2(b - a) / float(max(N - 1, 1));
  vec2 current(a);

  for (int i = 0; i < N; ++i) {
    result[i] = current;
    current += step;
  }
}

void DrawLineSDL(screen* screen, ivec2 a, ivec2 b, vec3 colour) {
  ivec2 delta = glm::abs(a - b);
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<ivec2> line(pixels);
  Interpolate(a, b, line);

  for (int px = 0; px < pixels; px++) {
    ivec2 pixel = line[px];
    if (isWithinBounds(pixel)) PutPixelSDL(screen, pixel.x, pixel.y, colour);
  }
}

void DrawPolygonEdges(screen* screen, const vector<vec4>& vertices, vec3 colour) {
  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D image position:
  vector<ivec2> projectedVertices(V);
  for (int i = 0; i < V; ++i) {
    // /* Perform transform on co-ordinate */
    // vec4 transformedCoord = transformMat * vertices[i];
    // /* Project co-ordinate to Homogenous space */
    // vec4 homogenousCoord = projection * transformedCoord;
    // std::cout << glm::to_string(homogenousCoord) << std::endl;
    // /* Perform homogenous divide (projects to plane at w = 1) */
    // vec4 homogenousDivide = (1/homogenousCoord.w) * homogenousCoord;
    /* Get position on screen */
    VertexShader(transformMat * vertices[i], projectedVertices[i]);
  }
  // Loop over all vertices and draw the edge from it to the next vertex:
  for (int i = 0; i < V; ++i) {
    int j = (i + 1) % V; // The next vertex
    DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], colour);
  }
}

void getProjectionMatrix(mat4 &mat) {
  mat = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 1/focalLength), vec4(0, 0, 0, 0));
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

void TransformationMatrix(vec4 camPos, mat3 rot, mat4 &T) {
  mat4 m1 = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), camPos);
  mat4 m2 = mat4(rot);
  mat4 m3 = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), -camPos);
  m1[3][3] = 1; m2[3][3] = 1; m3[3][3] = 1;
  // T = m1 * m2 * m3;
  T = m2 * m3;
}

void moveCameraRight(int direction, float distance) {
	vec4 right(rotation[0][0], rotation[0][1], rotation[0][2], 0);
	cameraPos += direction * distance * right;
}

void moveCameraUp(int direction, float distance) {
	vec4 up(rotation[1][0], rotation[1][1], rotation[1][2], 0);
	cameraPos += direction * distance * up;
}

void moveCameraForward(int direction, float distance) {
	vec4 forward(rotation[2][0], rotation[2][1], rotation[2][2], 0);
	cameraPos += direction * distance * forward;
}

/* Place updates of parameters here */
bool Update() {
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  std::cout << "Render time: " << dt << " ms. --------------------------------------------" << std::endl;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
  if (e.type == SDL_QUIT) {
    return false;
  } else if (e.type == SDL_KEYDOWN) {
    int key_code = e.key.keysym.sym;
    switch(key_code) {
      case SDLK_UP:
        pitch += M_PI / 18;
        break;
      case SDLK_DOWN:
        pitch -= M_PI / 18;
        break;
      case SDLK_LEFT:
        yaw -= M_PI / 18;
        break;
      case SDLK_RIGHT:
        yaw += M_PI / 18;
        break;
      case SDLK_r:
        /* Look-At function, points camera to 0,0,0 */
        // lookAt(ctw);
        break;
      case SDLK_t:
        // Reset camera position
        cameraPos = defaultCameraPos;
        pitch = 0;
        yaw = 0;
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
        return false;
      }
    }
  }

  return true;
}
