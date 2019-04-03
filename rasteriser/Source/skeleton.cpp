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
#define CLIP_OFFSET 50
#define FULLSCREEN_MODE false

std::vector<Triangle> triangles;

const float focalLength = SCREEN_HEIGHT;
const vec4 defaultCameraPos(0.0, 0.0, -3.001, 1.0);

vec4 cameraPos(0, 0, -3.001, 1.0);
mat3 rotationMat;
mat4 transformMat;
mat4 projectionMat;
float yaw = 0;
float pitch = 0;

enum Axis {X = 0, Y = 1, Z = 2, W = 3};

struct Pixel {
  int x;
  int y;
};

struct Vertex {
  vec4 position;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);

bool isWithinScreenBounds(Pixel p);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result );

void clipHomogenousVerticesToPlane(vector<Vertex>& clipBuffer, Axis axis, float maxVal, bool pos);
void clipTriangle(Triangle triangle, vector<Vertex>& inVertices);
bool isInsidePos(Vertex vertex, Axis axis, float maxVal);
bool isInsideNeg(Vertex vertex, Axis axis, float maxVal);
void calcIntersection(Vertex a, Vertex b, Vertex& c, Axis axis, float maxVal);
void homogenousFlatten(Vertex homogenousVertex, Vertex& flatVertex);
void homogenousDivide(Vertex homogenousVertex, Vertex& projectedVertex);
void projectToHomogenous(Vertex vertex, mat4 proj, Vertex& projectedVertex);

void TransformationMatrix(vec4 camPos, mat3 rot, mat4 &T);
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void getProjectionMatrix(mat4 &mat);

void VertexShader( const Vertex& v, Pixel& p );
void DrawClipOffset(screen* screen, bool fillOutline);
void DrawLineSDL(screen* surface, Pixel a, Pixel b, vec3 colour);
void DrawPolygonEdges(screen* screen, const vector<Vertex>& vertices, vec3 colour);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 colour);
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec3 colour);

int main(int argc, char* argv[]) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  LoadTestModel(triangles);
  // LoadTestTriangle(triangles);

  while (Update()) {
    Draw(screen);
    DrawClipOffset(screen, false);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

void VertexShader(const Vertex& v, Pixel& p) {
  p.x = focalLength * (v.position.x / v.position.z) + SCREEN_WIDTH / 2;
  p.y = focalLength * (v.position.y / v.position.z) + SCREEN_HEIGHT / 2;
}

void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  getRotationMatrix(pitch, yaw, 0, rotationMat);
  TransformationMatrix(cameraPos, rotationMat, transformMat);
  getProjectionMatrix(projectionMat);

  for (size_t i = 0; i < triangles.size(); i++) {
    /* Transform the triangle to camera space */
    vec4 tv0 = transformMat * triangles[i].v0;
    vec4 tv1 = transformMat * triangles[i].v1;
    vec4 tv2 = transformMat * triangles[i].v2;
    Triangle transformedTriangle = Triangle(tv0, tv1, tv2, triangles[i].color);

    vector<Vertex> clippedVertices;

    clipTriangle(transformedTriangle, clippedVertices);

    // TODO: Convert polygons to triangles

    DrawPolygonEdges(screen, clippedVertices, triangles[i].color); // Draw outlines
    // DrawPolygon(screen, clippedVertices, triangles[i].color); // Draw filled
  }
}

void clipTriangle(Triangle triangle, vector<Vertex>& inVertices) {
  float maxX = SCREEN_WIDTH/2 - CLIP_OFFSET, maxY = SCREEN_HEIGHT/2 - CLIP_OFFSET;
  vector<Vertex> clipBuffer;
  inVertices.clear();

  vector<Vertex> vertices(3);
  vertices[0].position = triangle.v0; vertices[1].position = triangle.v1; vertices[2].position = triangle.v2;

  /* Project (already transformed) vertices into Homogenous space */
  for (size_t i = 0; i < vertices.size(); i++) {
    Vertex homogenousCoord = vertices[i];
    projectToHomogenous(vertices[i], projectionMat, homogenousCoord);
    clipBuffer.push_back(homogenousCoord);
  }

  /* Clip to each plane in turn, maintaining a buffer of clipped vertices */
  clipHomogenousVerticesToPlane(clipBuffer, X, maxX, true);  /* Right  */
  clipHomogenousVerticesToPlane(clipBuffer, X, maxX, false); /* Left   */
  clipHomogenousVerticesToPlane(clipBuffer, Y, maxY, true);  /* Top    */
  clipHomogenousVerticesToPlane(clipBuffer, Y, maxY, false); /* Bottom */

  /* Flatten all vertices in the buffer (set their w component to 1) */
  for (size_t i = 0; i < clipBuffer.size(); i++) {
    Vertex flattenedCoord = clipBuffer[i];
    homogenousFlatten(clipBuffer[i], flattenedCoord);
    inVertices.push_back(flattenedCoord);
  }
}

void clipHomogenousVerticesToPlane(vector<Vertex>& clipBuffer, Axis axis, float maxVal, bool pos) {
  vector<Vertex> homogenousVertices = clipBuffer;
  size_t vs = homogenousVertices.size();
  clipBuffer.clear();

  for (size_t i = 0; i < vs; i++) {
    Vertex curVertex = homogenousVertices[i];
    Vertex nxtVertex = homogenousVertices[(i + 1) % vs];

    /* Check whether current and next vertices require clipping */
    bool isCurIn = pos ? isInsidePos(curVertex, axis, maxVal) : isInsideNeg(curVertex, axis, maxVal);
    bool isNxtIn = pos ? isInsidePos(nxtVertex, axis, maxVal) : isInsideNeg(nxtVertex, axis, maxVal);

    /* If the current vertex does not require clipping, push it unchanged */
    if (isCurIn) clipBuffer.push_back(curVertex);

    /* If either the current or the next vertex require clipping, push a new vertex on the boundary intersection */
    if ((isCurIn && !isNxtIn) || (!isCurIn && isNxtIn)) {
      Vertex newVertex;
      calcIntersection(curVertex, nxtVertex, newVertex, axis, pos ? maxVal : -maxVal);
      clipBuffer.push_back(newVertex);
    }
  }
}

/* Flattens a homogenous vertex (sets w to 1) */
void homogenousFlatten(Vertex homogenousVertex, Vertex& flatVertex) {
  flatVertex = homogenousVertex;
  flatVertex.position.w = 1;
}

/* Performs a homogenous divide on a vertex, projecting to the 3D plane w = 1 */
void homogenousDivide(Vertex homogenousVertex, Vertex& projectedVertex) {
  projectedVertex = homogenousVertex;
  projectedVertex.position = (1/homogenousVertex.position.w) * homogenousVertex.position;
}

/* Projects a 3D point into homogenous space */
void projectToHomogenous(Vertex vertex, mat4 proj, Vertex& projectedVertex) {
  projectedVertex = vertex;
  projectedVertex.position = proj * vertex.position;
}

void calcIntersection(Vertex a, Vertex b, Vertex& c, Axis axis, float maxVal) {
  float t = (a.position[axis] - maxVal * a.position.w) /
            ((maxVal * (b.position.w - a.position.w)) - (b.position[axis] - a.position[axis]));
  c.position = a.position + t * (b.position - a.position);
}

bool isInsidePos(Vertex vertex, Axis axis, float maxVal) {
  return (vertex.position[axis] <= (vertex.position.w * maxVal));
}

bool isInsideNeg(Vertex vertex, Axis axis, float maxVal) {
  return (vertex.position[axis] >= (vertex.position.w * -maxVal));
}

bool isWithinScreenBounds(Pixel p) {
  return p.x > 0 && p.x < SCREEN_WIDTH && p.y > 0 && p.y < SCREEN_HEIGHT;
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result ) {
  int N = result.size();

  float stepX = (b.x - a.x) / float(max(N - 1, 1));
  float stepY = (b.y - a.y) / float(max(N - 1, 1));

  float currentX(a.x);
  float currentY(a.y);

  for (int i = 0; i < N; ++i) {
    result[i].x = round(currentX);
    result[i].y = round(currentY);
    currentX += stepX;
    currentY += stepY;
  }
}

void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec3 colour) {
  int V = vertices.size();
  vector<Pixel> vertexPixels(V);

  for (int i = 0; i < V; i++) VertexShader(vertices[i], vertexPixels[i]);

  vector<Pixel> leftPixels; vector<Pixel> rightPixels;

  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(screen, leftPixels, rightPixels, colour);
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
  float maxY = -numeric_limits<float>::max();
  float minY = numeric_limits<float>::max();

  if (vertexPixels.size() == 0) return;

  for (size_t i = 0; i < vertexPixels.size(); i++) {
    if (vertexPixels[i].y >= maxY) maxY = vertexPixels[i].y;
    if (vertexPixels[i].y <= minY) minY = vertexPixels[i].y;
  }

  int numRows = static_cast<int>(maxY - minY + 1);

  leftPixels.resize(numRows); rightPixels.resize(numRows);

  for (int i = 0; i < numRows; i++) {
    leftPixels[i].x = numeric_limits<int>::max();
    leftPixels[i].y = i + minY;
    rightPixels[i].x = -numeric_limits<int>::max();
    rightPixels[i].y = i + minY;
  }

  for (size_t i = 0; i < vertexPixels.size(); i++) {
    Pixel curVertex = vertexPixels[i];
    Pixel nxtVertex = vertexPixels[(i + 1) % vertexPixels.size()];

    int deltaX = abs(curVertex.x - nxtVertex.x);
    int deltaY = abs(curVertex.y - nxtVertex.y);

    int lineLength = (deltaX > deltaY) ? deltaX : deltaY;

    vector<Pixel> pixels(lineLength + 1);

    Interpolate(curVertex, nxtVertex, pixels);

    for (size_t j = 0; j < pixels.size(); j++) {
      Pixel curPixel = pixels[j];
      int yindex = curPixel.y - minY;

      if (curPixel.x <= leftPixels[yindex].x) leftPixels[yindex] = curPixel;
      if (curPixel.x >= rightPixels[yindex].x) rightPixels[yindex] = curPixel;
    }
  }
}

void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 colour) {
  if (leftPixels.size() != rightPixels.size()) {
    std::cout << "Error in DrawPolygonRows: Left =/= Right" << std::endl;
    return;
  }

  for (size_t i = 0; i < leftPixels.size(); i++) {
    vector<Pixel> pixels(rightPixels[i].x - leftPixels[i].x + 1);
    Interpolate(leftPixels[i], rightPixels[i], pixels);

    for (size_t j = 0; j < pixels.size(); j++) {
      if (isWithinScreenBounds(pixels[j])) PutPixelSDL(screen, pixels[j].x, pixels[j].y, colour);
    }
  }
}

void DrawClipOffset(screen* screen, bool fillOutline) {
  vec3 colour = vec3(1, 0, 0);

  Pixel TL, TR, BL, BR;

  TL.x = CLIP_OFFSET; TL.y = CLIP_OFFSET;
  TR.x = SCREEN_WIDTH - CLIP_OFFSET; TR.y = CLIP_OFFSET;
  BL.x = CLIP_OFFSET; BL.y = SCREEN_HEIGHT - CLIP_OFFSET;
  BR.x = SCREEN_WIDTH - CLIP_OFFSET; BR.y = SCREEN_HEIGHT - CLIP_OFFSET;

  if (fillOutline) {
    DrawLineSDL(screen, TL, TR, colour); DrawLineSDL(screen, TR, BR, colour);
    DrawLineSDL(screen, BR, BL, colour); DrawLineSDL(screen, BL, TL, colour);
  } else {
    PutPixelSDL(screen, TL.x, TL.y, colour); PutPixelSDL(screen, TR.x, TR.y, colour);
    PutPixelSDL(screen, BL.x, BL.y, colour); PutPixelSDL(screen, BR.x, BR.y, colour);
  }
}

void DrawLineSDL(screen* screen, Pixel a, Pixel b, vec3 colour) {
  ivec2 delta = ivec2(glm::abs(a.x - b.x), glm::abs(a.y - b.y));
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<Pixel> line(pixels);
  Interpolate(a, b, line);

  for (int px = 0; px < pixels; px++) {
    Pixel pixel = line[px];
    if (isWithinScreenBounds(pixel)) PutPixelSDL(screen, pixel.x, pixel.y, colour);
  }
}

void DrawPolygonEdges(screen* screen, const vector<Vertex>& vertices, vec3 colour) {
  int V = vertices.size();
  /* Transform each vertex from 3D world position to 2D image position */
  vector<Pixel> projectedVertices(V);
  for (int i = 0; i < V; ++i) VertexShader(vertices[i], projectedVertices[i]);

  /* Loop over all vertices and draw the edge from it to the next vertex */
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
	vec4 right(rotationMat[0][0], rotationMat[0][1], rotationMat[0][2], 0);
	cameraPos += direction * distance * right;
}

void moveCameraUp(int direction, float distance) {
	vec4 up(rotationMat[1][0], rotationMat[1][1], rotationMat[1][2], 0);
	cameraPos += direction * distance * up;
}

void moveCameraForward(int direction, float distance) {
	vec4 forward(rotationMat[2][0], rotationMat[2][1], rotationMat[2][2], 0);
	cameraPos += direction * distance * forward;
}

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
