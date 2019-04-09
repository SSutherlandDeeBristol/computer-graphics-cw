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

#define FULLSCREEN_MODE false
#define SCREEN_WIDTH 600
#define SCREEN_HEIGHT 600
#define CLIP_OFFSET 50
#define NEAR_CLIP_THRESHOLD 2
#define FAR_CLIP_THRESHOLD 10
#define TRIANGULATE true
#define FAN false
#define FILL true

SDL_Event event;
vector<Triangle> triangles;

const float focalLength = SCREEN_HEIGHT;
const vec4 defaultCameraPos(0.0, 0.0, -3.0, 1.0);
const vec4 defaultLightPos(0, -0.6, -0.2, 1);
const vec3 defaultReflectance(1.2, 1.2, 1.2);

const vec3 lightPower = 10.0f * vec3(1, 1, 1);
const vec3 indirectLightPowerPerArea = 0.5f * vec3(1, 1, 1);

vec4 lightPos(0, -0.6, -0.2, 1);
vec4 cameraPos(0, 0, -3.001, 1.0);

mat3 rotationMat; mat4 transformMat; mat4 projectionMat;
float yaw = 0; float pitch = 0;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

enum Axis {X = 0, Y = 1, Z = 2, W = 3};

struct Pixel {
  int x;
  int y;
  float zinv;
  vec4 pos3d;
  vec3 colour3d;
  vec4 normal3d;
  vec3 reflectance3d;
};

struct Vertex {
  vec4 position;
};


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);

bool IsWithinScreenBounds(Pixel p);
void Interpolate(float a, float b, vector<float>& result);
void InterpolatePixels(Pixel a, Pixel b, vector<Pixel>& result);

void ClipHomogenousVerticesToPlane(vector<Vertex>& clipBuffer, Axis axis, float maxVal, bool pos);
void ClipTriangle(Triangle triangle, vector<Vertex>& inVertices);
bool IsInsidePos(Vertex vertex, Axis axis, float maxVal);
bool IsInsideNeg(Vertex vertex, Axis axis, float maxVal);
void CalculateIntersection(Vertex a, Vertex b, Vertex& c, Axis axis, float maxVal);
void HomogenousFlatten(Vertex homogenousVertex, Vertex& flatVertex);
void HomogenousDivide(Vertex homogenousVertex, Vertex& projectedVertex);
void ProjectToHomogenous(Vertex vertex, mat4 proj, Vertex& projectedVertex);
void TriangulateVerticesFan(vector<Vertex> vertices, vector<Triangle>& result, vec3 colour);
float TriangulateVertices(vector<Vertex> vertices, vector<Triangle>& result, vec3 colour);
float TriangleCost(Triangle triangle);

void TransformationMatrix(vec4 camPos, mat3 rot, mat4 &T);
void UpdateRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void UpdateProjectionMatrix(mat4 &mat);

void VertexShader(const Vertex& v, Pixel& p, vec3 colour, vec4 normal, vec3 reflectance);
void PixelShader(screen *screen, const Pixel& p);
void DrawClipOffset(screen* screen, bool fillOutline);
void DrawDepthBuffer(screen* screen);
void DrawLineSDL(screen* surface, Pixel a, Pixel b, vec3 colour);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec3 colour, vec4 normal, vec3 reflectance, bool fill);
void DrawTriangle(screen* screen, Triangle triangle, vec3 reflectance, bool fill);

void MoveCameraRight(int direction, float distance);
void MoveCameraUp(int direction, float distance);
void MoveCameraForward(int direction, float distance);
void LookAt(vec3 toPos);
void ResetCameraPosition();
/* FUNCTIONS END                                                               */
/* ----------------------------------------------------------------------------*/

int main(int argc, char* argv[]) {

  screen *depthScreen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
  screen *mainScreen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);

  // LoadTestTriangleZ(triangles);
  LoadTestModel(triangles);
  // LoadTestTriangle(triangles);

  lightPos = defaultLightPos;

  while (Update()) {
    Draw(mainScreen);
    DrawClipOffset(mainScreen, false);
    SDL_Renderframe(mainScreen);

    DrawDepthBuffer(depthScreen);
    SDL_Renderframe(depthScreen);
  }

  SDL_SaveImage(mainScreen, "screenshot.bmp");

  KillSDL(mainScreen);
  KillSDL(depthScreen);
  return 0;
}

void VertexShader(const Vertex& v, Pixel& p, vec3 colour, vec4 normal, vec3 reflectance) {
  p.x = focalLength * (v.position.x / v.position.z) + SCREEN_WIDTH / 2;
  p.y = focalLength * (v.position.y / v.position.z) + SCREEN_HEIGHT / 2;
  p.zinv = 1 / v.position.z;
  p.pos3d = v.position;
  p.normal3d = normal;
  p.colour3d = colour;
  p.reflectance3d = reflectance;
}

void PixelShader(screen *screen, const Pixel& p) {
  if (p.zinv > depthBuffer[p.y][p.x]) {
    depthBuffer[p.y][p.x] = p.zinv;

    /* Transform light into camera space */
    vec4 transformedLight = transformMat * lightPos;

    /* Calculate the amount of reflected light at this pixel */
    float radius = distance(transformedLight, p.pos3d);
    float area = 4 * M_PI * pow(radius, 2);
    vec4 rHat = normalize(transformedLight - p.pos3d);
    vec4 nHat = normalize(p.normal3d);

    vec3 directLightPower = lightPower * max(dot(rHat, nHat), 0.0f);
    vec3 directLightPowerPerArea = directLightPower / area;
    vec3 reflectedLight = p.reflectance3d * (directLightPowerPerArea + indirectLightPowerPerArea);

    /* Draw pixel to screen (if possible) */
    if (IsWithinScreenBounds(p)) PutPixelSDL(screen, p.x, p.y, reflectedLight * p.colour3d);
  }
}

void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  /* Reset Z-buffer values */
  for (int y = 0; y < SCREEN_HEIGHT; y++)
    for (int x = 0; x < SCREEN_WIDTH; x++)
      depthBuffer[y][x] = 0;

  UpdateRotationMatrix(pitch, yaw, 0, rotationMat);
  TransformationMatrix(cameraPos, rotationMat, transformMat);
  UpdateProjectionMatrix(projectionMat);

  for (size_t i = 0; i < triangles.size(); i++) {
    /* Transform the triangle to camera space */
    vec4 tv0 = transformMat * triangles[i].v0;
    vec4 tv1 = transformMat * triangles[i].v1;
    vec4 tv2 = transformMat * triangles[i].v2;
    Triangle transformedTriangle = Triangle(tv0, tv1, tv2, triangles[i].color);

    vector<Vertex> clippedVertices;

    ClipTriangle(transformedTriangle, clippedVertices);

    if (TRIANGULATE) {
      vector<Triangle> triangulatedVertices;
      if (FAN) TriangulateVerticesFan(clippedVertices, triangulatedVertices, triangles[i].color);
      else TriangulateVertices(clippedVertices, triangulatedVertices, triangles[i].color);

      for (size_t j = 0; j < triangulatedVertices.size(); j++)
        DrawTriangle(screen, triangulatedVertices[j], defaultReflectance, FILL);
    } else {
      DrawPolygon(screen, clippedVertices, triangles[i].color, triangles[i].normal, defaultReflectance, FILL);
    }
  }
}

void ClipTriangle(Triangle triangle, vector<Vertex>& inVertices) {
  float maxX = SCREEN_WIDTH/2 - CLIP_OFFSET, maxY = SCREEN_HEIGHT/2 - CLIP_OFFSET;
  vector<Vertex> clipBuffer;
  inVertices.clear();

  vector<Vertex> vertices(3);
  vertices[0].position = triangle.v0; vertices[1].position = triangle.v1; vertices[2].position = triangle.v2;

  /* Project (already transformed) vertices into Homogenous space */
  for (size_t i = 0; i < vertices.size(); i++) {
    Vertex homogenousCoord = vertices[i];
    ProjectToHomogenous(vertices[i], projectionMat, homogenousCoord);
    clipBuffer.push_back(homogenousCoord);
  }

  /* Clip to each plane in turn, maintaining a buffer of clipped vertices */
  ClipHomogenousVerticesToPlane(clipBuffer, X, maxX, true);                 /* Right  */
  ClipHomogenousVerticesToPlane(clipBuffer, X, maxX, false);                /* Left   */
  ClipHomogenousVerticesToPlane(clipBuffer, Y, maxY, true);                 /* Top    */
  ClipHomogenousVerticesToPlane(clipBuffer, Y, maxY, false);                /* Bottom */
  ClipHomogenousVerticesToPlane(clipBuffer, Z, NEAR_CLIP_THRESHOLD, false); /* Near   */
  ClipHomogenousVerticesToPlane(clipBuffer, Z, FAR_CLIP_THRESHOLD, true);   /* Far    */

  /* Flatten all vertices in the buffer (set their w component to 1) */
  for (size_t i = 0; i < clipBuffer.size(); i++) {
    Vertex flattenedCoord = clipBuffer[i];
    HomogenousFlatten(clipBuffer[i], flattenedCoord);
    inVertices.push_back(flattenedCoord);
  }
}

void ClipHomogenousVerticesToPlane(vector<Vertex>& clipBuffer, Axis axis, float maxVal, bool pos) {
  vector<Vertex> homogenousVertices = clipBuffer;
  size_t vs = homogenousVertices.size();
  clipBuffer.clear();

  for (size_t i = 0; i < vs; i++) {
    Vertex curVertex = homogenousVertices[i];
    Vertex nxtVertex = homogenousVertices[(i + 1) % vs];

    /* Check whether current and next vertices require clipping */
    bool isCurIn = pos ? IsInsidePos(curVertex, axis, maxVal) : IsInsideNeg(curVertex, axis, maxVal);
    bool isNxtIn = pos ? IsInsidePos(nxtVertex, axis, maxVal) : IsInsideNeg(nxtVertex, axis, maxVal);

    /* If the current vertex does not require clipping, push it unchanged */
    if (isCurIn) clipBuffer.push_back(curVertex);

    /* If either the current or the next vertex require clipping, push a new vertex on the boundary intersection */
    if ((isCurIn && !isNxtIn) || (!isCurIn && isNxtIn)) {
      Vertex newVertex;
      float val = (axis == Z) ? maxVal : (pos ? maxVal : -maxVal);
      CalculateIntersection(curVertex, nxtVertex, newVertex, axis, val);
      clipBuffer.push_back(newVertex);
    }
  }
}

/* Flattens a homogenous vertex (sets w to 1) */
void HomogenousFlatten(Vertex homogenousVertex, Vertex& flatVertex) {
  flatVertex = homogenousVertex;
  flatVertex.position.w = 1;
}

/* Performs a homogenous divide on a vertex, projecting to the 3D plane w = 1 */
void HomogenousDivide(Vertex homogenousVertex, Vertex& projectedVertex) {
  projectedVertex = homogenousVertex;
  projectedVertex.position = (1 / homogenousVertex.position.w) * homogenousVertex.position;
}

/* Projects a 3D point into homogenous space */
void ProjectToHomogenous(Vertex vertex, mat4 proj, Vertex& projectedVertex) {
  projectedVertex = vertex;
  projectedVertex.position = proj * vertex.position;
}

void CalculateIntersection(Vertex a, Vertex b, Vertex& c, Axis axis, float maxVal) {
  float t;

  if (axis == Z) {
    t = (maxVal - a.position.z) / (b.position.z - a.position.z);
  } else {
    t = (a.position[axis] - maxVal * a.position.w) /
        ((maxVal * (b.position.w - a.position.w)) - (b.position[axis] - a.position[axis]));
  }

  c.position = a.position + t * (b.position - a.position);
}

bool IsInsidePos(Vertex vertex, Axis axis, float maxVal) {
  return (axis == Z) ? vertex.position.z <= maxVal : (vertex.position[axis] <= (vertex.position.w * maxVal));
}

bool IsInsideNeg(Vertex vertex, Axis axis, float maxVal) {
  return (axis == Z) ? vertex.position.z >= maxVal : (vertex.position[axis] >= (vertex.position.w * -maxVal));
}

bool IsWithinScreenBounds(Pixel p) {
  return p.x > 0 && p.x < SCREEN_WIDTH && p.y > 0 && p.y < SCREEN_HEIGHT;
}

/* Triangulates using the fan approach */
void TriangulateVerticesFan(vector<Vertex> vertices, vector<Triangle>& result, vec3 colour) {
  for (size_t v = 1; v < vertices.size() - 1; v++) {
    vec4 tv0 = vertices[0].position;
    vec4 tv1 = vertices[v].position;
    vec4 tv2 = vertices[v+1].position;
    result.push_back(Triangle(tv0, tv1, tv2, colour));
  }
}

/* Triangulates using the optimum approach */
float TriangulateVertices(vector<Vertex> vertices, vector<Triangle>& result, vec3 colour) {
  int vertexCount = vertices.size();
  float cost = numeric_limits<float>::max();

  if (vertexCount < 3) return 0;

  for (int v = 1; v < vertexCount - 1; v++) {
    vec4 tv0 = vertices[0].position;
    vec4 tv1 = vertices[v].position;
    vec4 tv2 = vertices[vertexCount - 1].position;

    Triangle currentTriangle = Triangle(tv0, tv1, tv2, colour);

    vector<Vertex> previousVertices, nextVertices;
    vector<Triangle> previousResults, nextResults;

    for (int i = 0; i <= v; i++) previousVertices.push_back(vertices[i]);
    for (int i = v; i < vertexCount; i++) nextVertices.push_back(vertices[i]);

    float costPrevious = TriangulateVertices(previousVertices, previousResults, colour);
    float costNext = TriangulateVertices(nextVertices, nextResults, colour);
    float costCur = TriangleCost(currentTriangle);
    float netCost = costPrevious + costNext + costCur;

    if (netCost < cost) {
      cost = netCost;
      result.clear();
      result.insert(result.end(), previousResults.begin(), previousResults.end());
      result.push_back(currentTriangle);
      result.insert(result.end(), nextResults.begin(), nextResults.end());
    }
  }
  return cost;
}

float TriangleCost(Triangle triangle) {
  return distance(triangle.v0, triangle.v1) + distance(triangle.v1, triangle.v2) + distance(triangle.v2, triangle.v0);
}

void InterpolatePixels(Pixel a, Pixel b, vector<Pixel>& result) {
  int size = result.size();
  vector<float> x(size), y(size), zinv(size);
  vector<float> x3d(size), y3d(size), z3d(size);

  /* Interpolate 2D values */
  Interpolate(a.x, b.x, x); Interpolate(a.y, b.y, y); Interpolate(a.zinv, b.zinv, zinv);

  /* Interpolate 3D values with perspective correction */
  Interpolate(a.pos3d.x / a.pos3d.z, b.pos3d.x / b.pos3d.z, x3d);
  Interpolate(a.pos3d.y / a.pos3d.z, b.pos3d.y / b.pos3d.z, y3d);
  Interpolate(a.pos3d.z, b.pos3d.z, z3d);

  /* Populate results vector */
  for (int i = 0; i < size; i++) {
    result[i].x = round(x[i]); result[i].y = round(y[i]); result[i].zinv = zinv[i];

    result[i].pos3d = vec4(x3d[i] * z3d[i], y3d[i] * z3d[i], z3d[i], 1);
    result[i].colour3d = a.colour3d;
    result[i].normal3d = a.normal3d;
    result[i].reflectance3d = a.reflectance3d;
  }
}

void Interpolate(float a, float b, vector<float>& result) {
  int size = result.size();
  float dist = b - a, incr = dist / float(max(1, (size - 1)));

  for (int i = 0; i < size; i++) result[i] = a + (incr * i);
}

void DrawTriangle(screen* screen, Triangle triangle, vec3 reflectance, bool fill) {
  /* Convert triangle to list of vertices, then draw as usual */
  vector<Vertex> vertices;
  Vertex v0, v1, v2;
  v0.position = triangle.v0; v1.position = triangle.v1; v2.position = triangle.v2;
  vertices.push_back(v0); vertices.push_back(v1); vertices.push_back(v2);

  DrawPolygon(screen, vertices, triangle.color, triangle.normal, reflectance, fill);
}

void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec3 colour, vec4 normal, vec3 reflectance, bool fill) {
  int V = vertices.size();
  vector<Pixel> projectedVertices(V);

  for (int i = 0; i < V; i++) VertexShader(vertices[i], projectedVertices[i], colour, normal, reflectance);

  if (fill) {
    vector<Pixel> leftPixels; vector<Pixel> rightPixels;
    ComputePolygonRows(projectedVertices, leftPixels, rightPixels);
    DrawPolygonRows(screen, leftPixels, rightPixels);
  } else {
    /* Loop over all vertices and draw the edge from it to the next vertex */
    for (int i = 0; i < V; i++) {
      int j = (i + 1) % V; // The next vertex
      DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], colour);
    }
  }
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

    int deltaX = abs(nxtVertex.x - curVertex.x);
    int deltaY = abs(nxtVertex.y - curVertex.y);

    int lineLength = (deltaX > deltaY) ? deltaX : deltaY;

    vector<Pixel> pixels(lineLength + 1);

    InterpolatePixels(curVertex, nxtVertex, pixels);

    for (size_t j = 0; j < pixels.size(); j++) {
      Pixel curPixel = pixels[j];
      int yindex = curPixel.y - minY;

      if (curPixel.x <= leftPixels[yindex].x) leftPixels[yindex] = curPixel;
      if (curPixel.x >= rightPixels[yindex].x) rightPixels[yindex] = curPixel;
    }
  }
}

void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {
  if (leftPixels.size() != rightPixels.size()) {
    std::cout << "Error in DrawPolygonRows: Left =/= Right" << std::endl;
    return;
  }

  for (size_t i = 0; i < leftPixels.size(); i++) {
    vector<Pixel> pixels(rightPixels[i].x - leftPixels[i].x + 1);
    InterpolatePixels(leftPixels[i], rightPixels[i], pixels);

    for (size_t j = 0; j < pixels.size(); j++) {
      PixelShader(screen, pixels[j]);
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

void DrawDepthBuffer(screen* screen) {
  for (int y = 0; y < SCREEN_HEIGHT; y++) {
    for (int x = 0; x < SCREEN_WIDTH; x++) {
      vec3 colour = vec3(depthBuffer[y][x], depthBuffer[y][x], depthBuffer[y][x]);
      PutPixelSDL(screen, x, y, colour);
    }
  }
}

void DrawLineSDL(screen* screen, Pixel a, Pixel b, vec3 colour) {
  ivec2 delta = ivec2(glm::abs(a.x - b.x), glm::abs(a.y - b.y));
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<Pixel> line(pixels);
  InterpolatePixels(a, b, line);

  for (int px = 0; px < pixels; px++) {
    Pixel pixel = line[px];
    if (IsWithinScreenBounds(pixel)) PutPixelSDL(screen, pixel.x, pixel.y, colour);
  }
}

void UpdateProjectionMatrix(mat4 &mat) {
  mat = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 1/focalLength), vec4(0, 0, 0, 0));
}

void UpdateRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R) {
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

void MoveCameraRight(int direction, float distance) {
	vec4 right(rotationMat[0][0], rotationMat[1][0], rotationMat[2][0], 0);
	cameraPos += direction * distance * right;
}

void MoveCameraUp(int direction, float distance) {
	vec4 up(rotationMat[0][1], rotationMat[1][1], rotationMat[2][1], 0);
	cameraPos += direction * distance * up;
}

void MoveCameraForward(int direction, float distance) {
	vec4 forward(rotationMat[0][2], rotationMat[1][2], rotationMat[2][2], 0);
	cameraPos += direction * distance * forward;
}

void LookAt(vec3 toPos) {
	vec3 fromPos = vec3(cameraPos);

	vec3 forward = normalize(toPos - fromPos);

  pitch = asin(-forward.y);
  yaw = atan2(forward.x, forward.z);

  UpdateRotationMatrix(pitch, yaw, 0, rotationMat);
}

void ResetCameraPosition() {
  cameraPos = defaultCameraPos;
  pitch = 0;
  yaw = 0;
}

bool Update() {
  static int t = SDL_GetTicks();

  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2 - t);
  t = t2;

  std::cout << "Render time: " << dt << " ms." << std::endl;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
  if (e.type == SDL_QUIT) {
    return false;
  } else if (e.type == SDL_KEYDOWN) {
    int key_code = e.key.keysym.sym;
    switch(key_code) {
      case SDLK_UP:     pitch += M_PI / 18; break;
      case SDLK_DOWN:   pitch -= M_PI / 18; break;
      case SDLK_LEFT:   yaw -= M_PI / 18; break;
      case SDLK_RIGHT:  yaw += M_PI / 18; break;
      case SDLK_r:      LookAt(vec3(0, 0, 0)); break;
      case SDLK_t:      ResetCameraPosition(); break;
      case SDLK_w:      MoveCameraUp(-1, 0.25); break;
      case SDLK_s:      MoveCameraUp(1, 0.25); break;
      case SDLK_a:      MoveCameraRight(-1, 0.25); break;
      case SDLK_d:      MoveCameraRight(1, 0.25); break;
      case SDLK_EQUALS: MoveCameraForward(1, 0.25); break;
      case SDLK_MINUS:  MoveCameraForward(-1, 0.25); break;
      case SDLK_i:      lightPos.z += 0.2; break;
      case SDLK_k:      lightPos.z -= 0.2; break;
      case SDLK_j:      lightPos.x -= 0.2; break;
      case SDLK_l:      lightPos.x += 0.2; break;
      case SDLK_o:      lightPos.y -= 0.2; break;
      case SDLK_p:      lightPos.y += 0.2; break;
      case SDLK_ESCAPE: return false;
      }
    }
  }

  return true;
}
