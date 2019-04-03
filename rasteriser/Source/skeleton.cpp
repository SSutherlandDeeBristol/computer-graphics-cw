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
mat3 rotation;
mat4 transformMat;
mat4 projection;
float yaw = 0;
float pitch = 0;

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
void VertexShader( const Vertex& v, Pixel& p );
bool isWithinBounds(ivec2 v);
bool isWithinBounds(Pixel p);
void TransformationMatrix(vec4 camPos, mat3 rot, mat4 &T);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result );
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void DrawLineSDL(screen* surface, Pixel a, Pixel b, vec3 colour);
void DrawPolygonEdges(screen* screen, const vector<Vertex>& vertices, vec3 colour);
void getProjectionMatrix(mat4 &mat);
void clipTriangle(vector<vec4> vertices, vector<vec4>& inVertices);
bool isInsidePos(vec4 vertex, int axis, float maxVal);
bool isInsideNeg(vec4 vertex, int axis, float maxVal);
void DrawClipOffset(screen* screen);
void calcIntersection(vec4 a, vec4 b, vec4& c, int axis, float maxVal);
void homogenousFlatten(vec4 homogenousVertex, vec4& flatVertex);
void homogenousDivide(vec4 homogenousVertex, vec4& projectedVertex);
void projectToHomogenous(vec4 vertex, mat4 proj, vec4& projectedVertex);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(screen* screen, const vector<Vertex>& vertices);

int main( int argc, char* argv[] ) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  LoadTestModel(triangles);
  // LoadTestTriangle(triangles);

  while ( Update() ) {
    Draw(screen);
    DrawClipOffset(screen);
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

void VertexShader(const Vertex& v, Pixel& p) {
  p.x = focalLength * (v.position.x / v.position.z) + SCREEN_WIDTH / 2;
  p.y = focalLength * (v.position.y / v.position.z) + SCREEN_HEIGHT / 2;
}

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  getRotationMatrix(pitch, yaw, 0, rotation);
  TransformationMatrix(cameraPos, rotation, transformMat);
  getProjectionMatrix(projection);

  for( uint32_t i=0; i<triangles.size(); ++i ) {
    vector<vec4> transformedVertices(3);
    transformedVertices[0] = transformMat * triangles[i].v0;
    transformedVertices[1] = transformMat * triangles[i].v1;
    transformedVertices[2] = transformMat * triangles[i].v2;

    vector<vec4> clippedVertices;

    clipTriangle(transformedVertices, clippedVertices);

    // TODO: Convert polygon to triangles and draw

    vector<Vertex> cvs;
    for (int i = 0; i < clippedVertices.size(); i++) {
      Vertex v;
      v.position = clippedVertices[i];
      cvs.push_back(v);
    }

    vec3 colour = vec3(1, 1, 1);
    DrawPolygonEdges(screen, cvs, colour);
    //DrawPolygon(screen, cvs);
  }
}

void clipTriangle(vector<vec4> vertices, vector<vec4>& inVertices) {
  float maxX = SCREEN_WIDTH/2 - CLIP_OFFSET, maxY = SCREEN_HEIGHT/2 - CLIP_OFFSET;
  vector<vec4> homogenousVertices, clipBuffer;

  /* Project co-ordinates into Homogenous space */
  int vs = vertices.size();
  for (int i = 0; i < vs; i++) {
    vec4 homogenousCoord;
    projectToHomogenous(vertices[i], projection, homogenousCoord);
    homogenousVertices.push_back(homogenousCoord);
  }

  vs = homogenousVertices.size();

  /* Clip Right */
  for (int i = 0; i < vs; i++) {
    vec4 curCoord = homogenousVertices[i];
    vec4 nextCoord = homogenousVertices[(i == vs - 1) ? 0 : i + 1];

    bool isCurInX = isInsidePos(curCoord, 0, maxX);
    bool isPrevInX = isInsidePos(nextCoord, 0, maxX);

    if (isCurInX) {
      clipBuffer.push_back(curCoord);
    }

    if ((isCurInX && !isPrevInX) || (!isCurInX && isPrevInX)) {
      /* Calculate Intersection I */
      vec4 newCoord;
      calcIntersection(curCoord, nextCoord, newCoord, 0, maxX);
      clipBuffer.push_back(newCoord);
    }
  }

  homogenousVertices = clipBuffer;
  vs = homogenousVertices.size();
  clipBuffer.clear();

  /* Clip Left */
  for (int i = 0; i < vs; i++) {
    vec4 curCoord = homogenousVertices[i];
    vec4 nextCoord = homogenousVertices[(i == vs -1 ) ? 0 : i + 1];

    bool isCurInX = isInsideNeg(curCoord, 0, maxX);
    bool isPrevInX = isInsideNeg(nextCoord, 0, maxX);

    if (isCurInX) {
      clipBuffer.push_back(curCoord);
    }

    if ((isCurInX && !isPrevInX) || (!isCurInX && isPrevInX)) {
      /* Calculate Intersection I */
      vec4 newCoord;
      calcIntersection(curCoord, nextCoord, newCoord, 0, -maxX);
      clipBuffer.push_back(newCoord);
    }
  }

  homogenousVertices = clipBuffer;
  vs = homogenousVertices.size();
  clipBuffer.clear();

  /* Clip Top */
  for (int i = 0; i < vs; i++) {
    vec4 curCoord = homogenousVertices[i];
    vec4 nextCoord = homogenousVertices[(i == vs -1 ) ? 0 : i + 1];

    bool isCurIn = isInsidePos(curCoord, 1, maxY);
    bool isPrevIn = isInsidePos(nextCoord, 1, maxY);

    if (isCurIn) {
      clipBuffer.push_back(curCoord);
    }

    if ((isCurIn && !isPrevIn) || (!isCurIn && isPrevIn)) {
      /* Calculate Intersection I */
      vec4 newCoord;
      calcIntersection(curCoord, nextCoord, newCoord, 1, maxY);
      clipBuffer.push_back(newCoord);
    }
  }

  homogenousVertices = clipBuffer;
  vs = homogenousVertices.size();
  clipBuffer.clear();

  /* Clip Bot */
  for (int i = 0; i < vs; i++) {
    vec4 curCoord = homogenousVertices[i];
    vec4 nextCoord = homogenousVertices[(i == vs -1 ) ? 0 : i + 1];

    bool isCurIn = isInsideNeg(curCoord, 1, maxY);
    bool isPrevIn = isInsideNeg(nextCoord, 1, maxY);

    if (isCurIn) {
      clipBuffer.push_back(curCoord);
    }

    if ((isCurIn && !isPrevIn) || (!isCurIn && isPrevIn)) {
      /* Calculate Intersection I */
      vec4 newCoord;
      calcIntersection(curCoord, nextCoord, newCoord, 1, -maxY);
      clipBuffer.push_back(newCoord);
    }
  }

  homogenousVertices = clipBuffer;
  vs = homogenousVertices.size();
  clipBuffer.clear();
  inVertices.clear();

  /* Flatten!! */
  for (int i = 0; i < vs; i++) {
    vec4 flattenedCoord;
    homogenousFlatten(homogenousVertices[i], flattenedCoord);
    inVertices.push_back(flattenedCoord);
  }
}

void homogenousFlatten(vec4 homogenousVertex, vec4& flatVertex) {
  flatVertex = homogenousVertex;
  flatVertex.w = 1;
}

void homogenousDivide(vec4 homogenousVertex, vec4& projectedVertex) {
  projectedVertex = (1/homogenousVertex.w) * homogenousVertex;
}

void projectToHomogenous(vec4 vertex, mat4 proj, vec4& projectedVertex) {
  projectedVertex = proj * vertex;
}

/* axis 0 = x, 1 = y, 2 = z, 3 = w */
void calcIntersection(vec4 a, vec4 b, vec4& c, int axis, float maxVal) {
  float t = (a[axis] - maxVal * a.w) / ((maxVal * (b.w - a.w)) - (b[axis] - a[axis]));
  c = a + t * (b - a);
}

/* axis 0 = x, 1 = y, 2 = z, 3 = w */
bool isInsidePos(vec4 vertex, int axis, float maxVal) {
  return (vertex[axis] <= (vertex.w * maxVal));
}

/* axis 0 = x, 1 = y, 2 = z, 3 = w */
bool isInsideNeg(vec4 vertex, int axis, float maxVal) {
  return (vertex[axis] >= (vertex.w * -maxVal));
}

bool isWithinBounds(ivec2 v) {
  return v.x > 0 && v.x < SCREEN_WIDTH && v.y > 0 && v.y < SCREEN_HEIGHT;
}

bool isWithinBounds(Pixel p) {
  return p.x > 0 && p.x < SCREEN_WIDTH && p.y > 0 && p.y < SCREEN_HEIGHT;
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

///////////////////////

void DrawPolygon(screen* screen, const vector<Vertex>& vertices) {
  int V = vertices.size();
  vector<Pixel> vertexPixels(V);

  for(int i=0; i<V; ++i) VertexShader(vertices[i], vertexPixels[i]);

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;

  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(screen, leftPixels, rightPixels);
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
  // 1. Find max and min y-value of the polygon
  //    and compute the number of rows it occupies.
  float maxY = -numeric_limits<float>::max();
  float minY = numeric_limits<float>::max();
  float numRows;

  if (vertexPixels.size() == 0) return;

  for (int i = 0; i < vertexPixels.size(); i++) {
    if (vertexPixels[i].y >= maxY) maxY = vertexPixels[i].y;
    if (vertexPixels[i].y <= minY) minY = vertexPixels[i].y;
  }

  numRows = maxY - minY + 1;

  // 2. Resize leftPixels and rightPixels
  //    so that they have an element for each row.
  leftPixels.resize(numRows);
  rightPixels.resize(numRows);

  // 3. Initialize the x-coordinates in leftPixels
  //    to some really large value and the x-coordinates
  //    in rightPixels to some really small value.
  for (int i = 0; i < leftPixels.size(); i++) {
    leftPixels[i].x = numeric_limits<int>::max();
    leftPixels[i].y = i + minY;
    rightPixels[i].x = -numeric_limits<int>::max();
    rightPixels[i].y = i + minY;
  }

  // 4. Loop through all edges of the polygon and use
  //    linear interpolation to find the x-coordinate for
  //    each row it occupies. Update the corresponding
  //    values in rightPixels and leftPixels.
  for (int i = 0; i < vertexPixels.size(); i++) {
    int nextVertex = (i == vertexPixels.size() - 1) ? 0 : i + 1;

    int deltaX = abs(vertexPixels[i].x - vertexPixels[nextVertex].x);
    int deltaY = abs(vertexPixels[i].y - vertexPixels[nextVertex].y);

    int lineLength = (deltaX > deltaY) ? deltaX : deltaY;

    vector<Pixel> pixels(lineLength+1);

    Interpolate(vertexPixels[i], vertexPixels[nextVertex], pixels);

    for (int j = 0; j < pixels.size(); j++) {
      int yindex = pixels[j].y - minY;

      if (pixels[j].x <= leftPixels[yindex].x) leftPixels[yindex] = pixels[j];
      if (pixels[j].x >= rightPixels[yindex].x) rightPixels[yindex] = pixels[j];
    }
  }
}

void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {
  for (int i = 0; i < leftPixels.size(); i++) {
    vector<Pixel> pixels(rightPixels[i].x - leftPixels[i].x + 1);
    Interpolate(leftPixels[i], rightPixels[i], pixels);

    for (int j = 0; j < pixels.size(); j++) {
      if (isWithinBounds(pixels[j])) PutPixelSDL(screen, pixels[j].x, pixels[j].y, vec3(1,1,1));
    }
  }
}
///////////////////////

void DrawClipOffset(screen* screen) {
  vec3 colour = vec3(1, 0, 0);

  Pixel TL, TR, BL, BR;

  TL.x = CLIP_OFFSET; TL.y = CLIP_OFFSET;
  TR.x = SCREEN_WIDTH - CLIP_OFFSET; TR.y = CLIP_OFFSET;
  BL.x = CLIP_OFFSET; BL.y = SCREEN_HEIGHT - CLIP_OFFSET;
  BR.x = SCREEN_WIDTH - CLIP_OFFSET; BR.y = SCREEN_HEIGHT - CLIP_OFFSET;

  DrawLineSDL(screen, TL, TR, colour);
  DrawLineSDL(screen, TR, BR, colour);
  DrawLineSDL(screen, BR, BL, colour);
  DrawLineSDL(screen, BL, TL, colour);
}

void DrawLineSDL(screen* screen, Pixel a, Pixel b, vec3 colour) {
  ivec2 delta = ivec2(glm::abs(a.x - b.x), glm::abs(a.y - b.y));
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<Pixel> line(pixels);
  Interpolate(a, b, line);

  for (int px = 0; px < pixels; px++) {
    Pixel pixel = line[px];
    if (isWithinBounds(pixel)) PutPixelSDL(screen, pixel.x, pixel.y, colour);
  }
}

void DrawPolygonEdges(screen* screen, const vector<Vertex>& vertices, vec3 colour) {
  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D image position:
  vector<Pixel> projectedVertices(V);
  for (int i = 0; i < V; ++i) {
    // /* Perform transform on co-ordinate */
    // vec4 transformedCoord = transformMat * vertices[i];
    // /* Project co-ordinate to Homogenous space */
    // vec4 homogenousCoord = projection * transformedCoord;
    // std::cout << glm::to_string(homogenousCoord) << std::endl;
    // /* Perform homogenous divide (projects to plane at w = 1) */
    // vec4 homogenousDivide = (1/homogenousCoord.w) * homogenousCoord;
    /* Get position on screen */
    // vec4 homogenousDivide = (1/vertices[i].w) * vertices[i]; // [u, v, f, 1]
    VertexShader(vertices[i], projectedVertices[i]);
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
