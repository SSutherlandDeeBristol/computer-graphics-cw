#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include "glm/gtx/string_cast.hpp" // std::cout<<glm::to_string(hello)<<std::endl;
#include <math.h>
#include "limits"

using namespace std;
using glm::ivec2;
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::mat3;
using glm::mat4;

struct Pixel {
  int x;
  int y;
  float zinv;
  vec4 pos3d;
};

struct Vertex {
  vec4 position;
};

SDL_Event event;

#define SCREEN_WIDTH 600
#define SCREEN_HEIGHT 600
#define FULLSCREEN_MODE false

#define CLIP_OFFSET 50

std::vector<Triangle> triangles;

const float focalLength = SCREEN_HEIGHT;
const vec4 defaultCameraPos(0.0, 0.0, -3.001, 1.0);

vec4 cameraPos(0, 0, -3.001, 1.0);
mat3 rotation;
mat4 transformMat;
float yaw = 0;
float pitch = 0;

vec3 currentColor;
vec4 currentNormal;
vec3 currentReflectance;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

vec4 defaultLightPos(0,-0.4,-0.3, 1);
vec4 lightPos(0,-0.3,-0.3, 1);
vec3 lightPower = 11.0f * vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader(const Vertex& v, Pixel& p);
void PixelShader(screen *screen, const Pixel& p);
bool isWithinBounds(ivec2 v);
bool isWithinBounds(Pixel v);
vector<Vertex> clipTriangleToScreen(vector<Vertex> vertices);
unsigned char ComputeOutcode(Pixel pixel);
void TransformationMatrix(vec4 camPos, mat3 rot, mat4 &T);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void Interpolate(Pixel a, Pixel b, vector<ivec2>& result);
void getRotationMatrix(float thetaX, float thetaY, float thetaZ, mat3 &R);
void DrawLineSDL(screen* surface, ivec2 a, ivec2 b, vec3 colour);
void DrawPolygonEdges(screen* screen, const vector<vec4>& vertices, vec3 colour);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(screen* screen, const vector<Vertex>& vertices);
void updateRotation();
void moveCameraForward(int direction, float distance);
void moveCameraUp(int direction, float distance);
void moveCameraRight(int direction, float distance);

unsigned char ComputeOutcode(Pixel p);

int main( int argc, char* argv[] ) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  /* Test outcodes */
  Pixel a, b, c, d, e, f, g, h, i;
  a.x = -1; a.y = -1;                              // 1001 = 9
  b.x = -1; b.y = SCREEN_HEIGHT / 2;               // 0001 = 1
  c.x = SCREEN_WIDTH / 2; c.y = -1;                // 1000 = 8
  d.x = SCREEN_WIDTH / 2; d.y = SCREEN_HEIGHT / 2; // 0000 = 0
  e.x = SCREEN_WIDTH + 1; e.y = -1;                // 1010 = 10
  f.x = SCREEN_WIDTH + 1; f.y = SCREEN_HEIGHT / 2; // 0010 = 2
  g.x = SCREEN_WIDTH + 1; g.y = SCREEN_HEIGHT + 1; // 0110 = 6
  h.x = SCREEN_WIDTH / 2; h.y = SCREEN_HEIGHT + 1; // 0100 = 4
  i.x = -1; i.y = SCREEN_HEIGHT + 1;               // 0101 = 5

  unsigned char acode = ComputeOutcode(a);
  unsigned char bcode = ComputeOutcode(b);
  unsigned char ccode = ComputeOutcode(c);
  unsigned char dcode = ComputeOutcode(d);
  unsigned char ecode = ComputeOutcode(e);
  unsigned char fcode = ComputeOutcode(f);
  unsigned char gcode = ComputeOutcode(g);
  unsigned char hcode = ComputeOutcode(h);
  unsigned char icode = ComputeOutcode(i);

  cout << "a:" << (int) acode << ", b:" << (int) bcode << ", c:" << (int) ccode << ", d:" << (int) dcode << "." << endl;
  cout << "e:" << (int) ecode << ", f:" << (int) fcode << ", g:" << (int) gcode << ", h:" << (int) hcode << ", i:" << (int) icode << "." << endl;
  /* End outcodes test */


  LoadTestModel(triangles);

  lightPos = defaultLightPos;

  while ( Update() ) {
    Draw(screen);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

void VertexShader(const Vertex& v, Pixel& p) {
  vec4 pos = transformMat*v.position;
  p.zinv = 1/pos.z;
  p.x = focalLength * (pos.x * p.zinv) + SCREEN_WIDTH / 2;
  p.y = focalLength * (pos.y * p.zinv) + SCREEN_HEIGHT / 2;
  p.pos3d = v.position;
}

void PixelShader(screen *screen, const Pixel& p) {
  int x = p.x;
  int y = p.y;

  if (p.zinv > depthBuffer[y][x]) {
    depthBuffer[y][x] = p.zinv;

    vec4 rVec = lightPos - p.pos3d;
    float r = distance(lightPos, p.pos3d);
    vec4 rHat = normalize(rVec);
    vec4 nHat = normalize(currentNormal);

    float A = float(4 * M_PI) * pow(r,2);

    vec3 B = (lightPower * max(dot(rHat, nHat), 0.0f));

    vec3 D = B / A;

    vec3 R = currentReflectance * (D + indirectLightPowerPerArea);

    PutPixelSDL(screen, x, y, R * currentColor);
  }
}

/* Place your drawing here */
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  getRotationMatrix(pitch, yaw, 0, rotation);
  TransformationMatrix(cameraPos, rotation, transformMat);

  for(int y=0; y<SCREEN_HEIGHT; ++y)
    for(int x=0; x<SCREEN_WIDTH; ++x)
      depthBuffer[y][x] = 0;

  for( uint32_t i=0; i<triangles.size(); ++i ) {
    if (dot(transformMat * cameraPos, triangles[i].normal) > 0) {
      vector<Vertex> vertices(3);
      vertices[0].position = triangles[i].v0;
      vertices[1].position = triangles[i].v1;
      vertices[2].position = triangles[i].v2;

      currentColor = triangles[i].color;
      currentNormal = triangles[i].normal;
      currentReflectance = vec3(1.2,1.2,1.2);

      DrawPolygon(screen, vertices);
    }
  }
}

bool isWithinBounds(ivec2 v) {
  return v.x > 0 && v.x < SCREEN_WIDTH && v.y > 0 && v.y < SCREEN_HEIGHT;
}

bool isWithinBounds(Pixel v) {
  return v.x > 0 && v.x < SCREEN_WIDTH && v.y > 0 && v.y < SCREEN_HEIGHT;
}

vector<Pixel> clipTriangleToScreen(vector<Pixel> vertexPixels) {
  int xmin = CLIP_OFFSET, xmax = SCREEN_WIDTH - CLIP_OFFSET,
      ymin = CLIP_OFFSET, ymax = SCREEN_HEIGHT - CLIP_OFFSET;

  int V = vertexPixels.size();
  vector<Pixel> invertexPixels;

  for (int i = 0; i < V; i++) {
    Pixel vertex = vertexPixels[i];
    Pixel nextVertex = vertexPixels[(i == V - 1) ? 0 : i + 1];

    // calculate outcodes for each projected vertex
    unsigned char p1code = ComputeOutcode(vertex);
    unsigned char p2code = ComputeOutcode(nextVertex);

    unsigned char outcode = p1code | p2code;

    printf("p1code: %d\n", p1code);
    printf("p2code: %d\n", p2code);

    // check if on-screen 0000 V 0000 = 0000
    if ((outcode) == 0) {
      printf("appending vertex\n");
      // add the first vertex to invertexPixels
      invertexPixels.push_back(vertex);

    // check if the whole line is off the screen and only continue if not
    } else if ((p1code & p2code) == 0) {
      // check top plane
      if ((outcode | 1<<3) == 1<<3) {
        // need to clip to top

        // calculate x value
        int x = vertex.x + (nextVertex.x - vertex.x) * ((ymin - vertex.y) / (nextVertex.y - vertex.y));

        // BROKEN: calculate pos3d
        vec4 pos = vertex.pos3d + (nextVertex.pos3d - vertex.pos3d) * ((ymin - vertex.pos3d.y) / (nextVertex.pos3d.y - vertex.pos3d.y));

        Pixel p;
        p.x = x;
        p.y = ymin;
        p.pos3d = pos;
        p.zinv = 1/pos.z;

        // add pixel to invertexPixels
        invertexPixels.push_back(p);

        // if the first vertex is in then add it to invertexPixels
        if (p1code == 0) invertexPixels.push_back(vertex);
      }
    }
  }

  return invertexPixels;
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result ) {
  int N = result.size();
  vec2 step = glm::vec2(b - a) / float(max(N - 1, 1));
  vec2 current(a);

  for (int i = 0; i < N; ++i) {
    result[i] = round(current);
    current += step;
  }
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result ) {
  int N = result.size();

  float stepX = (b.x - a.x) / float(max(N - 1, 1));
  float stepY = (b.y - a.y) / float(max(N - 1, 1));
  float stepZInv = (b.zinv - a.zinv) / float(max(N - 1, 1));
  vec4 stepPos = ((b.pos3d * b.zinv) - (a.pos3d * a.zinv)) / float(max(N - 1, 1));

  float currentX(a.x);
  float currentY(a.y);
  float currentZInv(a.zinv);
  vec4 currentPos(a.pos3d * a.zinv);

  for (int i = 0; i < N; ++i) {
    result[i].x = round(currentX);
    result[i].y = round(currentY);
    result[i].zinv = currentZInv;
    result[i].pos3d = vec4(currentPos.x / currentZInv, currentPos.y / currentZInv, currentPos.z / currentZInv, 1);
    currentX += stepX;
    currentY += stepY;
    currentZInv += stepZInv;
    currentPos += stepPos;
  }
}

void DrawLineSDL(screen* screen, Pixel a, Pixel b, vec3 colour) {
  ivec2 delta(abs(a.x - b.x), abs(a.y - b.y));
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<Pixel> line(pixels);
  Interpolate(a, b, line);

  for (int px = 0; px < pixels; px++) {
    Pixel pixel = line[px];
    if (isWithinBounds(pixel)) PixelShader(screen, pixel);
  }
}

void DrawPolygonEdges(screen* screen, const vector<Vertex>& vertices, vec3 colour) {
  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D image position:
  vector<Pixel> projectedVertices(V);
  for (int i = 0; i < V; ++i) {
    VertexShader(vertices[i], projectedVertices[i]);
  }
  // Loop over all vertices and draw the edge from it to the next vertex:
  for (int i = 0; i < V; ++i) {
    int j = (i + 1) % V; // The next vertex
    DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], colour);
  }
}

void DrawPolygon(screen* screen, const vector<Vertex>& vertices) {
  int V = vertices.size();
  vector<Pixel> vertexPixels(V);

  for(int i=0; i<V; ++i) {
    VertexShader(vertices[i], vertexPixels[i]);
    cout << "x: " << vertexPixels[i].x << ", y:" << vertexPixels[i].y << endl;
  }

  vector<Pixel> invertexPixels = clipTriangleToScreen(vertexPixels);

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;

  if (invertexPixels.size() > 0) {
    ComputePolygonRows(invertexPixels, leftPixels, rightPixels);

    //ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
    DrawPolygonRows(screen, leftPixels, rightPixels);
  }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
  // 1. Find max and min y-value of the polygon
  //    and compute the number of rows it occupies.
  float maxY = -numeric_limits<float>::max();
  float minY = numeric_limits<float>::max();
  float numRows;

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
      if (isWithinBounds(pixels[j])) PixelShader(screen, pixels[j]);
    }
  }
}

unsigned char ComputeOutcode(Pixel pixel) {
  unsigned char outcode = 0;
  int xmin = CLIP_OFFSET, xmax = SCREEN_WIDTH - CLIP_OFFSET,
      ymin = CLIP_OFFSET, ymax = SCREEN_HEIGHT - CLIP_OFFSET;

  if (pixel.x  < xmin) outcode = outcode | 1;
  if (pixel.x > xmax) outcode = outcode | 2;
  if (pixel.y < ymin) outcode = outcode | 8;
  if (pixel.y > ymax) outcode = outcode | 4;

  return outcode;
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

void TransformationMatrix(vec4 camPos, mat3 rot, mat4&T) {
  //mat4 m1 = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), camPos);
  //m1[3][3] = 1;

  mat4 m2 = mat4(rot);

  mat4 m3 = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), -camPos);
  m3[3][3] = 1;

  T = m2 * m3;
}

void updateRotation() {
  getRotationMatrix(pitch, yaw, 0, rotation);

  TransformationMatrix(cameraPos, rotation, transformMat);
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

  std::cout << "Render time: " << dt << " ms." << std::endl;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
      return false;
    } else if (e.type == SDL_KEYDOWN) {
      int key_code = e.key.keysym.sym;
      switch(key_code) {
        case SDLK_UP:
          pitch += M_PI / 18;
          updateRotation();
          break;
        case SDLK_DOWN:
          pitch -= M_PI / 18;
          updateRotation();
          break;
        case SDLK_LEFT:
          yaw -= M_PI / 18;
          updateRotation();
          break;
        case SDLK_RIGHT:
          yaw += M_PI / 18;
          updateRotation();
          break;
        case SDLK_r:
          /* Look-At function, points camera to 0,0,0 */
          // lookAt(ctw);
          break;
        case SDLK_t:
          // Reset camera position
          cameraPos = defaultCameraPos;
          lightPos = defaultLightPos;
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
        case SDLK_i:
          lightPos.z += 0.2;
          break;
        case SDLK_k:
          lightPos.z -= 0.2;
          break;
        case SDLK_j:
          lightPos.x -= 0.2;
          break;
        case SDLK_l:
          lightPos.x += 0.2;
          break;
        case SDLK_o:
          lightPos.y -= 0.2;
          break;
        case SDLK_p:
          lightPos.y += 0.2;
          break;
        case SDLK_ESCAPE:
          return false;
        }
      }
  }

  return true;
}
