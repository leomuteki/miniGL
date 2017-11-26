/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

// Global structs
struct Vertex;
struct Triangle;
struct Pixel;

// Global variables
vector<Vertex> listOfVertices;
vector<Triangle> listOfTriangles;
vector<mat4> modelViewStack;
vector<mat4> projectionStack;
MGLpoly_mode drawMode;
MGLmatrix_mode matMode;
vec3 MGLcolor = { 255.0f, 255.0f, 255.0f };

float areaRatio(int ax, int ay, int bx, int by, int cx, int cy) {
  return ax*(by-cy) + ay*(cx-bx) + (bx*cy-by*cx);
}

/*
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
  printf("%s\n", description);
  exit(1);
}

vector<mat4> currentStack() {
  if (matMode == MGL_PROJECTION) {
    return projectionStack;
  }
  else if (matMode == MGL_MODELVIEW) {
    return modelViewStack;
  }
  else {
    MGL_ERROR("Matrix mode was not set.");
    return NULL;
  }
}

// Global data structures
struct Vertex {
  vec3 color;
  vec4 position;
  Vertex() {
    color = { 0, 0, 0 };
    position = { 0, 0, 0, 0 };
  }
  Vertex(float x, float y, float z) {
    color = vec3(1.0f, 1.0f, 1.0f);
    position = { x, y, z, 1 };
  }
  Vertex(float r, float g, float b, float x, float y) {
    color = { r, g, b };
    position = { x, y, 0, 1 };
  }
  Vertex(float r, float g, float b, float x, float y, float z, float w) {
    color = { r, g, b };
    position = { x, y, z, w };
  }
};

struct Triangle {
  Vertex a, b, c;
  Triangle(Vertex aIn, Vertex bIn, Vertex cIn) {
    a = aIn;
    b = bIn;
    c = cIn;
  }
};

struct Pixel {
  vec3 color;
  vec3 position;
  Pixel() {
    color = { 0, 0, 0 };
    position = { 0, 0, 0 };
  }
  Pixel(float x, float y, float z) {
    color = vec3(1.0f, 1.0f, 1.0f);
    position = { x, y, z };
  }
  Pixel(float r, float g, float b, float x, float y) {
    color = { r, g, b };
    position = { x, y, 0};
  }
  Pixel(float r, float g, float b, float x, float y, float z) {
    color = { r, g, b };
    position = { x, y, z };
  }
};

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
  // Get every triangle from the global list
  unsigned listOfTrianglesSize = listOfTriangles.size();
  // Fill zBuffer with null Pixels for zValue comparison
  Pixel zBuffer[width][height];
  for (unsigned i = 0; i < width; ++i) {
    for (unsigned j = 0; j < height; ++j) {
      zBuffer[i][j] = Pixel();
    }
  }
  for (unsigned i = 0; i < listOfTrianglesSize; ++i) {
    Triangle curTri = listOfTriangles.at(i);
    // Transform to pixel coordinates
    int ax = (curTri.a.position[0] + 1) * width / 2;
    int ay = (curTri.a.position[1] + 1) * height / 2;
    int bx = (curTri.b.position[0] + 1) * width / 2;
    int by = (curTri.b.position[1] + 1) * height / 2;
    int cx = (curTri.c.position[0] + 1) * width / 2;
    int cy = (curTri.c.position[1] + 1) * height / 2;
    // Find rectangle tightly fit around triangle
    int minX = min(min(ax, bx), cx);
    int maxX = max(max(ax, bx), cx);
    int minY = min(min(ay, by), cy);
    int maxY = max(max(ay, by), cy);
    for (int i = minX; i < maxX; ++i) {
      for (int j = minY; j < maxY; ++j) {
        // Check if image pixel in within triangle
        float areaABC = areaRatio(ax, ay, bx, by, cx, cy);
        float alpha = areaRatio(i, j, bx, by, cx, cy) / areaABC;
        float beta = areaRatio(ax, ay, i, j, cx, cy) / areaABC;
        float gamma = 1 - alpha - beta;
        if (alpha >= 0 && beta >= 0 && gamma >= 0) {
          // push MGLpixel onto zBuffer
          if (zBuffer[i][j].position[2] < (curTri.a.position[2]*alpha+curTri.b.position[2]*beta+curTri.c.position[2]*gamma)) {
            zBuffer[i][j] = Pixel(i, j,(curTri.a.position[2]*alpha+curTri.b.position[2]*beta+curTri.c.position[2]*gamma));
          }
        }
      }
    }
  }
  // Push highest in each zBuffer to *data
  for (unsigned i = 0; i < width; ++i) {
    for (unsigned j = 0; j < height; ++j){
      // Fill in passed-in array
      *(data + (unsigned int)(zBuffer[i][j].position[0] + zBuffer[i][j].position[1]) * width) =
        Make_Pixel(255,255,255);

    }
  }
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */  
void mglBegin(MGLpoly_mode mode)
{
  drawMode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */

void mglEnd()
{
  if (drawMode == MGL_TRIANGLES) {
    for (unsigned i = 0; i < listOfVertices.size(); ++i) {
      if (i + 2 < listOfVertices.size()) {
        Vertex a, b, c;
        a = listOfVertices.at(i);
        b = listOfVertices.at(i+1);
        c = listOfVertices.at(i+2);
        listOfTriangles.push_back(Triangle(a, b, c));
      }
    }
  }
  else if (drawMode == MGL_QUADS) {
    for (unsigned i = 0; i < listOfVertices.size(); ++i) {
      if (i + 3 < listOfVertices.size()) {
        Vertex a, b, c, d;
        a = listOfVertices.at(i);
        b = listOfVertices.at(i+1);
        c = listOfVertices.at(i+2);
        d = listOfVertices.at(i+3);
        listOfTriangles.push_back(Triangle(a, b, c));
        listOfTriangles.push_back(Triangle(b, c, d));
      }
    }
  }
}


/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
  listOfVertices.push_back(Vertex(x, y, 0.0f));
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
  listOfVertices.push_back(Vertex(x, y, z));
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
  matMode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
  
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
  
}
