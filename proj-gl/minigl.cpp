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
#include <math.h>
#include <vector>
#include <cstdio>
#include <stdlib.h>
#include <limits>

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
vec3 currentColor = { 255, 255, 255 };

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

vector<mat4>& currentStack() {
  if (matMode == MGL_PROJECTION) {
    return projectionStack;
  }
  else if (matMode == MGL_MODELVIEW) {
    return modelViewStack;
  }
  else {
    MGL_ERROR("Matrix mode was not set.");
    exit(1);
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
  MGLfloat z;
  Pixel() {
    color = { 0, 0, 0 };
    z = INFINITY;
  }
  Pixel(float r, float g, float b, float z) {
    color = { r, g, b };
    this->z = z;
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
        if (alpha >= -0 && beta >= -0 && gamma >= -0) {
          // push MGLpixel onto zBuffer
          MGLfloat curZValue = curTri.a.position[2]*alpha
            +curTri.b.position[2]*beta
            +curTri.c.position[2]*gamma;
          // Check if the current z value is smaller than the smallest
          // in the zbuffer, then update it if necessary
          if (zBuffer[i][j].z == INFINITY || curZValue < zBuffer[i][j].z) {
            zBuffer[i][j] =
              Pixel(currentColor[0], currentColor[1], currentColor[2], curZValue);
            zBuffer[i][j].color[0] = currentColor[0];
            zBuffer[i][j].color[1] = currentColor[1];
            zBuffer[i][j].color[2] = currentColor[2];
          }
        }
      }
    }
  }
  // Push highest in each zBuffer to *data
  for (MGLsize i = 0; i < width; ++i) {
    for (MGLsize j = 0; j < height; ++j){
      MGLfloat red = zBuffer[i][j].color[0];
      MGLfloat green = zBuffer[i][j].color[1];
      MGLfloat blue = zBuffer[i][j].color[2];
      *(data + i + j * width) = Make_Pixel(red, green, blue);
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
  listOfVertices.push_back(Vertex(x, y, INFINITY));
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
  currentStack().push_back(currentStack().back());
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
  currentStack().pop_back();
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
  // Load with identity matrix
  currentStack().push_back({ 1,0,0,0,
                             0,1,0,0,
                             0,0,1,0,
                             0,0,0,1 });
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
  mat4 newMat;
  for (unsigned int i = 0; i < 4; ++i) {
    for (unsigned int j = 0; j < 4; ++j) {
      newMat(i, j) = *(matrix + i + j*4);
    }
  }
  currentStack().push_back(newMat);
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
  mat4 newMat;
  for (unsigned int i = 0; i < 4; ++i) {
    for (unsigned int j = 0; j < 4; ++j) {
      newMat(i, j) = *(matrix + i + j*4);
    }
  }
  currentStack().back() = currentStack().back() * newMat;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
  mat4 translateMat = {1,0,0,0,
                    0,1,0,0,
                    0,0,1,0,
                    x,y,z,1};
  mglMultMatrix(translateMat.values);
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
  vec3 rotationAxis = vec3(x,y,z);
	rotationAxis = rotationAxis.normalized();
	MGLfloat cosTheta = cos(angle * (M_PI/180));
	MGLfloat sinTheta = sin(angle * (M_PI/180));
	MGLfloat xNorm= rotationAxis[0];
	MGLfloat yNorm = rotationAxis[1];
	MGLfloat zNorm = rotationAxis[2];
	mat4 rotationMat = {xNorm*xNorm*(1-cosTheta)+cosTheta, yNorm*xNorm*(1-cosTheta)+(zNorm*sinTheta), xNorm*zNorm*(1-cosTheta)-(yNorm*sinTheta), 0,
                 xNorm*yNorm*(1-cosTheta)-(zNorm*sinTheta), yNorm*yNorm*(1-cosTheta)+cosTheta, yNorm*zNorm*(1-cosTheta)+(xNorm*sinTheta), 0,
                 xNorm*zNorm*(1-cosTheta)+(yNorm*sinTheta), yNorm*zNorm*(1-cosTheta)-(xNorm*sinTheta), zNorm*zNorm*(1-cosTheta)+cosTheta, 0,
                 0,0,0,1};
	mglMultMatrix(rotationMat.values);
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
  mat4 scaleMat = {x,0,0,0,
                   0,y,0,0,
                   0,0,z,0,
                   0,0,0,1};
	mglMultMatrix(scaleMat.values);
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
  if (near == 0)
    {
      near = numeric_limits<float>::min(); 
    }
	mat4 frust = {(2*near)/(right-left), 0, 0, 0,
                0, (2 * near)/(top-bottom), 0, 0,
                (right+left)/(right-left), (top+bottom)/(top-bottom), -((far + near)/(far-near)), -1,
                0,0, -((2*far * near)/(far - near)), 0};
  //                 {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
	//topofactivestack() = topofactivestack() * frust;
	mglMultMatrix(frust.values);

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
  mat4 ortho = {2/(right-left), 0, 0, 0, 0, 2/(top-bottom), 0, 0, 0, 0, -2/(far-near), 0,-((right+left)/(right-left)),-((top+bottom)/(top-bottom)), -((far + near)/(far - near)),1 };
	//                 {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
	
	//cout << topofactivestack() << endl;
	mglMultMatrix(ortho.values);

}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
  currentColor[0] = red*255;
  currentColor[1] = green*255;
  currentColor[2] = blue*255;
}
