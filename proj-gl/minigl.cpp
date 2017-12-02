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
//#define INFINITY 0

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
  Vertex(MGLfloat x, MGLfloat y, MGLfloat z) {
    color = vec3(1.0f, 1.0f, 1.0f);
    position = { x, y, z, 1 };
  }
  Vertex(MGLfloat r, MGLfloat g, MGLfloat b, MGLfloat x, MGLfloat y) {
    color = { r, g, b };
    position = { x, y, 0, 1 };
  }
  Vertex(MGLfloat r, MGLfloat g, MGLfloat b, MGLfloat x, MGLfloat y, MGLfloat z) {
    color = { r, g, b };
    position = { x, y, z, 1 };
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
  Pixel(MGLfloat r, MGLfloat g, MGLfloat b, MGLfloat z) {
    color = { r, g, b };
    this->z = z;
  }
};

float areaRatio(int ax, int ay, int bx, int by, int cx, int cy) {
  return ax*(by-cy) + ay*(cx-bx) + (bx*cy-by*cx);
}

void transform(vec4 &vec) {
  mat4 mod=modelViewStack.back();
  mat4 proj=projectionStack.back();
  vec = proj * mod * vec;
}

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
  MGLsize listOfTrianglesSize = listOfTriangles.size();
  // Fill zBuffer with null Pixels for zValue comparison
  Pixel zBuffer[width][height];
  for (MGLsize i = 0; i < width; ++i) {
    for (MGLsize j = 0; j < height; ++j) {
      zBuffer[i][j] = Pixel();
    }
  }
  for (MGLsize i = 0; i < listOfTrianglesSize; ++i) {
    Triangle curTri = listOfTriangles.at(i);
    MGLfloat wa = curTri.a.position[3];
    curTri.a.position[0] /= wa;
    curTri.a.position[1] /= wa;
    curTri.a.position[2] /= wa;
    MGLfloat wb = curTri.b.position[3];
    curTri.b.position[0] /= wb;
    curTri.b.position[1] /= wb;
    curTri.b.position[2] /= wb;
    MGLfloat wc = curTri.c.position[3];
    curTri.c.position[0] /= wc;
    curTri.c.position[1] /= wc;
    curTri.c.position[2] /= wc;
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
    if (minX < 0) minX = 0;
    if (maxX >= (int)width) maxX = width-1;
    if (minY < 0) minY = 0;
    if (maxY >= (int)height) maxY = height-1;
    for (int i = minX; i < maxX; ++i) {
      for (int j = minY; j < maxY; ++j) {
        // Check if image pixel in within triangle
        MGLfloat areaABC = areaRatio(ax, ay, bx, by, cx, cy);
        MGLfloat alpha0 = areaRatio(i, j, bx, by, cx, cy) / areaABC;
        MGLfloat beta0 = areaRatio(ax, ay, i, j, cx, cy) / areaABC;
        MGLfloat gamma0 = 1 - alpha0 - beta0;
        // Apply perspective transformations on alpha,beta,gamma
        MGLfloat baryPerspectiveSum=((alpha0/wa)+(beta0/wb)+(gamma0/wc));
        MGLfloat alpha=(alpha0/wa)/baryPerspectiveSum;
        MGLfloat beta=(beta0/wb)/baryPerspectiveSum;
        MGLfloat gamma=1-alpha-beta;
        // Calculate interpolated z
        MGLfloat interpolatedZ=curTri.a.position[2]*alpha+curTri.b.position[2]*beta+curTri.c.position[2]*gamma;
        // Check alpha, beta, gamma to see whether to draw the pixel
        if (alpha >= 0 && beta >= 0 && gamma >= 0 && (interpolatedZ==1||interpolatedZ==-1||(interpolatedZ<1&&interpolatedZ>-1))) {
          // in the zbuffer, then update it if necessary
          if (zBuffer[i][j].z == INFINITY || interpolatedZ <= zBuffer[i][j].z) {
            // Calulate color at pixel position
            MGLfloat r = alpha*curTri.a.color[0]+beta*curTri.b.color[0]+gamma*curTri.c.color[0];
            MGLfloat g = alpha*curTri.a.color[1]+beta*curTri.b.color[1]+gamma*curTri.c.color[1];
            MGLfloat b = alpha*curTri.a.color[2]+beta*curTri.b.color[2]+gamma*curTri.c.color[2];

            zBuffer[i][j] =
              Pixel(r,g,b,interpolatedZ);
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
  for ( MGLsize i = 0; i <listOfVertices.size(); ++i) {
  }
  if (drawMode == MGL_TRIANGLES) {
    for (MGLsize i = 0; i < listOfVertices.size(); i+=3) {
        Vertex a, b, c;
        a = listOfVertices.at(i);
        b = listOfVertices.at(i+1);
        c = listOfVertices.at(i+2);
        listOfTriangles.push_back(Triangle(a, b, c));
    }
  }
  else if (drawMode == MGL_QUADS) {
    for (MGLsize i = 0; i < listOfVertices.size(); i+=4) {
      //if (i + 3 < listOfVertices.size()) {
        Vertex a, b, c, d;
        a = listOfVertices.at(i);
        b = listOfVertices.at(i+1);
        c = listOfVertices.at(i+2);
        d = listOfVertices.at(i+3);
        listOfTriangles.push_back(Triangle(a, b, c));
        listOfTriangles.push_back(Triangle(c, d, a));
        // }
    }
  }
  listOfVertices.clear();
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
  Vertex newVertex=(Vertex(currentColor[0],currentColor[1],currentColor[2],x,y,0.0f));
  transform(newVertex.position);
  listOfVertices.push_back(newVertex);//Vertex(x,y,INFINITY));
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
  Vertex newVertex=(Vertex(currentColor[0],currentColor[1],currentColor[2],x,y,z));
    transform(newVertex.position);
    listOfVertices.push_back(newVertex);//Vertex(x,y,z));
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
  currentStack().push_back({
                             1,0,0,0,
                             0,1,0,0,
                             0,0,1,0,
                             0,0,0,1
                           });
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
  for (MGLsize i = 0; i < 4; ++i) {
    for (MGLsize j = 0; j < 4; ++j) {
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
  for (MGLsize i = 0; i < 4; ++i) {
    for (MGLsize j = 0; j < 4; ++j) {
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
  MGLfloat A = (right + left) / (right - left);
  MGLfloat B = (top + bottom) / (top - bottom);
  MGLfloat C = -((far + near) / (far - near));
  MGLfloat D = -((2.0f*far*near) / (far - near));
  MGLfloat E = (2.0f*near) / (right - left);
  MGLfloat F = (2.0f*near) / (top - bottom);

  mat4 frustum = {
    E,0,0,0,
    0,F,0,0,
    A,B,C,-1,
    0,0,D,0
                 };
  mglMultMatrix(frustum.values);

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
  MGLfloat x = -((right + left) / (right - left));
  MGLfloat y = -((top + bottom) / (top - bottom));
  MGLfloat z = -((far + near) / (far - near));

  MGLfloat w = 2 / (right - left);
  MGLfloat v = 2 / (top - bottom);
  MGLfloat t = -2 / (far - near);

  mat4 orthoMat = {
    w,0,0,0,
    0,v,0,0,
    0,0,t,0,
    x,y,z,1
                  };
  mglMultMatrix(orthoMat.values);
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
