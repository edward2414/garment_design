//
//  Filename         : glutils.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Various GL related utility functions
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "glutils.h"
#include "math.h"

namespace glutils {

  void fadedLine2D(float a1,
		   float a2,
		   float b1,
		   float b2)
  {
    float k1 = (b1 - a1) / 4;
    float k2 = (b2 - a2) / 4;
    glBegin(GL_LINES);
    glColor4f(0.6, 0.6, 0.6, 1);
    glVertex2f(a1, a2);
    glColor4f(0.6, 0.6, 0.6, 0);
    glVertex2f(a1 + k1, a2 + k2);
    glVertex2f(b1 - k1, b2 - k2);
    glColor4f(0.6, 0.6, 0.6, 1);
    glVertex2f(b1, b2);
    glEnd();
  }

  void squareX(float a1,
	       float a2,
	       float b1,
	       float b2,
	       float x)
  {
    glBegin(GL_TRIANGLE_STRIP);
    glVertex3f(x, a1, a2);
    glVertex3f(x, a1, b2);
    glVertex3f(x, b1, a2);
    glVertex3f(x, b1, b2);
    glEnd();
  }
  
  void circle(float cx, float cy, float radius) 
  {
    float vectorY1=cy+radius;
    float vectorX1=cx;
    float vectorX=0.0;
    float vectorY=0.0;
    glBegin(GL_LINE_STRIP);         
    for(float angle=0.0f;angle<=(2.0f*M_PI);angle+=0.01f)
    {
      vectorX=cx+(radius*(float)sin(angle));
      vectorY=cy+(radius*(float)cos(angle));     
      glVertex2d(vectorX1,vectorY1);
      vectorY1=vectorY;
      vectorX1=vectorX;
    }
    glEnd();
  }

  void squareY(float a1,
	       float a2,
	       float b1,
	       float b2,
	       float y)
  {
    glBegin(GL_TRIANGLE_STRIP);
    glVertex3f(a1, y, a2);
    glVertex3f(a1, y, b2);
    glVertex3f(b1, y, a2);
    glVertex3f(b1, y, b2);
    glEnd();
  }

  void squareZ(float a1,
	       float a2,
	       float b1,
	       float b2,
	       float z)
  {
    glBegin(GL_TRIANGLE_STRIP);
    glVertex3f(a1, a2, z);
    glVertex3f(a1, b2, z);
    glVertex3f(b1, a2, z);
    glVertex3f(b1, b2, z);
    glEnd();
  }
  
  void fadedBox2D(float a1,
		  float a2,
		  float b1,
		  float b2)
  {
    fadedLine2D(a1, a2, a1, b2);
    fadedLine2D(a1, b2, b1, b2);
    fadedLine2D(b1, b2, b1, a2);
    fadedLine2D(b1, a2, a1, a2);
  }

  void line3D(float a1,
	      float a2,
	      float a3,
	      float b1,
	      float b2,
	      float b3) {
    glBegin(GL_LINES);
    glVertex3f(a1, a2, a3);
    glVertex3f(b1, b2, b3);
    glEnd();
  }

  void fadedLine3D(float a1,
		   float a2,
		   float a3,
		   float b1,
		   float b2,
		   float b3) {
    float k1 = (b1 - a1) / 4;
    float k2 = (b2 - a2) / 4;
    float k3 = (b3 - a3) / 4;
    glBegin(GL_LINES);
    glColor4f(0.6, 0.6, 0.6, 1);
    glVertex3f(a1, a2, a3);
    glColor4f(0.6, 0.6, 0.6, 0);
    glVertex3f(a1 + k1, a2 + k2, a3 + k3);
    glVertex3f(b1 - k1, b2 - k2, b3 - k3);
    glColor4f(0.6, 0.6, 0.6, 1);
    glVertex3f(b1, b2, b3);
    glEnd();
  }
  
  void box3D(float a1,
	     float a2,
	     float a3,
	     float b1,
	     float b2,
	     float b3) {
    line3D(a1, a2, a3, a1, b2, a3);
    line3D(a1, b2, a3, b1, b2, a3);
    line3D(b1, b2, a3, b1, a2, a3);
    line3D(b1, a2, a3, a1, a2, a3);
    line3D(a1, a2, b3, a1, b2, b3);
    line3D(a1, b2, b3, b1, b2, b3);
    line3D(b1, b2, b3, b1, a2, b3);
    line3D(b1, a2, b3, a1, a2, b3);
    line3D(a1, a2, a3, a1, a2, b3);
    line3D(a1, b2, a3, a1, b2, b3);
    line3D(b1, b2, a3, b1, b2, b3);
    line3D(b1, a2, a3, b1, a2, b3);
  }

  void fadedBox3D(float a1,
		  float a2,
		  float a3,
		  float b1,
		  float b2,
		  float b3) {
    fadedLine3D(a1, a2, a3, a1, b2, a3);
    fadedLine3D(a1, b2, a3, b1, b2, a3);
    fadedLine3D(b1, b2, a3, b1, a2, a3);
    fadedLine3D(b1, a2, a3, a1, a2, a3);
    fadedLine3D(a1, a2, b3, a1, b2, b3);
    fadedLine3D(a1, b2, b3, b1, b2, b3);
    fadedLine3D(b1, b2, b3, b1, a2, b3);
    fadedLine3D(b1, a2, b3, a1, a2, b3);
    fadedLine3D(a1, a2, a3, a1, a2, b3);
    fadedLine3D(a1, b2, a3, a1, b2, b3);
    fadedLine3D(b1, b2, a3, b1, b2, b3);
    fadedLine3D(b1, a2, a3, b1, a2, b3);
  }

  void cube(float a1,
	    float a2,
	    float a3,
	    float b1,
	    float b2,
	    float b3) {
    squareX(a2, a3, b2, b3, a1);
    squareX(b2, a3, a2, b3, b1);
    squareY(a1, b3, b1, a3, a2);
    squareY(a1, a3, b1, b3, b2);
    squareZ(a1, a2, b1, b2, a3);
    squareZ(a1, b2, b1, a2, b3);
  }
  
} // end of namespace glutils


