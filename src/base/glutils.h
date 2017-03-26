//
//  Filename         : glutils.h
//  Author           : Emmanuel Turquin
//  Purpose          : Various GL related utility functions
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  GLUTILS_H
# define GLUTILS_H

# include <GL/gl.h>

namespace glutils {

  void squareX(float a1,
	       float a2,
	       float b1,
	       float b2,
	       float x);
  
  void squareY(float a1,
	       float a2,
	       float b1,
	       float b2,
	       float y);
  
  void squareZ(float a1,
	       float a2,
	       float b1,
	       float b2,
	       float z);
  
  void circle(float cx, float cy, float radius);

  void fadedLine2D(float a1,
		   float a2,
		   float b1,
		   float b2);
  
  void fadedBox2D(float a1,
		  float a2,
		  float b1,
		  float b2);

  void line3D(float a1,
	      float a2,
	      float a3,
	      float b1,
	      float b2,
	      float b3);

  void fadedLine3D(float a1,
		   float a2,
		   float a3,
		   float b1,
		   float b2,
		   float b3);
  
  void box3D(float a1,
	     float a2,
	     float a3,
	     float b1,
	     float b2,
	     float b3);

  void fadedBox3D(float a1,
		  float a2,
		  float a3,
		  float b1,
		  float b2,
		  float b3);

  void cube(float a1,
	    float a2,
	    float a3,
	    float b1,
	    float b2,
	    float b3);

} // end of namespace glutils

#endif // GLUTILS_H
