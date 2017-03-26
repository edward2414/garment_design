//
//  Filename         : vectypes.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Vectors instantiated types for easier use.
//  Date of creation : 12/06/2003
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  VECTYPES_H
# define VECTYPES_H

# include <vecmat.h>

// 2D.

typedef vecmat::Vec2<char>		Vec2c;
typedef vecmat::Vec2<unsigned char>	Vec2uc;
typedef vecmat::Vec2<short>		Vec2s;
typedef vecmat::Vec2<unsigned short>	Vec2us;
typedef vecmat::Vec2<int>		Vec2i;
typedef vecmat::Vec2<unsigned int>	Vec2ui;
typedef vecmat::Vec2<float>		Vec2f;
typedef vecmat::Vec2<double>		Vec2d;

// 3D.

typedef vecmat::Vec3<char>		Vec3c;
typedef vecmat::Vec3<unsigned char>	Vec3uc;
typedef vecmat::Vec3<short>		Vec3s;
typedef vecmat::Vec3<unsigned short>	Vec3us;
typedef vecmat::Vec3<int>		Vec3i;
typedef vecmat::Vec3<unsigned int>	Vec3ui;
typedef vecmat::Vec3<float>		Vec3f;
typedef vecmat::Vec3<double>		Vec3d;

typedef vecmat::HVec3<char>		HVec3c;
typedef vecmat::HVec3<unsigned char>	HVec3uc;
typedef vecmat::HVec3<short>		HVec3s;
typedef vecmat::HVec3<unsigned short>	HVec3us;
typedef vecmat::HVec3<int>		HVec3i;
typedef vecmat::HVec3<unsigned int>	HVec3ui;
typedef vecmat::HVec3<float>		HVec3f;
typedef vecmat::HVec3<double>		HVec3d;

#endif // VECTYPES_H
