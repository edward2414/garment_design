//
// C++ Implementation: buckling
//
// Author: Philippe Decaudin
//

#include "buckling.h"
#include <cassert>
#include <QGLViewer/qglviewer.h>
#include <GL/gl.h>
#include <GL/glu.h>

using namespace std;

#define COL_TAG -7581	// for debug purpose
#define TOP_TAG -7582	// for debug purpose


//	Curved patch (parametric).
//	Subtraction between inner and outer control points define tangents.
struct CurvedPatch
{
	Vec3d	A[16];		// control points
						//	   A12-A13--A14-A15
						//      |   :    :   |
						//     A8..A9   A10.A11
						//      |            |
						//     A4..A5   A6..A7
						//      |   :    :   |
						//     A0--A1---A2--A3
	bool	inverseNormal;
	static inline void	evalCurve(Vec3d *p, const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d, double t)
	{
		double t2  = t*t;
		double t3  = t2*t;
		double it  = 1.0f-t;
		double it2 = it*it;
		double it3 = it2*it;
		*p = it3*a + (3.0f*it2*t)*b + (3.0f*it*t2)*c + t3*d;
	}
	static inline void	evalCurveTangent(Vec3d *tgt, const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d, double t)
	{
		double t2  = t*t;
		double it  = 1.0f-t;
		double it2 = it*it;
		*tgt = (3.0f*it2)*(b-a)+(3.0f*it*t)*(c-b)+(3.0f*t2)*(d-c);
	}
	inline void	eval(Vec3d *p, double u, double v) const
	{
		Vec3d a, b, c, d;
		evalCurve(&a, A[0] , A[1] , A[2] , A[3] , u);
		evalCurve(&b, A[4] , A[5] , A[6] , A[7] , u);
		evalCurve(&c, A[8] , A[9] , A[10], A[11], u);
		evalCurve(&d, A[12], A[13], A[14], A[15], u);
		evalCurve(p, a, b, c, d, v);
	}
	inline void	evalNormal(Vec3d *n, double u, double v) const
	{
		Vec3d a, b, c, d, tu, tv;
		evalCurveTangent(&a, A[0] , A[1] , A[2] , A[3] , u);
		evalCurveTangent(&b, A[4] , A[5] , A[6] , A[7] , u);
		evalCurveTangent(&c, A[8] , A[9] , A[10], A[11], u);
		evalCurveTangent(&d, A[12], A[13], A[14], A[15], u);
		evalCurve(&tu, a, b, c, d, v);
		evalCurveTangent(&a, A[0] , A[4] , A[8] , A[12], v);
		evalCurveTangent(&b, A[1] , A[5] , A[9] , A[13], v);
		evalCurveTangent(&c, A[2] , A[6] , A[10], A[14], v);
		evalCurveTangent(&d, A[3] , A[7] , A[11], A[15], v);
		evalCurve(&tv, a, b, c, d, u);
		if( inverseNormal )
			*n = tu ^ tv;
		else
			*n = tv ^ tu;
		n->normalizeSafe();
	}
	
	static inline void computeCP(Vec3d *out, const Vec3d& p1, const Vec3d& p2, const Vec3d& n1, const Vec3d& n2, const double tgtFactor=2.0/3.0)
	{
		(void)n2;
		Vec3d p12 = p2-p1;
		const double t  = tgtFactor;
		const double t1 = 1.0f-t;
		*out = t*p1 + t1*( p2 - (p12 * n1)*n1 );
	}	
	
	static inline double baryInterp(Vec3d *out, const Vec3d& p0, const Vec3d& pu, const Vec3d& pv, double u, double v)
	{
		double w = 1.0f-u-v;
		*out = w*p0 + u*pu + v*pv;
		return w;
	}
};


//	The 4 parts of a patch:
// 	          y
//	       <----+
//	6----7----8 |
//	|',  |  ,'| | x
//	|  ',|,'  | v
//	3----4----5
//	|  ,'|',  |
//	|,'  |  ',|
//	0----1----2
const int PARTS[4][4] = { {4,1,3,0},	// part0: bottom/left
						  {5,2,4,1},	// part1: bottom/right
						  {7,4,6,3},	// part2: top/left
						  {8,5,7,4} };	// part3: top/right

// Colors for patchs
const double COLORS[][3] = { {0.9,0.8,0.5}, {0.9,0.8,0.3}, {0.9,0.8,0.1}, {0.8,0.7,0.5}, {0.8,0.7,0.3}, {0.8,0.7,0.1}, {0.9,0.8,0.7} };

						  
ControlMesh::Patch::Patch()
{
	int i;
	for( i=0; i<9; ++i )
	{
		cpIndices[i] = -1;
		texCoords[i] = Vec2d(0, 0);
	}	
	textureID = -1;
	for( i=0; i<6; ++i )
		length0[i] = 0;
	for( i=0; i<4; ++i )
		_adjacencies[i] = -1;		
}


ControlMesh::ControlMesh()
{
	patchSubdiv = 4;				// 4
	patchTangentFactor = 0.6;		// 0.5 for a cylinder with 2 buckling patchs per section, else 2.0/3.0
	twistFactor = 3.0;				// 3.0
	diamondFactor = 0.2; 			// 0.3
	defaultCollisionMaxDistFactor = 0.1;
	verticalFoldsFactor = 1.0;
	twistVersusVerticalFactor = 100.0;
	diamondVersusVerticalFactor = 100.0;
    diamondBucklingEnabled = false; // causes non-intuitive square folding
	twistBucklingEnabled = true;
	verticalBucklingEnabled = true;
	smoothFilterEnabled = true;
	collisionEnabled = true;
	doCollisionsAtTheEnd = false;
	colorPatchsEnabled = true;
	defaultFrontColor = Vec3d(0.98,0.98,0.98);		//for video
	//defaultFrontColor = Vec3d(0.98,0.48,0.28);	//for figures
	defaultBackColor = 0.2*defaultFrontColor;
	textureEnabled = true;
	textureScale = Vec2d(1, 1);
	
    rows = 5;
    cols = 5;
	
	_bboxMin = Vec3d(-1,-1,-1);
	_bboxMax = Vec3d(1,1,1);
	_bboxRadius = (_bboxMax-_bboxMin).norm()/2;
	_moveOutsideCB = NULL;
	_moveOutsideClientData = NULL;
}


void ControlMesh::reset()
{
	patchs.resize(0);
	controlPoints.resize(0);
	cpSeamIDs.resize(0);
	_bboxMin = Vec3d(-1,-1,-1);
	_bboxMax = Vec3d(1,1,1);
	_bboxRadius = (_bboxMax-_bboxMin).norm()/2;
	_moveOutsideCB = NULL;
	_moveOutsideClientData = NULL;
}

 
int ControlMesh::addPatch(const Vec3d *cpArray, double cpInsideEps, const int *cpSeamIdArray, double cpSeamEps, const double *length0Array, const Vec2d *texCoordArray, int textureID, const int *cpAtlasIndices, const std::string atlasName, const Vec3d* interpPatternPoints)
{
	assert( cpArray!=NULL );
	Patch p;
    if(atlasName != NULL)
      p.atlasName = atlasName;	
	int i;
	for( i=0; i<9; ++i )	
	{
		if( cpSeamIdArray==NULL )
			p.cpIndices[i] = addControlPoint(cpArray[i], 0, cpInsideEps, p.cpIndices);
		else
			p.cpIndices[i] = addControlPoint(cpArray[i], cpSeamIdArray[i], (cpSeamIdArray[i]!=0) ? cpSeamEps : cpInsideEps, p.cpIndices);
			
	  	if( cpAtlasIndices!=NULL )
    		p.cpAtlasIndices[i] = cpAtlasIndices[i];	
    	
        if( interpPatternPoints!=NULL)
          p.interpPatternPoints[i] = interpPatternPoints[i];
	}
	
	if( length0Array!=NULL )
		for( i=0; i<6; ++i )
			p.length0[i] = length0Array[i];
	else
	{
		p.length0[0] = 	(cpArray[1]-cpArray[0]).norm() + (cpArray[2]-cpArray[1]).norm();
		p.length0[1] = 	(cpArray[7]-cpArray[6]).norm() + (cpArray[8]-cpArray[7]).norm();
		p.length0[2] = 	(cpArray[6]-cpArray[3]).norm() + (cpArray[3]-cpArray[0]).norm();
		p.length0[3] = 	(cpArray[8]-cpArray[5]).norm() + (cpArray[5]-cpArray[2]).norm();
		p.length0[4] = 	(cpArray[6]-cpArray[4]).norm() + (cpArray[4]-cpArray[2]).norm();
		p.length0[5] = 	(cpArray[8]-cpArray[4]).norm() + (cpArray[4]-cpArray[0]).norm();
	}
	if( texCoordArray!=NULL )
		for( i=0; i<9; ++i )
			p.texCoords[i] =  texCoordArray[i];
	p.textureID = textureID;
	int patchIdx = (int)patchs.size();
	patchs.push_back(p);
	return patchIdx;
}


int ControlMesh::addControlPoint(const Vec3d& cp, int seamID, double eps, int *currentPatchCPIndices)
{
	int i, j;
	int closestCP = -1;
	double dist, closestDist = eps;
	for( i=0; i<(int)controlPoints.size(); ++i )
	{
		bool belongsToSamePatch = false;
		if( currentPatchCPIndices!=NULL )	// a control point cannot be weld with a control point of the same patch
			for( j=0; j<9; ++j )
				if( currentPatchCPIndices[j]==i )
				{
					belongsToSamePatch = true;
					break;
				}
		dist = (controlPoints[i]-cp).norm();
		if( cpSeamIDs[i]==seamID && (!belongsToSamePatch) && dist<closestDist )
		{
			closestDist = dist;
			closestCP = i;
		}
	}
	if( closestCP>=0 )
		return closestCP;
	else
	{
		i = (int)controlPoints.size();
		controlPoints.push_back(cp);
		cpSeamIDs.push_back(seamID);
		assert( cpSeamIDs.size()==controlPoints.size() );
		return i;
	}
}


void ControlMesh::resetPatchsLength0()
{
	for( int i=0; i<(int)patchs.size(); ++i )
	{
		Patch& p = patchs[i];
		p.length0[0] = 	(controlPoints[p.cpIndices[1]]-controlPoints[p.cpIndices[0]]).norm() + (controlPoints[p.cpIndices[2]]-controlPoints[p.cpIndices[1]]).norm();
		p.length0[1] = 	(controlPoints[p.cpIndices[7]]-controlPoints[p.cpIndices[6]]).norm() + (controlPoints[p.cpIndices[8]]-controlPoints[p.cpIndices[7]]).norm();
		p.length0[2] = 	(controlPoints[p.cpIndices[6]]-controlPoints[p.cpIndices[3]]).norm() + (controlPoints[p.cpIndices[3]]-controlPoints[p.cpIndices[0]]).norm();
		p.length0[3] = 	(controlPoints[p.cpIndices[8]]-controlPoints[p.cpIndices[5]]).norm() + (controlPoints[p.cpIndices[5]]-controlPoints[p.cpIndices[2]]).norm();
		p.length0[4] = 	(controlPoints[p.cpIndices[6]]-controlPoints[p.cpIndices[4]]).norm() + (controlPoints[p.cpIndices[4]]-controlPoints[p.cpIndices[2]]).norm();
		p.length0[5] = 	(controlPoints[p.cpIndices[8]]-controlPoints[p.cpIndices[4]]).norm() + (controlPoints[p.cpIndices[4]]-controlPoints[p.cpIndices[0]]).norm();		
	}
}


int ControlMesh::getNumPatchs() const
{
	return (int)patchs.size();
}


int ControlMesh::getNumControlPoints() const
{
	return (int)controlPoints.size();
}


void ControlMesh::getPatch(int patchIndex, Vec3d *outCPArray, int *outCPIndicesArray, int *outCPSeamIdArray, double *outLength0Array, Vec2d *outTexCoordArray, int *outTextureID, int* outAtlasIndices, std::string* outAtlasName, Vec3d* outInterpPatternPoints) const
{
	assert( patchIndex>=0 && patchIndex<(int)patchs.size() );
	int i;
	if( outCPArray )
		for( i=0; i<9; ++i )
			outCPArray[i] = controlPoints[patchs[patchIndex].cpIndices[i]];
	if( outCPIndicesArray )
		for( i=0; i<9; ++i )
			outCPIndicesArray[i] = patchs[patchIndex].cpIndices[i];
	if( outCPSeamIdArray )
		for( i=0; i<9; ++i )
			outCPSeamIdArray[i] = cpSeamIDs[patchs[patchIndex].cpIndices[i]];
	if( outLength0Array )
		for( i=0; i<6; ++i )
			outLength0Array[i] = patchs[patchIndex].length0[i];
	if( outTexCoordArray )
		for( i=0; i<9; ++i )
			outTexCoordArray[i] = patchs[patchIndex].texCoords[i];
	if( outTextureID )
		*outTextureID = patchs[patchIndex].textureID;
    if(outAtlasIndices)
      for( i=0; i<9; ++i )
        outAtlasIndices[i] = patchs[patchIndex].cpAtlasIndices[i];
    if(outInterpPatternPoints)
      for( i=0; i<9; ++i )
        outInterpPatternPoints[i] = patchs[patchIndex].interpPatternPoints[i];
    if(outAtlasName)
      *outAtlasName = patchs[patchIndex].atlasName;		
}


const Vec3d& ControlMesh::getControlPoint(int cpIndex) const
{
	assert( cpIndex>=0 && cpIndex<(int)controlPoints.size() );
	return controlPoints[cpIndex];
}


void ControlMesh::setControlPoint(int cpIndex, const Vec3d& cp)
{
	assert( cpIndex>=0 && cpIndex<(int)controlPoints.size() );
	controlPoints[cpIndex] = cp;
}


double ControlMesh::getBBox(Vec3d *outBBoxMin, Vec3d *outBBoxMax) const
{
	const_cast<ControlMesh *>(this)->computeBBox();
	if( outBBoxMin!=NULL )
		*outBBoxMin = _bboxMin;
	if( outBBoxMax!=NULL )
		*outBBoxMax = _bboxMax;
	return _bboxRadius;
}


void ControlMesh::buildCylinder(int n, int m, const Vec3d& bbMin, const Vec3d& bbMax, double rBottom, double rTop, int textureID)
{
	Vec3d orig( (bbMin[0]+bbMax[0])/2, bbMin[1], (bbMin[2]+bbMax[2])/2 );
	double height = bbMax[1]-bbMin[1];
	double rx = (bbMax[0]-bbMin[0])/2;
	double rz = (bbMax[2]-bbMin[2])/2;
	if( rx>rz )
		rz = rx;
	else if( rz>rx )
		rx = rz;

	int i, j;
	Vec3d center, cp;
	double angle, rf;
	controlPoints.resize(0);
	
	const int NUM_RADII = 100;
	double radii[NUM_RADII];
	for( int rIdx=0; rIdx<NUM_RADII; ++rIdx )
		radii[rIdx] = rBottom + (double(rIdx)/NUM_RADII) * (rTop-rBottom);
	double *b = new double[2*m];
	computeBValues(b, 2*m, radii, NUM_RADII);
	
	center = orig;
	for( j=0; j<=m*2; ++j )
	{
		for( i=0; i<n*2; ++i )
		{
			angle = ((double)i)*(2.0f*M_PI)/(n*2);
			rf = rBottom + (double(j)/(m*2)) * (rTop-rBottom);
			cp = center + Vec3d( rx*rf*cos(angle), 0, rz*rf*sin(angle) );
			controlPoints.push_back(cp);
			cpSeamIDs.push_back(0);
		}
		if( j<m*2 )
			//center += Vec3d(0, double(height)/(m*2), 0);
			center += Vec3d(0, b[j]*double(height), 0);
	}
	assert( (int)controlPoints.size()==(n*2)*(m*2+1) );

	delete[] b;
	b = NULL;
	
	patchs.resize(0);
	Patch p;
	double u0, u1, u2, v0, v1, v2;
	for( j=0; j<m; ++j )
	{
		for( i=0; i<n; ++i )
		{
			p.cpIndices[0] = ((2*i+2)%(n*2)) + (2*j)*(n*2);		
			p.cpIndices[1] = 2*i+1           + (2*j)*(n*2);
			p.cpIndices[2] = 2*i             + (2*j)*(n*2);			

			p.cpIndices[3] = ((2*i+2)%(n*2)) + (2*j+1)*(n*2);						
			p.cpIndices[4] = 2*i+1           + (2*j+1)*(n*2);
			p.cpIndices[5] = 2*i             + (2*j+1)*(n*2);			
			
			p.cpIndices[6] = ((2*i+2)%(n*2)) + (2*j+2)*(n*2);			
			p.cpIndices[7] = 2*i+1           + (2*j+2)*(n*2);
			p.cpIndices[8] = 2*i             + (2*j+2)*(n*2);			
			
			p.length0[0] = 	(controlPoints[p.cpIndices[1]]-controlPoints[p.cpIndices[0]]).norm()
							+ (controlPoints[p.cpIndices[2]]-controlPoints[p.cpIndices[1]]).norm();
			p.length0[1] = 	(controlPoints[p.cpIndices[7]]-controlPoints[p.cpIndices[6]]).norm()
							+ (controlPoints[p.cpIndices[8]]-controlPoints[p.cpIndices[7]]).norm();
			p.length0[2] = 	(controlPoints[p.cpIndices[6]]-controlPoints[p.cpIndices[3]]).norm()
							+ (controlPoints[p.cpIndices[3]]-controlPoints[p.cpIndices[0]]).norm();
			p.length0[3] = 	(controlPoints[p.cpIndices[8]]-controlPoints[p.cpIndices[5]]).norm()
							+ (controlPoints[p.cpIndices[5]]-controlPoints[p.cpIndices[2]]).norm();
			p.length0[4] = 	(controlPoints[p.cpIndices[6]]-controlPoints[p.cpIndices[4]]).norm()
							+ (controlPoints[p.cpIndices[4]]-controlPoints[p.cpIndices[2]]).norm();
			p.length0[5] = 	(controlPoints[p.cpIndices[8]]-controlPoints[p.cpIndices[4]]).norm()
							+ (controlPoints[p.cpIndices[4]]-controlPoints[p.cpIndices[0]]).norm();
							
			u0 = 1.0-(double(i)+0.0)/n;
			u1 = 1.0-(double(i)+0.5)/n;
			u2 = 1.0-(double(i)+1.0)/n;
			v0 = 1.0-(double(j)+0.0)/m;
			v1 = 1.0-(double(j)+0.5)/m;
			v2 = 1.0-(double(j)+1.0)/m;
			p.texCoords[0] = Vec2d(u2, v0);
			p.texCoords[1] = Vec2d(u1, v0);
			p.texCoords[2] = Vec2d(u0, v0);
			p.texCoords[3] = Vec2d(u2, v1);
			p.texCoords[4] = Vec2d(u1, v1);
			p.texCoords[5] = Vec2d(u0, v1);
			p.texCoords[6] = Vec2d(u2, v2);
			p.texCoords[7] = Vec2d(u1, v2);
			p.texCoords[8] = Vec2d(u0, v2);
			p.textureID = textureID;
			
			patchs.push_back(p);
		}
	}
	assert( (int)patchs.size()==n*m );

}


void ControlMesh::computeNormals()
{
	size_t i, j;
	_cpNormals.resize(controlPoints.size());
	for( i=0; i<_cpNormals.size(); ++i )
		_cpNormals[i] = Vec3d(0,0,0);

	const int PATCH_SECTOR[8][3] = { {0,4,3}, {0,1,4}, {1,2,4}, {2,5,4}, {4,5,8}, {4,8,7}, {4,7,6}, {4,6,3} };
	Vec3d normal;
	for( j=0; j<patchs.size(); ++j )
	{
		Patch& p = patchs[j];
		for( i=0; i<8; ++i )
		{
			normal =   ( controlPoints[p.cpIndices[PATCH_SECTOR[i][2]]] - controlPoints[p.cpIndices[PATCH_SECTOR[i][0]]] ) 
					 ^ ( controlPoints[p.cpIndices[PATCH_SECTOR[i][0]]] - controlPoints[p.cpIndices[PATCH_SECTOR[i][1]]] );
			normal.normalizeSafe();
			_cpNormals[p.cpIndices[PATCH_SECTOR[i][0]]] += normal;
			_cpNormals[p.cpIndices[PATCH_SECTOR[i][1]]] += normal;
			_cpNormals[p.cpIndices[PATCH_SECTOR[i][2]]] += normal;
		}
	}

	for( i=0; i<_cpNormals.size(); ++i )
		_cpNormals[i].normalizeSafe();
}


void ControlMesh::computeAdjacencies()
{
	size_t i, j;

	for( i=0; i<patchs.size(); ++i )
		for( j=0; j<4; ++j )
			patchs[i]._adjacencies[j] = -1;
		
	for( i=0; i<patchs.size(); ++i )
	{
		int *indI = &(patchs[i].cpIndices[0]);
		for( j=0; j<patchs.size(); ++j )
			if( i!=j )
			{
				int *indJ = &(patchs[j].cpIndices[0]);			
				if( indI[0]==indJ[6] && indI[1]==indJ[7] && indI[2]==indJ[8] )
				{
					patchs[i]._adjacencies[0] = j;
					patchs[j]._adjacencies[1] = i;
				}
				else if( indI[6]==indJ[0] && indI[7]==indJ[1] && indI[8]==indJ[2] )
				{
					patchs[i]._adjacencies[1] = j;
					patchs[j]._adjacencies[0] = i;
				}
				else if( indI[0]==indJ[2] && indI[3]==indJ[5] && indI[6]==indJ[8] )
				{
					patchs[i]._adjacencies[2] = j;
					patchs[j]._adjacencies[3] = i;
				}	
				else if( indI[2]==indJ[0] && indI[5]==indJ[3] && indI[8]==indJ[6] )
				{
					patchs[i]._adjacencies[3] = j;
					patchs[j]._adjacencies[2] = i;
				}				
			}
	}
}


void ControlMesh::computeBBox()
{
	if( controlPoints.size()>0 )
	{
		_bboxMin = _bboxMax = controlPoints[0];
		for( size_t i=1; i<controlPoints.size(); ++i )
		{
			const Vec3d& cp = controlPoints[i];
			for( int j=0; j<3; ++j )
			{
				if( cp[j]<_bboxMin[j] ) 
					_bboxMin[j] = cp[j];
				if( cp[j]>_bboxMax[j] ) 
					_bboxMax[j] = cp[j];
			}
		}
		_bboxRadius = (_bboxMax-_bboxMin).norm()/2;
	}
	else
	{
		_bboxMin = _bboxMax = Vec3d(0,0,0);
		_bboxRadius = 0;
	}
}


static Vec3d rotate(const Vec3d& v, const Vec3d& axis, double angle)
{
	double axisLen = axis.norm();
	double x = axis[0]/axisLen;
	double y = axis[1]/axisLen;
	double z = axis[2]/axisLen;
	double c = cos(angle);
	double s = sin(angle);
	Vec3d r;
	r[0] = (x*x+c*(1.0-x*x))*v[0] + (x*y*(1.0-c)-z*s)*v[1] + (x*z*(1.0-c)+y*s)*v[2];
	r[1] = (x*y*(1.0-c)+z*s)*v[0] + (y*y+c*(1.0-y*y))*v[1] + (y*z*(1.0-c)-x*s)*v[2];
	r[2] = (x*z*(1.0-c)-y*s)*v[0] + (y*z*(1.0-c)+x*s)*v[1] + (z*z+c*(1.0-z*z))*v[2];
	return r;
}


void ControlMesh::moveControlPointsOutside()
{
	if( _moveOutsideCB==NULL )
		return;

	for( size_t i=1; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		Vec3d movedCP = _moveOutsideCB(cp, _moveOutsideClientData);
		if( (movedCP-cp).norm()>_bboxRadius )
			cpSeamIDs[i] = COL_TAG;
		else
			cp = movedCP;
	}
}


void ControlMesh::compressVert(double scale)
{
	computeBBox();
	Vec3d center = (_bboxMin+_bboxMax)/2;
	for( size_t i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		cp[1] = center[1] + scale*(cp[1]-center[1]);
	}
}


void ControlMesh::twistVert(double angleBottom, double angleTop)
{
std::cout<<"TWISTANGLE_IN="<<angleBottom<<std::endl;
	computeBBox();
	Vec3d center = (_bboxMin+_bboxMax)/2;	
	for( size_t i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		Vec3d c(center[0], cp[1], center[2]);
		Vec3d v = cp - c;
		double t = (cp[1]-_bboxMin[1])/(_bboxMax[1]-_bboxMin[1]);
		double cosa = cos((1.0-t)*angleBottom + t*angleTop);
		double sina = sin((1.0-t)*angleBottom + t*angleTop);
		Vec3d vrot = v;
		vrot[0] = cosa*v[0] - sina*v[2];
		vrot[2] = sina*v[0] + cosa*v[2];
		cp = c + vrot;
	}
}

 
void ControlMesh::compressHoriz(double scale0X, double scale0Z, double scale1X, double scale1Z)
{
	computeBBox();
	Vec3d center = (_bboxMin+_bboxMax)/2;
	for( size_t i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		double t = (cp[1]-_bboxMin[1])/(_bboxMax[1]-_bboxMin[1]);
		double sx = (1.0-t)*scale0X + t*scale1X;
		double sz = (1.0-t)*scale0Z + t*scale1Z;
		cp[0] = center[0] + sx*(cp[0]-center[0]);
		cp[2] = center[2] + sz*(cp[2]-center[2]);
	}
}


void ControlMesh::pushBackAlongAxis(const Vec3d& startAxis, const Vec3d& endAxis, double shrinkedLength, double startRigidLength, double endRigidLength, double influenceRadius)
{
	Vec3d axis = endAxis - startAxis;
	double origLength = axis.norm();
	axis /= origLength;
	Vec3d h;
	double s, fs, d;
	for( size_t i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		s = (cp-startAxis)*axis;
		h = startAxis + s*axis;
		if( (cp-h).norm()<influenceRadius && s<=origLength*1.0001 )
		{
			if( s<startRigidLength )
				fs = s;
			else if( s>origLength-endRigidLength )
				fs = s-origLength+shrinkedLength;
			else
			{
				d = (s-startRigidLength)/(origLength-endRigidLength-startRigidLength);
				fs = (1.0-d)*s + d*(s-origLength+shrinkedLength);
			}
			cp = startAxis + fs*axis + (cp-h);
		}
	}
}


void ControlMesh::twistAlongAxis(const Vec3d& startAxis, const Vec3d& endAxis, double startAngle, double endAngle, double influenceRadius)
{
	Vec3d axis = endAxis - startAxis;
	double axisLen = axis.norm();
	axis /= axisLen;
	Vec3d h, v;
	double s;
	for( size_t i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		s = (cp-startAxis)*axis;
		h = startAxis + s*axis;
		v = cp - h;
		if( v.norm()<influenceRadius && s>=0 && s<=axisLen )
			cp = h + rotate(v, axis, (1.0-s/axisLen)*startAngle + (s/axisLen)*endAngle);
	}
}


void ControlMesh::scaleAlongAxis(const Vec3d& startAxis, const Vec3d& endAxis, double startScale, double endScale, double influenceRadius)
{
	Vec3d axis = endAxis - startAxis;
	double axisLen = axis.norm();
	axis /= axisLen;
	Vec3d h, v;
	double s;
	for( size_t i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		s = (cp-startAxis)*axis;
		h = startAxis + s*axis;
		v = cp - h;
		if( v.norm()<influenceRadius && s>=0 && s<=axisLen )
			cp = h + ((1.0-s/axisLen)*startScale + (s/axisLen)*endScale)*v;
	}
}

/*
void ControlMesh::bendBone(const Vec3d& pivot, const Vec3d& endBone, const Vec3d& rotAxis, double rotAngle, double blendLength, double influenceRadius)
{
	Vec3d bone = pivot - endBone;
	double boneLen = bone.norm();
	bone /= boneLen;
	Vec3d w = (bone ^ rotAxis).normalizeSafe();
	Vec3d h, v;
	double s, t, fs, a;
	for( size_t i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		s = (cp-pivot)*bone;
		h = pivot + s*bone;
		v = cp - h;
		a = - ((cp-pivot)*w)/((cp-pivot).norm()) * (rotAngle/(4.0*M_PI));		
		//a = - ((cp-pivot)*w)/((cp-pivot).norm()) * (M_PI/2.0/(M_PI/8.0)); //(rotAngle/(4.0*M_PI));
		if( a>0 )
			a = 0;
		if( v.norm()<influenceRadius )
		{
			if( s>-blendLength && s<blendLength )
			{
				// stretch/shrink
				t = s/blendLength;
				fs = s + a * s * (1.0+cos(M_PI*t))/2.0;
				cp = pivot + fs*bone + (cp-h);
		
				// bend
				const double BEND_FACTOR = 0.2;
				t = (s+BEND_FACTOR*blendLength)/(2.0*(BEND_FACTOR*blendLength));
				if( t>0 && t<1 )
					cp = (1.0-t)*cp + t*(pivot + rotate(cp-pivot ,rotAxis, rotAngle));
				else if( t>0 )
					cp = pivot + rotate(cp-pivot ,rotAxis, rotAngle);
			}
			else if( s>0 && s<=boneLen*1.0001 )
				cp = pivot + rotate(cp-pivot, rotAxis, rotAngle);
		}
	}
}
*/


void ControlMesh::bendBone(const Vec3d& pivot, const Vec3d& endBone, const Vec3d& rotAxis, double rotAngle, double blendLength, double roundPow, double influenceRadius)
{
	Vec3d eu = endBone - pivot;
	double boneLen = eu.norm();
	eu /= boneLen;
	Vec3d ev = rotAxis/rotAxis.norm();
	Vec3d ew = (ev ^ eu).normalizeSafe();
	
	double u, v, w, t;
	Vec3d euRot = rotate(eu, rotAxis, rotAngle);
	Vec3d ewRot = rotate(ew, rotAxis, rotAngle);
	Vec3d ewRot2 = rotate(ew, rotAxis, rotAngle/2);
	double k, ewRot2Coeff = 1.0/cos(rotAngle/2);
	for( size_t i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		u = (cp-pivot)*eu;
		v = (cp-pivot)*ev;
		w = (cp-pivot)*ew;
		if( sqrt(v*v+w*w)<influenceRadius && u<=1.0001*boneLen && u>-blendLength )
		{
			if( u>=blendLength )
				cp = pivot + rotate(cp-pivot, rotAxis, rotAngle);
			else if( u<=0 )
			{
				t = -u/blendLength;
				if( w>0 )
					k = (1.0-pow(t,roundPow))*1.0 + ewRot2Coeff*pow(t,roundPow);
				else
					k = ewRot2Coeff;
				cp = pivot + u*eu + v*ev + w*( (1.0-t)*k*ewRot2 + t*ew );	
			}
			else
			{
				t = u/blendLength;
				if( w>0 )
					k = (1.0-pow(t,roundPow))*1.0 + ewRot2Coeff*pow(t,roundPow);
				else
					k = ewRot2Coeff;
				cp = pivot + u*euRot + v*ev + w*( (1.0-t)*k*ewRot2 + t*ewRot );
			}
		}
	}
}


void ControlMesh::jitterAlongAxis(const Vec3d& startAxis, const Vec3d& endAxis, double axisRandom, double radiusRandom, double angleRandom, double influenceRadius)
{
	srand(999);
	Vec3d axis = endAxis - startAxis;
	double axisLen = axis.norm();
	axis /= axisLen;
	Vec3d h, v;
	double s;
	vector<Vec3d> jitters;
	jitters.resize(controlPoints.size());
	size_t i;	
	for( i=0; i<controlPoints.size(); ++i )
	{
		jitters[i][0] = axisRandom*(double(rand())/RAND_MAX-0.5)*axisLen/4;
		jitters[i][1] = radiusRandom*(double(rand())/RAND_MAX-0.5)/2;
		jitters[i][2] = angleRandom*(double(rand())/RAND_MAX-0.5)*M_PI/2;
	}	
	for( i=0; i<patchs.size(); ++i )
	{
		const Patch& p = patchs[i];
		jitters[p.cpIndices[1]] = 0.5*(jitters[p.cpIndices[0]] + jitters[p.cpIndices[2]]);
		jitters[p.cpIndices[3]] = 0.5*(jitters[p.cpIndices[0]] + jitters[p.cpIndices[6]]);
		jitters[p.cpIndices[5]] = 0.5*(jitters[p.cpIndices[2]] + jitters[p.cpIndices[8]]);
		jitters[p.cpIndices[7]] = 0.5*(jitters[p.cpIndices[6]] + jitters[p.cpIndices[8]]);
		jitters[p.cpIndices[4]] = 0.25*(jitters[p.cpIndices[1]] + jitters[p.cpIndices[3]] + jitters[p.cpIndices[5]] + jitters[p.cpIndices[7]]);
	}
	for( i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		s = (cp-startAxis)*axis;
		h = startAxis + s*axis;
		v = cp - h;
		if( v.norm()<influenceRadius && s>=0 && s<=axisLen )
		{
			h += jitters[i][0]*axis;
			v *= 1.0 + jitters[i][1];
			cp = h + rotate(v, axis, jitters[i][3]);
		}
	}
	/*
	for( i=0; i<controlPoints.size(); ++i )
	{
		Vec3d& cp = controlPoints[i];
		s = (cp-startAxis)*axis;
		h = startAxis + s*axis;
		v = cp - h;
		if( v.norm()<influenceRadius && s>=0 && s<=axisLen )
		{
			h += (axisRandom*(double(rand())/RAND_MAX-0.5)*axisLen/4)*axis;
			v *= 1.0 + radiusRandom*(double(rand())/RAND_MAX-0.5)/2;
			disp[i] = ( h + rotate(v, axis, angleRandom*(double(rand())/RAND_MAX-0.5)*M_PI/2) ) - cp;
		}
	}
	*/	
}


void ControlMesh::moveTopControlPoints(double distMin, DistanceToColliderCallback distCB)
{
	if( distCB==NULL )
		return;
	computeBBox();
	computeNormals();
	assert( _cpNormals.size()==controlPoints.size() );
	assert( _cpNormals.size()==cpSeamIDs.size() );
	size_t i;
	const double COS_MIN = cos(0.3*M_PI);
	// detect "top" control points
	for( i=0; i<controlPoints.size(); ++i )
	{
		if( _cpNormals[i][1]>=COS_MIN )
		{
			Vec3d p = controlPoints[i];
			const double NO_DIST = -999;
			double dist = NO_DIST;
			distCB(p[0], p[1], p[2], dist);
			if( dist!=NO_DIST && dist>distMin && dist<0.2*(_bboxMax[1]-_bboxMin[1]) )
			{
				cpSeamIDs[i] = TOP_TAG;
				
				p[1] -= dist;
				double newDist = NO_DIST;
				distCB(p[0], p[1], p[2], newDist);
				if( newDist!=NO_DIST && newDist<dist )
				{
					controlPoints[i] = p;
				}
			}
		}
	}
}


void ControlMesh::update(MoveOutsideCallback moveOutsideCB, void *clientData)
{
	_moveOutsideCB = moveOutsideCB;
	_moveOutsideClientData = clientData;
	if( collisionEnabled ) // && !doCollisionsAtTheEnd ) --> no! control points should always be moved first
		moveControlPointsOutside();
	computeBBox();
	computeNormals();
	computeAdjacencies();
	initRefinedMesh();
}


void ControlMesh::initRefinedMesh()
{
	_refinedVertices.resize(0);
	_refinedNormals.resize(0);

	size_t patchIndex;
	int x, y;
	double tx, ty;
	Vec3d q;
		
	for( patchIndex=0; patchIndex<patchs.size(); ++patchIndex )
	{
		Patch& p = patchs[patchIndex];
		for( int part=0; part<4; ++part )
		{
			p._vtxIndices[part].resize((patchSubdiv+1)*(patchSubdiv+1));
			const Vec3d& p0 = controlPoints[p.cpIndices[PARTS[part][0]]];
			const Vec3d& p1 = controlPoints[p.cpIndices[PARTS[part][1]]];
			const Vec3d& p2 = controlPoints[p.cpIndices[PARTS[part][2]]];
			const Vec3d& p3 = controlPoints[p.cpIndices[PARTS[part][3]]];
		
			for( y=0; y<=patchSubdiv; ++y )			
			{
				ty = double(y)/patchSubdiv;	
				for( x=0; x<=patchSubdiv; ++x )
				{
					tx = double(x)/patchSubdiv;	
					q = (1.0-ty)*( (1.0-tx)*p0+tx*p1 ) + ty*( (1.0-tx)*p2+tx*p3 );
					p._vtxIndices[part][x+y*(patchSubdiv+1)] = addRefinedVertex(q, Vec3d(0,0,0));
				}
			}
		}
	}
	 
	computeRefinedAdjacencies();
}


int ControlMesh::addRefinedVertex(const Vec3d& vertex, const Vec3d& normal)
{
	const double EPS = 1.0e-12;
	size_t i;
	for( i=0; i<_refinedVertices.size(); ++i )
	{
		const Vec3d& v = _refinedVertices[i];
		if( fabs(v[0]-vertex[0])<EPS && fabs(v[1]-vertex[1])<EPS && fabs(v[2]-vertex[2])<EPS )
			return (int)i;
	}
	i = _refinedVertices.size();
	_refinedVertices.push_back(vertex);
	_refinedNormals.push_back(normal);
	assert( _refinedVertices.size()==_refinedNormals.size() );
	return (int)i;
}


void ControlMesh::VtxAdj::append(int table, int i)
{
	int j;
	for( j=0; j<count[table]; ++j )
		if( indices[table][j]==i )
			return;	// already done
	//assert( Count[table]<MAX );
	if( count[table]<MAX )
		indices[table][count[table]++] = i;
}


void ControlMesh::computeRefinedAdjacencies()
{
	int x, y;
	int part, a, b, c, d;
	_refinedAdjacencies.clear();
	_refinedAdjacencies.resize(_refinedVertices.size());
	
	for( size_t patchIndex=0; patchIndex<patchs.size(); ++patchIndex )
	{
		const Patch& patch = patchs[patchIndex];
		for( part=0; part<4; ++part )
		{
			for( y=0; y<patchSubdiv; ++y )
			{
				for( x=0; x<patchSubdiv; ++x )
				{
					a = patch._vtxIndices[part][(x  )+(y  )*(patchSubdiv+1)];
					b = patch._vtxIndices[part][(x+1)+(y  )*(patchSubdiv+1)];
					c = patch._vtxIndices[part][(x  )+(y+1)*(patchSubdiv+1)];
					d = patch._vtxIndices[part][(x+1)+(y+1)*(patchSubdiv+1)];

					_refinedAdjacencies[a].append(0, b);
					_refinedAdjacencies[a].append(0, c);
					_refinedAdjacencies[a].append(1, d);
					
					_refinedAdjacencies[b].append(0, a);
					_refinedAdjacencies[b].append(0, d);
					_refinedAdjacencies[b].append(1, c);
					
					_refinedAdjacencies[c].append(0, a);
					_refinedAdjacencies[c].append(0, d);
					_refinedAdjacencies[c].append(1, b);
					
					_refinedAdjacencies[d].append(0, c);
					_refinedAdjacencies[d].append(0, b);
					_refinedAdjacencies[d].append(1, a);					
				}
			}
		}
	}
}


void ControlMesh::smoothRefinedMesh()
{
	// This method applies a square filter to the refined mesh
	// Coeffs of the filter are (1/4, 1/8, 1/16).

	size_t i, j;
	
	static vector<Vec3d> s_CopyOfRefinedVertices;	// static to improve perf
	s_CopyOfRefinedVertices.resize(_refinedVertices.size());
	for( i=0; i<_refinedVertices.size(); ++i )
		s_CopyOfRefinedVertices[i] = _refinedVertices[i];
		
	for( i=0; i<_refinedVertices.size(); ++i )
	{
		const VtxAdj& adj = _refinedAdjacencies[i];
		Vec3d& v = _refinedVertices[i];
		if( adj.count[0]==4 && adj.count[1]==4 )
		{
			v = 0.25 * s_CopyOfRefinedVertices[i];
			for( j=0; j<4; ++j )
			{
				v += 0.125  * s_CopyOfRefinedVertices[ adj.indices[0][j] ];
				v += 0.0625 * s_CopyOfRefinedVertices[ adj.indices[1][j] ];
			}
		}
		else if ( adj.count[0]==3 )
		{
			v = 0.5 * s_CopyOfRefinedVertices[i];
			for( j=0; j<3; ++j )
				v += (1.0/6.0)  * s_CopyOfRefinedVertices[ adj.indices[0][j] ];
		}
	}	
}


void ControlMesh::drawControlMesh(bool drawControlPoints, bool drawNormals, bool drawAdjacencies)
{
	cout << "DRAW CONTROL MESH: " << patchs.size() << " patchs, " << controlPoints.size() << " control points" << endl;

	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	
	size_t patchIndex, j;
	for( patchIndex=0; patchIndex<patchs.size(); ++patchIndex )
	{
		const Patch& p = patchs[patchIndex];
				
  
		// Draw control points
		//if( drawControlPoints )
        if(0)
		{
			glPointSize(8.0f);
			glBegin(GL_POINTS);	
			for( j=0; j<9; ++j )
			{
				if( cpSeamIDs[p.cpIndices[j]]==COL_TAG )
					glColor3d(1, 1, 1);
				else if( cpSeamIDs[p.cpIndices[j]]==TOP_TAG )
					glColor3d(1, 0, 1);					
				else if( cpSeamIDs[p.cpIndices[j]]!=0 )
					glColor3d(1, 0, 0);
				else
					glColor3d(0, 0.9, 0.9);
				glVertex3dv( controlPoints[p.cpIndices[j]].address() );
			}
			glEnd();	
		}
	
		// Draw edges
		const int EDGES[8][3] = { {0,1,2}, {3,4,5}, {6,7,8}, {0,3,6}, {1,4,7}, {2,5,8}, {2,4,6}, {0,4,8} };	
        glLineWidth(4.0);
		glBegin(GL_LINES);
		for( j=0; j<8; ++j )
		{
    
            if(j==4) continue;
			if( cpSeamIDs[p.cpIndices[EDGES[j][0]]]!=0 && cpSeamIDs[p.cpIndices[EDGES[j][1]]]!=0 && cpSeamIDs[p.cpIndices[EDGES[j][2]]]!=0 )
				//glColor3d(0.8, 0, 0);
                glColor3d(0.7, 0.3, 0.3);
			else if( j>=3 && j<=5 )
				//glColor3d(0.5, 0.5, 0.5);
     
                 glColor3d(0.0,0.0,0.7); // blue verticals
			else if( j>=6 )
				//glColor3d(0, 0, 0.8);
                  glColor3d(0.0,0.7,0.0); // green diagonals
			else
				glColor3d(0.7, 0.7, 0.7);
			double eps = patchs[0].length0[0]/200;
			Vec3d a = controlPoints[p.cpIndices[EDGES[j][0]]] + eps*_cpNormals[p.cpIndices[EDGES[j][0]]];
			glVertex3dv( a.address() );
			a = controlPoints[p.cpIndices[EDGES[j][1]]] + eps*_cpNormals[p.cpIndices[EDGES[j][1]]];			
			glVertex3dv( a.address() );
			glVertex3dv( a.address() );
			a = controlPoints[p.cpIndices[EDGES[j][2]]] + eps*_cpNormals[p.cpIndices[EDGES[j][2]]];
			glVertex3dv( controlPoints[p.cpIndices[EDGES[j][2]]].address() );	
   
		}
		glEnd();
  glLineWidth(1.0);
		
		// Draw faces
		glEnable(GL_LIGHTING);
		if( colorPatchsEnabled )
			glColor3dv(COLORS[patchIndex%(sizeof(COLORS)/sizeof(COLORS[0]))]);
		else
			glColor3d(0.98, 0.98, 0.98);
		if( textureEnabled && p.textureID>=0 )
		{
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, p.textureID);
			glMatrixMode(GL_TEXTURE);
			glLoadIdentity();
			glScaled(textureScale[0], textureScale[1], 1);
			glMatrixMode(GL_MODELVIEW);
		}
		const int FACES[8][3] = { {0,4,3}, {0,1,4}, {1,2,4}, {2,5,4}, {4,5,8}, {4,8,7}, {4,7,6}, {4,6,3} };
		glBegin(GL_TRIANGLES);
		for( j=0; j<8; ++j )
		{
			glTexCoord2dv( p.texCoords[FACES[j][0]].address() );
			glNormal3dv( _cpNormals[p.cpIndices[FACES[j][0]]].address() );		
			glVertex3dv( controlPoints[p.cpIndices[FACES[j][0]]].address() );
			glTexCoord2dv( p.texCoords[FACES[j][1]].address() );			
			glNormal3dv( _cpNormals[p.cpIndices[FACES[j][1]]].address() );			
			glVertex3dv( controlPoints[p.cpIndices[FACES[j][1]]].address() );
			glTexCoord2dv( p.texCoords[FACES[j][2]].address() );			
			glNormal3dv( _cpNormals[p.cpIndices[FACES[j][2]]].address() );			
			glVertex3dv( controlPoints[p.cpIndices[FACES[j][2]]].address() );
		}
		glEnd();		
		glDisable(GL_TEXTURE_2D);		
		glDisable(GL_LIGHTING);
				
		// Draw normals
		if( drawNormals )
		{
			glColor3f(0, 0.9f, 0);
			glBegin(GL_LINES);
			for( j=0; j<9; ++j )
			{
				const Vec3d& a = controlPoints[p.cpIndices[j]];
				Vec3d b = a + (p.length0[0]/6) * _cpNormals[p.cpIndices[j]];
				glVertex3dv( a.address() );
				glVertex3dv( b.address() );
			}
			glEnd();
		}
		
		// Draw adjacencies
		if( drawAdjacencies )
		{
			const int MID_A[4] = { 1, 7, 3, 5 };
			const int MID_B[4] = { 7, 1, 5, 3 };	
			glColor3f(0.9f, 0, 0.9f);
			glBegin(GL_LINES);
			for( j=0; j<4; ++j )
				if( p._adjacencies[j]>=0 )
				{
					double eps = patchs[0].length0[0]/100;
					Vec3d a = ( controlPoints[p.cpIndices[MID_A[j]]] + controlPoints[p.cpIndices[4]] )/2 + eps*_cpNormals[p.cpIndices[4]];
					const Patch& q = patchs[p._adjacencies[j]];
					Vec3d b = ( controlPoints[q.cpIndices[MID_B[j]]] + controlPoints[q.cpIndices[4]] )/2 + eps*_cpNormals[q.cpIndices[4]];
					glVertex3dv( a.address() );
					glVertex3dv( b.address() );
				}
			glEnd();
		}
	}
	
	glPopAttrib();	
}


void ControlMesh::drawCurvedSurface()
{
	size_t i;
	for( i=0; i<patchs.size(); ++i )
		updateRefinedCurvedPatch(i);	
	
	for( i=0; i<patchs.size(); ++i )
		drawRefinedPatch(i);
}


void ControlMesh::updateRefinedCurvedPatch(int patchIndex)
{
	const Patch& patch = patchs[patchIndex];
	
	//	4 curved patchs are displayed for each patch
	//		+---+---+
	//		| 2 | 3 |
	//		+---X---+
	//		| 0 | 1 |
	//		+---+---+
	int x, y, borderX, borderY, parseBorderX, parseBorderY;
	double tx, ty;
	Vec3d onCurve, onCurveNorm;
	CurvedPatch curvedPatch;
	curvedPatch.inverseNormal = false;
	int part, vtxIndex;
	for( part=0; part<4; ++part )
	{
		const Vec3d& n0 = _cpNormals[patch.cpIndices[PARTS[part][0]]];
		const Vec3d& n1 = _cpNormals[patch.cpIndices[PARTS[part][1]]];
		const Vec3d& n2 = _cpNormals[patch.cpIndices[PARTS[part][2]]];
		const Vec3d& n3 = _cpNormals[patch.cpIndices[PARTS[part][3]]];
		const Vec3d& p0 = controlPoints[patch.cpIndices[PARTS[part][0]]];
		const Vec3d& p1 = controlPoints[patch.cpIndices[PARTS[part][1]]];
		const Vec3d& p2 = controlPoints[patch.cpIndices[PARTS[part][2]]];
		const Vec3d& p3 = controlPoints[patch.cpIndices[PARTS[part][3]]];

		curvedPatch.A[0] = p0;
		CurvedPatch::computeCP(&curvedPatch.A[1], p0, p1, n0, n1, patchTangentFactor);
		CurvedPatch::computeCP(&curvedPatch.A[4], p0, p2, n0, n2, patchTangentFactor);
		curvedPatch.A[5] = curvedPatch.A[1] + curvedPatch.A[4] - p0;

		curvedPatch.A[3] = p1;
		CurvedPatch::computeCP(&curvedPatch.A[2], p1, p0, n1, n0, patchTangentFactor);
		CurvedPatch::computeCP(&curvedPatch.A[7], p1, p3, n1, n3, patchTangentFactor); 
		curvedPatch.A[6] = curvedPatch.A[2] + curvedPatch.A[7] - p1;

		curvedPatch.A[12] = p2;
		CurvedPatch::computeCP(&curvedPatch.A[8], p2, p0, n2, n0, patchTangentFactor);
		CurvedPatch::computeCP(&curvedPatch.A[13], p2, p3, n2, n3, patchTangentFactor); 
		curvedPatch.A[9] = curvedPatch.A[8] + curvedPatch.A[13] - p2;

		curvedPatch.A[15] = p3;
		CurvedPatch::computeCP(&curvedPatch.A[11], p3, p1, n3, n1, patchTangentFactor);
		CurvedPatch::computeCP(&curvedPatch.A[14], p3, p2, n3, n2, patchTangentFactor); 
		curvedPatch.A[10] = curvedPatch.A[11] + curvedPatch.A[14] - p3;

		borderY = 1; //( i==0 || i==2 ) ? 1 : 0;
		for( y=0; y<patchSubdiv + borderY; ++y )
		{
			parseBorderY = ( y==patchSubdiv ) ? 1 : 0;
			ty = double(y)/patchSubdiv;
			borderX = 1; //( i==0 || i==1 ) ? 1 : 0;
			for( x=0; x<patchSubdiv + borderX; ++x )
			{
				parseBorderX = ( x==patchSubdiv ) ? 1 : 0;
				tx = double(x)/patchSubdiv;

				//	point on the curved surface
				curvedPatch.eval(&onCurve, tx, ty);
				curvedPatch.evalNormal(&onCurveNorm, tx, ty);
								
				vtxIndex = patch._vtxIndices[part][x+y*(patchSubdiv+1)];
				_refinedVertices[vtxIndex] = onCurve;
				_refinedNormals[vtxIndex] = onCurveNorm;
			}
		}
	}
}


void ControlMesh::drawBucklingSurface()
{
	cout << "DRAW BUCKLING SURFACE: " << patchs.size() << " patchs, " << controlPoints.size() << " control points, ";
	cout << _refinedVertices.size() << " refined vertices, ";
	cout << patchs.size()*4*(patchSubdiv*patchSubdiv*2) << " triangles." << endl;

	size_t i;
	for( i=0; i<patchs.size(); ++i )
		updateRefinedBucklingPatch(i);	
	
	if( smoothFilterEnabled )
		smoothRefinedMesh();
		
	if( collisionEnabled && doCollisionsAtTheEnd )
		moveRefinedVerticesOutside();
		
	updateRefinedNormals();
		
	for( i=0; i<patchs.size(); ++i )
		drawRefinedPatch(i);
}


void ControlMesh::updateRefinedBucklingPatch(int patchIndex)
{
	const Patch& patch = patchs[patchIndex];

	// control points positions
	Vec3d pos[9];
	for( int ii=0; ii<9; ++ii )
		pos[ii] = controlPoints[patch.cpIndices[ii]];
		
	// external control points
	Vec3d extPosLeft, extPosRight;
	if( patch._adjacencies[2]>=0 )	// 2=left
		extPosLeft = controlPoints[patchs[patch._adjacencies[2]].cpIndices[4]];
	else
		extPosLeft = pos[3];
	if( patch._adjacencies[3]>=0 )	// 3=right
		extPosRight = controlPoints[patchs[patch._adjacencies[3]].cpIndices[4]];
	else
		extPosRight = pos[3];
		

	//	0. Coumpute buckling factors
	
	double compressLeft  = 0;
	double compressRight = 0;
	double compressVert = 0;
	double length0[6];
	Vec3d vl0, vl1;
	vl0 = pos[0]-pos[1];
	vl1 = pos[1]-pos[2];
	length0[0] = vl0.norm() + vl1.norm();	// bottom
	vl0 = pos[6]-pos[7];
	vl1 = pos[7]-pos[8];
	length0[1] = vl0.norm() + vl1.norm();	// top
	vl0 = pos[0]-pos[3];
	vl1 = pos[3]-pos[6];
	length0[2] = vl0.norm() + vl1.norm();	// left
	vl0 = pos[2]-pos[5];
	vl1 = pos[5]-pos[8];
	length0[3] = vl0.norm() + vl1.norm();	// right
	vl0 = pos[2]-pos[4];
	vl1 = pos[4]-pos[6];
	length0[4] = vl0.norm() + vl1.norm();	// diag topleft/bottomright
	vl0 = pos[0]-pos[4];
	vl1 = pos[4]-pos[8];
	length0[5] = vl0.norm() + vl1.norm();	// diag bottomleft/topright
	compressLeft = 1.0f - length0[2] / patch.length0[2];
	if( compressLeft<0 )
		compressLeft = 0;
	compressRight = 1.0f - length0[3] / patch.length0[3];
	if( compressRight<0 )
		compressRight = 0;
	double compressLeftDiff = patch.length0[2] - length0[2];
	if( compressLeftDiff<0 )
		compressLeftDiff = 0;
	double compressRightDiff = patch.length0[3] - length0[3];
	if( compressRightDiff<0 )
		compressRightDiff = 0;
	double compressTop = 1.0f - length0[1] / patch.length0[1];
	if( compressTop<0 )
		compressTop = 0;	
	double compressTopDiff = patch.length0[1] - length0[1];
	if( compressTopDiff<0 )
		compressTopDiff = 0;
	double compressBottom = 1.0f - length0[0] / patch.length0[0];
	if( compressBottom<0 )
		compressBottom = 0;			
	double compressBottomDiff = patch.length0[0] - length0[0];
	if( compressBottomDiff<0 )
		compressBottomDiff = 0;
	double twist4 = 0;
	double twist5 = 0;
	double twist;
	twist4 = 1.0f - length0[4] / patch.length0[4];
	if( twist4<0 )
		twist4 = 0;
	twist5 = 1.0f - length0[5] / patch.length0[5];
	if( twist5<0 )
		twist5 = 0;
	twist = twist5 - twist4;
	//if( twist4>twist5 )
	//	twist = - (twist4 - twist5);
	//else
	//	twist = twist5 - twist4;
	double twistVersusVerticalFolds;
	if( twistBucklingEnabled && verticalBucklingEnabled )
		twistVersusVerticalFolds = 1.0/(1.0 + twistVersusVerticalFactor*(compressTop+compressBottom));
	else if( twistBucklingEnabled )
		twistVersusVerticalFolds = 1.0;
	else
		twistVersusVerticalFolds = 0.0;
	twist *= twistVersusVerticalFolds;
	
	double diamondVersusVerticalFolds;
	if( diamondBucklingEnabled && verticalBucklingEnabled )
		diamondVersusVerticalFolds = 1.0/(1.0 + diamondVersusVerticalFactor*(compressTop+compressBottom));
	else if( diamondBucklingEnabled )
		diamondVersusVerticalFolds = 1.0;
	else
		diamondVersusVerticalFolds = 0.0;
	compressRight *= diamondVersusVerticalFolds;
	compressLeft *= diamondVersusVerticalFolds;
	
	
	//	1. Compute "pushed" positions: middle control points are moved onto the diamond triangles
	
	Vec3d pushedPos[9];
	pushedPos[0] = pos[0];
	pushedPos[1] = 0.5*(pos[0]+pos[2]);
	pushedPos[2] = pos[2];
	pushedPos[3] = 0.5*(extPosLeft+pos[4]);
	pushedPos[4] = pos[4];
	pushedPos[5] = 0.5*(extPosRight+pos[4]);
	pushedPos[6] = pos[6];
	pushedPos[7] = 0.5*(pos[6]+pos[8]);
	pushedPos[8] = pos[8];	
	if( diamondBucklingEnabled )
	{
		pushedPos[0] += diamondFactor * compressLeftDiff * _cpNormals[patch.cpIndices[0]];
		pushedPos[6] += diamondFactor * compressLeftDiff * _cpNormals[patch.cpIndices[6]];
		pushedPos[4] += diamondFactor * 0.5*(compressLeftDiff+compressRightDiff) * _cpNormals[patch.cpIndices[4]];
		pushedPos[2] += diamondFactor * compressRightDiff * _cpNormals[patch.cpIndices[2]];
		pushedPos[8] += diamondFactor * compressRightDiff * _cpNormals[patch.cpIndices[8]];
	}

	// Sinusoidal folds
	//pushedPos[7] += (compressTopDiff) * _cpNormals[patch.cpIndices[7]];
	//pushedPos[1] += (compressBottomDiff) * _cpNormals[patch.cpIndices[1]];
	//pushedPos[4] += (0.5*(compressTopDiff+compressBottomDiff)) * _cpNormals[patch.cpIndices[4]];
	
			
	//	Compute max threshold due to collision

	//double colDistCompress = 1.0e+30f;	// infinity
	double colDistMax = -1;
	/*
	double diamDistMax = 0;
	bool colFound;
	if( g_EnableThreshold )
	{
		float distMax0;
		D3DXVECTOR3 colOrig, colDir;
		float dist;
		int j, k;
		//for( j=0; j<4; ++j )
		//{
		//	k = 2*j+1;
		for( k=0; k<9; ++k )
		{
			colOrig = pos[k];
			colDir = -norm[k];
			//D3DXVECTOR3 vDistMax = pushedPos[k] - colOrig;
			distMax0 = 0.1f*BBoxRadius; //D3DXVec3Length(&vDistMax);
			//if( distMax0>distMax )
			//	distMax = distMax0;

			colFound = ( CollideRayDist(&dist, NULL, colOrig, colDir, distMax0, true)>0 );
			if( !colFound )
				dist = 0;

			#ifdef BUCKROSE
				//dist *= 0.5f;
				dist *= 0.9f;
			#else
				dist *= 1.1f;
				//dist *= 0.9f;
			#endif

			if( (k&1)==1 && dist<colDistCompress && colFound )
				colDistCompress = dist;
			if( dist>colDistMax && colFound )
				colDistMax = dist;

			if( (k&1)==1 )
			{
				D3DXVECTOR3 v = pushedPos[k] - colOrig;
				float vd = D3DXVec3Length(&v);
				if( vd>diamDistMax )
					diamDistMax = vd;
			}
		}
	}
	*/
	if( colDistMax<0 )	// used by onTwist
		colDistMax = defaultCollisionMaxDistFactor*_bboxRadius;	
	
	//	2. Interp between the curved surface and the diamonds
	
	//	4 curved patchs are displayed for each patch
	//		+---+---+
	//		| 2 | 3 |
	//		+---X---+
	//		| 0 | 1 |
	//		+---+---+
	int x, y, borderX, borderY, parseBorderX, parseBorderY;
	double tx, ty;
	Vec3d onCurve, onCurveNorm, onDiam;
	Vec3d onTwistDisp;
	Vec3d dispCompress, disp;
	//double dispCompressLen, dispLen;

	CurvedPatch curvedPatch;
	curvedPatch.inverseNormal = false;
	int part, vtxIndex;
	for( part=0; part<4; ++part )
	{
		const Vec3d& n0 = _cpNormals[patch.cpIndices[PARTS[part][0]]];
		const Vec3d& n1 = _cpNormals[patch.cpIndices[PARTS[part][1]]];
		const Vec3d& n2 = _cpNormals[patch.cpIndices[PARTS[part][2]]];
		const Vec3d& n3 = _cpNormals[patch.cpIndices[PARTS[part][3]]];
		Vec3d p0 = controlPoints[patch.cpIndices[PARTS[part][0]]];
		Vec3d p1 = controlPoints[patch.cpIndices[PARTS[part][1]]];
		Vec3d p2 = controlPoints[patch.cpIndices[PARTS[part][2]]];
		Vec3d p3 = controlPoints[patch.cpIndices[PARTS[part][3]]];

		if( diamondBucklingEnabled )
		{
			if( part==0 || part==2 )
			{
				p2 += diamondFactor * compressLeftDiff * n2;
				p3 += diamondFactor * compressLeftDiff * n3;
				p0 += diamondFactor * 0.5*(compressLeftDiff+compressRightDiff) * n0;
				p1 += diamondFactor * 0.5*(compressLeftDiff+compressRightDiff) * n1;
			}
			else // part 1 || 3
			{
				p0 += diamondFactor * compressRightDiff * n0;
				p1 += diamondFactor * compressRightDiff * n1;
				p2 += diamondFactor * 0.5*(compressLeftDiff+compressRightDiff) * n2;
				p3 += diamondFactor * 0.5*(compressLeftDiff+compressRightDiff) * n3;			
			}
		}
				
		//	Vertical sinusoidal folds (parallel to the pseudo-cylinder axis)
		if( verticalBucklingEnabled )
		{
			double topSin = (1.0-twistVersusVerticalFolds)*verticalFoldsFactor*compressTopDiff;
			double bottomSin = (1.0-twistVersusVerticalFolds)*verticalFoldsFactor*compressBottomDiff;
			//double topSin = 0.5*sqrt((2.0*patch.length0[1]+compressTopDiff)*compressTopDiff);
			//double bottomSin = 0.5*sqrt((2.0*patch.length0[0]+compressBottomDiff)*compressBottomDiff);		
			double middleSin = 0.5*(topSin+bottomSin);
			if( part==0 )
			{
				p1 += bottomSin * _cpNormals[patch.cpIndices[1]];
				p0 += middleSin * _cpNormals[patch.cpIndices[4]];
			}
			else if( part==1 )
			{
				p3 += bottomSin * _cpNormals[patch.cpIndices[1]];
				p2 += middleSin * _cpNormals[patch.cpIndices[4]];		
			}
			else if( part==2 )
			{
				p0 += topSin * _cpNormals[patch.cpIndices[7]];
				p1 += middleSin * _cpNormals[patch.cpIndices[4]];		
			}
			else if( part==3 )
			{
				p2 += topSin * _cpNormals[patch.cpIndices[7]];
				p3 += middleSin * _cpNormals[patch.cpIndices[4]];		
			}
		}
				
		curvedPatch.A[0] = p0;
		CurvedPatch::computeCP(&curvedPatch.A[1], p0, p1, n0, n1, patchTangentFactor);
		CurvedPatch::computeCP(&curvedPatch.A[4], p0, p2, n0, n2, patchTangentFactor);
		curvedPatch.A[5] = curvedPatch.A[1] + curvedPatch.A[4] - p0;

		curvedPatch.A[3] = p1;
		CurvedPatch::computeCP(&curvedPatch.A[2], p1, p0, n1, n0, patchTangentFactor);
		CurvedPatch::computeCP(&curvedPatch.A[7], p1, p3, n1, n3, patchTangentFactor); 
		curvedPatch.A[6] = curvedPatch.A[2] + curvedPatch.A[7] - p1;

		curvedPatch.A[12] = p2;
		CurvedPatch::computeCP(&curvedPatch.A[8], p2, p0, n2, n0, patchTangentFactor);
		CurvedPatch::computeCP(&curvedPatch.A[13], p2, p3, n2, n3, patchTangentFactor); 
		curvedPatch.A[9] = curvedPatch.A[8] + curvedPatch.A[13] - p2;

		curvedPatch.A[15] = p3;
		CurvedPatch::computeCP(&curvedPatch.A[11], p3, p1, n3, n1, patchTangentFactor);
		CurvedPatch::computeCP(&curvedPatch.A[14], p3, p2, n3, n2, patchTangentFactor); 
		curvedPatch.A[10] = curvedPatch.A[11] + curvedPatch.A[14] - p3;

		borderY = 1; //( i==0 || i==2 ) ? 1 : 0;
		for( y=0; y<patchSubdiv + borderY; ++y )
		{
			parseBorderY = ( y==patchSubdiv ) ? 1 : 0;
			ty = double(y)/patchSubdiv;
			borderX = 1; //( i==0 || i==1 ) ? 1 : 0;
			for( x=0; x<patchSubdiv + borderX; ++x )
			{
				parseBorderX = ( x==patchSubdiv ) ? 1 : 0;
				tx = double(x)/patchSubdiv;

				//	point onto the curved surface
				curvedPatch.eval(&onCurve, tx, ty);
				curvedPatch.evalNormal(&onCurveNorm, tx, ty);
					
				//	point onto the diamond triangles
				double crossFold = 0;
				if( part==3 )
				{
					if( y>x )
						crossFold = CurvedPatch::baryInterp(&onDiam, pushedPos[7], pushedPos[4], pushedPos[8], tx, 1.0f-ty);
					else
						crossFold = CurvedPatch::baryInterp(&onDiam, pushedPos[5], pushedPos[8], pushedPos[4], 1.0f-tx, ty);
				}
				else if( part==2 )
				{
					if( patchSubdiv-y>x )
						crossFold = CurvedPatch::baryInterp(&onDiam, pushedPos[7], pushedPos[4], pushedPos[6], tx, ty);
					else
						crossFold = CurvedPatch::baryInterp(&onDiam, pushedPos[3], pushedPos[6], pushedPos[4], 1.0f-tx, 1.0f-ty);
				}
				else if( part==1 )
				{
					if( patchSubdiv-y>x )
						crossFold = CurvedPatch::baryInterp(&onDiam, pushedPos[5], pushedPos[2], pushedPos[4], tx, ty);
					else
						crossFold = CurvedPatch::baryInterp(&onDiam, pushedPos[1], pushedPos[4], pushedPos[2], 1.0f-tx, 1.0f-ty);
				}
				else // if( part==0 )
				{
					if( y>x )
						crossFold = CurvedPatch::baryInterp(&onDiam, pushedPos[3], pushedPos[0], pushedPos[4], tx, 1.0f-ty);
					else
						crossFold = CurvedPatch::baryInterp(&onDiam, pushedPos[1], pushedPos[4], pushedPos[0], 1.0f-tx, ty);
				}
				
				//	point on the twisted surface
				onTwistDisp = Vec3d(0,0,0);
				if( twistBucklingEnabled )
				{
					bool onTwistFold = 		( twist<0 && (  ( (part==3 || part==0) && x==y )
															|| (part==1 || part==2) && ( (x==0 && y==patchSubdiv) || (x==patchSubdiv && y==0) ) ) ) 
										||  ( twist>0 && ( 	( (part==1 || part==2) && x==patchSubdiv-y ) 
															|| (part==0 || part==3) && ( (x==0 && y==0) || (x==patchSubdiv && y==patchSubdiv) ) ) ); 
					   
					if( !onTwistFold )
						onTwistDisp = (- 0.5 * twistFactor*colDistMax) * onCurveNorm;
					else
						onTwistDisp = (+ 0.5 * twistFactor*colDistMax) * onCurveNorm;
				}
								
				if( diamondBucklingEnabled )
				{
					if( part==0 || part==2 )
						compressVert = ty*compressLeft + (1.0f-ty)*0.5f*(compressLeft+compressRight);
					else // if( part==1 || part==3 )
						compressVert = ty*0.5f*(compressLeft+compressRight) + (1.0f-ty)*compressRight;
				}
				else
					compressVert = 0;
				compressVert = pow(compressVert, 0.15);
					
				//	blend 
				//dispCompress = compressVert*(onDiam - onCurve);
				dispCompress = (crossFold) * compressVert*(onDiam - onCurve);

				//dispCompress = onHorizSinDisp;
				//dispCompressLen = D3DXVec3Length(&dispCompress);
				/*
				float dispHorizSinLen = D3DXVec3Length(&onHorizSinDisp);
				if( diamDistMax<=dispHorizSinLen )
				{
					dispCompress += ((dispHorizSinLen-diamDistMax)/dispHorizSinLen)*onHorizSinDisp;
					dispCompressLen = D3DXVec3Length(&dispCompress);
				}
				*/
				/*
				if( dispCompressLen>colDistCompress )
					dispCompress *= colDistCompress/dispCompressLen;
				*/

				disp = dispCompress + fabs(twist)*onTwistDisp;
				/*
				if( collisionEnabled )
				{
					dispLen = disp.norm();
					if( dispLen>colDistMax )
						disp *= colDistMax/dispLen;
				}
				*/
							
				vtxIndex = patch._vtxIndices[part][x+y*(patchSubdiv+1)];
				_refinedVertices[vtxIndex] = onCurve + disp;
				if( collisionEnabled && _moveOutsideCB!=NULL && !doCollisionsAtTheEnd )
					_refinedVertices[vtxIndex] = _moveOutsideCB(_refinedVertices[vtxIndex], _moveOutsideClientData);
				//_refinedNormals[vtxIndex] = Vec3d(0,0,0);  //normals should be recomputed by calling updateRefinedNormals()
			}
		}
	}
}


void ControlMesh::updateRefinedNormals()
{
	size_t i;
	for( i=0; i<_refinedNormals.size(); ++i )
		_refinedNormals[i] = Vec3d(0,0,0);

	int part, x, y;		
	int a, b, c, d;
	Vec3d normal;
	for( i=0; i<patchs.size(); ++i )
	{
		const Patch& patch = patchs[i];
		for( part=0; part<4; ++part )
			for( y=0; y<patchSubdiv; ++y )
				for( x=0; x<patchSubdiv; ++x )
				{
					a = patch._vtxIndices[part][(x  )+(y  )*(patchSubdiv+1)];
					b = patch._vtxIndices[part][(x+1)+(y  )*(patchSubdiv+1)];
					c = patch._vtxIndices[part][(x  )+(y+1)*(patchSubdiv+1)];
					d = patch._vtxIndices[part][(x+1)+(y+1)*(patchSubdiv+1)];
					normal = (_refinedVertices[a]-_refinedVertices[d]) ^ (_refinedVertices[c]-_refinedVertices[b]);
					normal.normalizeSafe();
					_refinedNormals[a] += normal;
					_refinedNormals[b] += normal;
					_refinedNormals[c] += normal;
					_refinedNormals[d] += normal;
				}
	}
				
	for( i=0; i<_refinedNormals.size(); ++i )
		_refinedNormals[i].normalizeSafe();
}


void ControlMesh::moveRefinedVerticesOutside()
{
	if( _moveOutsideCB==NULL )
		return;

	for( size_t i=0; i<_refinedVertices.size(); ++i )
		_refinedVertices[i] = _moveOutsideCB(_refinedVertices[i], _moveOutsideClientData);
}


void ControlMesh::drawRefinedPatch(int patchIndex) const
{
	const Patch& patch = patchs[patchIndex];
 
	glPushAttrib(GL_ALL_ATTRIB_BITS);	
	glEnable(GL_LIGHTING);
	
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	float ambient[] = {0.05f, 0.05f, 0.05f, 1.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
  	//float specular[] = {0.1f, 0.1f, 0.1f, 1.0};
	float specular[] = {0, 0, 0, 0};
  	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);
	float diffuse[] = {defaultBackColor[0], defaultBackColor[1], defaultBackColor[2], 1.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glColor3fv(diffuse);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
	float lmAmbient[] = {0.05,0.05,0.05,1};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmAmbient);
	glEnable(GL_LIGHT0);
	float lightAmbient[] = {0,0,0,1};
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	float lightDiffuse[] = {1.05f,1.05f,1.05f,1};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	float lightSpecular[] = {0.5,0.5,0.5,1};
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	float lightDir[] = {0.4,0.2,0.6,0};
	glLightfv(GL_LIGHT0, GL_POSITION, lightDir);
	glPopMatrix();
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	
	if( textureEnabled && patch.textureID>=0 )
	{
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, patch.textureID);
		glMatrixMode(GL_TEXTURE);
		glLoadIdentity();
		glScaled(textureScale[0], textureScale[1], 1);
		glMatrixMode(GL_MODELVIEW);		
	}
	else
		glDisable(GL_TEXTURE_2D);
		
	if( colorPatchsEnabled )
		glColor3dv(COLORS[patchIndex%(sizeof(COLORS)/sizeof(COLORS[0]))]);
	else
		glColor3dv(defaultFrontColor.address());

	glBegin(GL_TRIANGLES);
	
	//	4 curved patchs are displayed for each patch
	int x, y;
	int part, a, b, c, d;
	double tx0, tx1, ty0, ty1;
	Vec2d ta, tb, tc, td;
	Vec2d tsa, tsb, tsc, tsd;	
	for( part=0; part<4; ++part )
	{
		ta = patch.texCoords[PARTS[part][0]];
		tb = patch.texCoords[PARTS[part][1]];
		tc = patch.texCoords[PARTS[part][2]];
		td = patch.texCoords[PARTS[part][3]];
		
		for( y=0; y<patchSubdiv; ++y )
		{
			ty0 = double(y)/patchSubdiv;
			ty1 = double(y+1)/patchSubdiv;			
			for( x=0; x<patchSubdiv; ++x )
			{
				tx0 = double(x)/patchSubdiv;
				tx1 = double(x+1)/patchSubdiv;
							
				tsa = (1.0-ty0)*( (1.0-tx0)*ta+tx0*tb ) + ty0*( (1.0-tx0)*tc+tx0*td );
				tsb = (1.0-ty0)*( (1.0-tx1)*ta+tx1*tb ) + ty0*( (1.0-tx1)*tc+tx1*td );
				tsc = (1.0-ty1)*( (1.0-tx0)*ta+tx0*tb ) + ty1*( (1.0-tx0)*tc+tx0*td );
				tsd = (1.0-ty1)*( (1.0-tx1)*ta+tx1*tb ) + ty1*( (1.0-tx1)*tc+tx1*td );
				
				a = patch._vtxIndices[part][(x  )+(y  )*(patchSubdiv+1)];
				b = patch._vtxIndices[part][(x+1)+(y  )*(patchSubdiv+1)];
				c = patch._vtxIndices[part][(x  )+(y+1)*(patchSubdiv+1)];
				d = patch._vtxIndices[part][(x+1)+(y+1)*(patchSubdiv+1)];
				
				if( part==0 || part==3 )
				{
					glTexCoord2dv( tsa.address() );
					glNormal3dv( _refinedNormals[a].address() );
					glVertex3dv( _refinedVertices[a].address() );
					glTexCoord2dv( tsd.address() );					
					glNormal3dv( _refinedNormals[d].address() );
					glVertex3dv( _refinedVertices[d].address() );
					glTexCoord2dv( tsb.address() );					
					glNormal3dv( _refinedNormals[b].address() );
					glVertex3dv( _refinedVertices[b].address() );
	
					glTexCoord2dv( tsa.address() );
					glNormal3dv( _refinedNormals[a].address() );
					glVertex3dv( _refinedVertices[a].address() );
					glTexCoord2dv( tsc.address() );					
					glNormal3dv( _refinedNormals[c].address() );
					glVertex3dv( _refinedVertices[c].address() );
					glTexCoord2dv( tsd.address() );
					glNormal3dv( _refinedNormals[d].address() );
					glVertex3dv( _refinedVertices[d].address() );												
				}
				else
				{
					glTexCoord2dv( tsa.address() );
					glNormal3dv( _refinedNormals[a].address() );
					glVertex3dv( _refinedVertices[a].address() );
					glTexCoord2dv( tsc.address() );
					glNormal3dv( _refinedNormals[c].address() );
					glVertex3dv( _refinedVertices[c].address() );
					glTexCoord2dv( tsb.address() );
					glNormal3dv( _refinedNormals[b].address() );
					glVertex3dv( _refinedVertices[b].address() );
	
					glTexCoord2dv( tsd.address() );
					glNormal3dv( _refinedNormals[d].address() );
					glVertex3dv( _refinedVertices[d].address() );
					glTexCoord2dv( tsb.address() );
					glNormal3dv( _refinedNormals[b].address() );
					glVertex3dv( _refinedVertices[b].address() );
					glTexCoord2dv( tsc.address() );
					glNormal3dv( _refinedNormals[c].address() );
					glVertex3dv( _refinedVertices[c].address() );																
				}
			}
		}
	}
	
	glEnd();

	glDisable(GL_COLOR_MATERIAL);		
	glDisable(GL_TEXTURE_2D);
	glPopAttrib();	
}


void ControlMesh::computeBValues(double *outBValues, int numOfB, const double *radii, int numOfRadii)
{
	assert( numOfRadii>0 );
	assert( radii[0]>0 );
	assert( numOfB>0 );
	const int NUM_SAMPLES = 500;
	const double B_FACTOR = 5.0;
	double b0, sumOfB, lambda, deltaIdx, r;
	double closestSumOfB = -1.0;
	int closestSample = int(double(NUM_SAMPLES)/B_FACTOR);
	int rIdx0, rIdx1;
	for( int sample=1; sample<NUM_SAMPLES; ++sample )
	{
		b0 = (double(sample)/NUM_SAMPLES) * (B_FACTOR/double(numOfB));
		lambda = b0 / sqrt(radii[0]);
		sumOfB = 0;
		for( int i=0; i<numOfB; ++i )
		{
			rIdx0 = (i*numOfRadii)/numOfB;
			rIdx1 = (i*numOfRadii)/numOfB + 1;
			deltaIdx = double(i*numOfRadii)/numOfB - rIdx0;
			if( rIdx1>=numOfB )
				r = radii[rIdx0];
			else
				r = (1.0-deltaIdx)*radii[rIdx0] + deltaIdx*radii[rIdx1];
			sumOfB += lambda * sqrt(r);
		}
		if( closestSumOfB<=0 || fabs(sumOfB-1.0)<fabs(closestSumOfB-1.0) )
		{
			closestSumOfB = sumOfB;
			closestSample = sample;
		}
	}
	
	b0 = (double(closestSample)/NUM_SAMPLES) * (B_FACTOR/double(numOfB));
	lambda = b0 / sqrt(radii[0]);
	sumOfB = 0;
	for( int i=0; i<numOfB; ++i )
	{
		rIdx0 = (i*numOfRadii)/numOfB;
		rIdx1 = (i*numOfRadii)/numOfB + 1;
		deltaIdx = double(i*numOfRadii)/numOfB - rIdx0;
		if( rIdx1>=numOfB )
			r = radii[rIdx0];
		else
			r = (1.0-deltaIdx)*radii[rIdx0] + deltaIdx*radii[rIdx1];
		outBValues[i] = ( lambda * sqrt(r) ) / closestSumOfB;
		sumOfB += outBValues[i];
	}
	cout << "sumOfB=" << closestSumOfB << "  sample=" << closestSample << endl;
}


// for testing purpose
struct MoveOutsideParam
{
	Vec3d bboxMin, bboxMax;
	double radiusRatio;
};
Vec3d moveOutsideCB(const Vec3d& point, void *param)
{
	if( param==NULL )
		return point;
	MoveOutsideParam *p = static_cast<MoveOutsideParam *>(param);
	Vec3d center = (p->bboxMin+p->bboxMax)/2;
	double rx = (p->bboxMax[0]-p->bboxMin[0])/2;
	double rz = (p->bboxMax[2]-p->bboxMin[2])/2;
	double r = p->radiusRatio * ( (rx>rz) ? rx : rz );
	Vec3d v = point - Vec3d(center[0], point[1], center[2]);
	double len = v.norm();
	double newLen = (len<r) ? r : len;
	return (newLen/len)*v + Vec3d(center[0], point[1], center[2]);
}

// To be called in pattern.cpp: int Pattern::prepareModel() ->  if(!_disp_pattern) { /*...*/ HERE }
//		In Pattern::init() :
//			setBackgroundColor(QColor(255,255,255));	// for paper figures (Phil)
//			setBackgroundColor(QColor(0x80,0xa7,0xd0));	// for video  0xff80a7d0 (Phil)
//		In Pattern::draw() : 
//			resize(640, 480);	// for video (Phil)
void bucklingTest(const Vec3d& bbMin, const Vec3d& bbMax, QGLViewer *viewer)
{
	static ControlMesh s_CtrlMesh;
	static bool s_Initialized = false;

	static int s_ViewerInitialized = 0;
	if( s_ViewerInitialized<3 )
	{
		viewer->resize(640, 480);	
		//viewer->setBackgroundColor(QColor(255,255,255));		// for paper figures (Phil)
		viewer->setBackgroundColor(QColor(0x80,0xa7,0xd0));		// for video  0xff80a7d0 (Phil)
		viewer->setSnapshotFormat("PNG");		
		viewer->setSnapshotFileName("MOVIES/buckling2.png");
		++s_ViewerInitialized;
	}
	
	if( !s_Initialized )
	{
		static int s_Counter = 0;
	
		// Create a texture by hand
		GLuint textureID = 0;
		/*
		GLuint texData[64] = { 	0xffff00ff, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0,
								0xffff00ff, 0xff00ffff, 0xff00ffff, 0xff00ffff, 0xff00ffff, 0xff00ffff, 0xff00ffff, 0xffa0a0a0,
								0xffff00ff, 0xff00ffff, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0,
								0xffff00ff, 0xff00ffff, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0,
								0xffff00ff, 0xff00ffff, 0xff00ffff, 0xff00ffff, 0xff00ffff, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0,
								0xffff00ff, 0xff00ffff, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0,
								0xffff00ff, 0xff00ffff, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0,
								0xff0000ff, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0, 0xffa0a0a0 };
		glGenTextures(1, &textureID);
		glBindTexture(GL_TEXTURE_2D, textureID);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, 8, 8, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
		gluBuild2DMipmaps(GL_TEXTURE_2D, 4, 8, 8, GL_RGBA, GL_UNSIGNED_BYTE, texData);
    	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		s_CtrlMesh.textureScale = Vec2d(10,10);
		*/

		// Create a cylindrical control mesh		
		s_CtrlMesh.buildCylinder(8, 5, bbMin, bbMax, 0.5, 0.5, textureID);
		//s_CtrlMesh.jitterAlongAxis(Vec3d((bbMin[0]+bbMax[0])/2, bbMax[1], (bbMin[2]+bbMax[2])/2), Vec3d((bbMin[0]+bbMax[0])/2, bbMin[1], (bbMin[2]+bbMax[2])/2), 0.1, 0.1, 0.3, 1000);
		s_CtrlMesh.resetPatchsLength0();

		//double comp = 1 - 0.8*  (1-cos(M_PI*double(s_Counter)/50))/2;
		//s_CtrlMesh.compressVert(comp);
		//s_CtrlMesh.compressVert(0.5);
		//s_CtrlMesh.twistVert(0, 0.5*M_PI);
		s_CtrlMesh.compressHoriz(0.8, 0.8, 1, 1);
		//s_CtrlMesh.pushBackAlongAxis(Vec3d((bbMin[0]+bbMax[0])/2, bbMax[1], (bbMin[2]+bbMax[2])/2), Vec3d((bbMin[0]+bbMax[0])/2, bbMin[1], (bbMin[2]+bbMax[2])/2), 0.7*(bbMax[1]-bbMin[1]), 0.2*(bbMax[1]-bbMin[1]), 0.2*(bbMax[1]-bbMin[1]), 1000);
		//s_CtrlMesh.twistAlongAxis(Vec3d(0.7*bbMin[0]+0.3*bbMax[0], bbMin[1], (bbMin[2]+bbMax[2])/2), Vec3d(0.7*bbMin[0]+0.3*bbMax[0], bbMax[1], (bbMin[2]+bbMax[2])/2), 0, 0.4*M_PI, 1000);
		//s_CtrlMesh.scaleAlongAxis(Vec3d((bbMin[0]+bbMax[0])/2, bbMax[1], (bbMin[2]+bbMax[2])/2), Vec3d((bbMin[0]+bbMax[0])/2, bbMin[1], (bbMin[2]+bbMax[2])/2), 1.0, 0.5, 1000);
		//double angle = 0.5*M_PI*  (cos(M_PI*double(s_Counter)/50)-1)/2;
		//s_CtrlMesh.bendBone(Vec3d(0.47*bbMin[0]+0.53*bbMax[0], 0.5*bbMin[1]+0.5*bbMax[1], 0.5*bbMin[2]+0.5*bbMax[2]), Vec3d(0.47*bbMin[0]+0.53*bbMax[0], 1.1*bbMin[1]-0.1*bbMax[1], 0.5*bbMin[2]+0.5*bbMax[2]), Vec3d(0,0,1), angle, 0.3*(bbMax[1]-bbMin[1]), 0.35, 1000);
		
		static MoveOutsideParam s_MoveOutsideParam = { bbMin, bbMax, 0.39 };
		//s_CtrlMesh.update(moveOutsideCB, &s_MoveOutsideParam);
		s_CtrlMesh.update();
		
//		s_Initialized = true;
		++s_Counter;		
	}
	 
	// Draw
	//s_CtrlMesh.drawControlMesh();	
	//s_CtrlMesh.drawCurvedSurface();
	s_CtrlMesh.drawBucklingSurface();
	
	if( 0 ) // s_ViewerInitialized>=3 )
	{
		viewer->saveSnapshot(true, true);
		cout<<"SNAP "<<viewer->snapshotFileName().ascii()<<" "<<viewer->snapshotCounter()<<endl;
	}
}
