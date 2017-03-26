//
// C++ Interface: buckling
//
// Author: Philippe Decaudin 
//

 
#ifndef BUCKLING_H
#define BUCKLING_H

//#define BUCKLING_TEST		// uncomment to do tests
 
#include <vector>
#include <vectypes.h>
     
struct ControlMesh
{
	// 	Values that change the behaviour of the refined surface processing
	int		patchSubdiv;				// patch subdivision (controls the number of facets used to display a patch)
	double	patchTangentFactor;			// = 2/3 (or less)
	double	twistFactor;				// = 1 (or more)
	double	diamondFactor;				// = 0 (or more)
	double	defaultCollisionMaxDistFactor;	// = 0.1
	double	verticalFoldsFactor;		// = 1 (or less)
	double  twistVersusVerticalFactor;	// = 100
	double  diamondVersusVerticalFactor;// = 100	
	bool	diamondBucklingEnabled;
	bool	twistBucklingEnabled;
	bool	verticalBucklingEnabled;
	bool	smoothFilterEnabled;
	bool	collisionEnabled;
	bool	doCollisionsAtTheEnd;
	bool	colorPatchsEnabled;
	Vec3d	defaultFrontColor;
	Vec3d	defaultBackColor;
	bool	textureEnabled;
	Vec2d	textureScale;
	
	// 	JDW added so we can promote meshes up the garment hierarchy
    int rows;
    int cols;

	//	Reset the control mesh
	void reset();
	
	//	addPatch(..): Adds a new patch to the control mesh
	//	cpArray is an array of 9 control points. 
	//	cpSeamIdArray is an array of 9 integers telling if the corresponding control point is on a seam line or not: if id==0 it is not a seam point, if id>0 it is a seam point and it can be weld only with other seam points having the same id.
	//	If a control points is not a seam and is closer than cpInsideEps from an existing non-seam control point, they will be welded.		
	//	If a control points is a seam and is closer than cpSeamEps from an existing seam control point with same id, they will be welded.
	//	length0Array is an array of 6 double that corresponds to the intial edges lengths of the patch (if NULL they are computed from the control points positions).
	//	The order of 9 control points and the 6 lengths is important (see struct Patch)
 	int addPatch(const Vec3d *cpArray, double cpInsideEps, const int *cpSeamIdArray=NULL, double cpSeamEps=0, const double *length0Array=NULL, const Vec2d *texCoordArray=NULL, int textureID=-1, const int *cpAtlasIndices=NULL, const std::string atlasName=NULL, const Vec3d* interpPatternPoints=NULL);
	
	int getNumPatchs() const;
	int getNumControlPoints() const;
	
	// getPatch(..): outCPArray=array of 9 Vec3d or NULL, outCPIndicesArray=array of 9 int or NULL, outCPIsSeamArray=array of 9 bool or NULL, outLength0Array=array of 6 double or NULL, outTexCoordArray=array of 9 Vec2d or NULL, outTextureID=pointer to an int or NULL.
	void getPatch(int patchIndex, Vec3d *outCPArray, int *outCPIndicesArray, int *outCPSeamIdArray, double *outLength0Array, Vec2d *outTexCoordArray, int *outTextureID, int *outAtlasIndices, std::string *outAtlasName, Vec3d *outInterpPatternPoints) const;
	
	const Vec3d& getControlPoint(int cpIndex) const;
	void setControlPoint(int cpIndex, const Vec3d& cp);
	
	//	getBBox: Computes the control points bounding box. Returns the bounding box radius, and min and max points.
	double getBBox(Vec3d *outBBoxMin=NULL, Vec3d *outBBoxMax=NULL) const;
	
	//	recompute patchs intial lenghts from current positions of control points:
	void resetPatchsLength0();
	   	
	// 	update(..) should be called after populating the control mesh. Then, call one of the drawing methods.
	typedef Vec3d (* MoveOutsideCallback)(const Vec3d& point, void *clientData);
	void update(MoveOutsideCallback moveOutsideCB=NULL, void *clientData=NULL);
	
	//	drawing methods
	void drawControlMesh(bool drawControlPoints=true, bool drawNormals=false, bool drawAdjacencies=false);
	void drawCurvedSurface();
	void drawBucklingSurface();	
	
	// 	compute numOfB values of B (B is the height of a buckling patch) so that the sum of the Bs is 1.
	static void computeBValues(double *outBValues, int numOfB, const double *radii, int numOfRadii);	// 'radii' is an array of an arbitrary number of radii taken at equal distances along the vertical axis from bottom to top.
		
	//	Procedural moves of control points
	
	void compressVert(double scale);
	void twistVert(double angleBottom, double angleTop);
	void compressHoriz(double scaleBottomX, double scaleBottomZ, double scaleTopX, double scaleTopZ);
	
	// 	push back control points along an axis:
	//	Total length is shrinked from length(endAxis-startAxis) to shrinkedLength along the axis. Rigid parts are uncompressed. Points further than 'influenceRadius' are not moved.
	//	 startAxis|-------------------------------------->endAxis
	//	          |<---shrinkedLength---->|- - - - - - - |
	//	          |<--a-->|       |<--b-->|                    startRigidLength=a, endRigidLength=b
	void pushBackAlongAxis(const Vec3d& startAxis, const Vec3d& endAxis, double shrinkedlength, double startRigidLength, double endRigidLength, double influenceRadius);
	
	//	twist control points around an arbitrary axis:
	//	the axis is a segment from point startAxis to point endAxis. Points further than 'influenceRadius' from the axis are not moved.
	void twistAlongAxis(const Vec3d& startAxis, const Vec3d& endAxis, double startAngle, double endAngle, double influenceRadius);

	//	scale control points along an arbitrary axis:
	//	the axis is a segment from point startAxis to point endAxis. Points further than 'influenceRadius' from the axis are not moved.	
	void scaleAlongAxis(const Vec3d& startAxis, const Vec3d& endAxis, double startScale, double endScale, double influenceRadius);
	
	//	bend a section of control points (between pivot and endBone, around rotAxis)
	//	roundPow makes the bend interpolation more or less rounded.
	void bendBone(const Vec3d& pivot, const Vec3d& endBone, const Vec3d& rotAxis, double rotAngle, double blendLength, double roundPow, double influenceRadius);

	//	jitter control points along an arbitrary axis:
	void jitterAlongAxis(const Vec3d& startAxis, const Vec3d& endAxis, double axisJitter, double radiusJitter, double angleJitter, double influenceRadius);

	typedef void (* DistanceToColliderCallback)(double x, double y, double z, double& dist);
	void moveTopControlPoints(double distanceMin, DistanceToColliderCallback distCB);
		
		
	// 	for testing purpose
	void buildCylinder(int n, int m, const Vec3d& bbMin, const Vec3d& bbMax, double rMin, double rMax, int textureID);

	ControlMesh();	
		
protected:
	struct Patch
	{
		int cpIndices[9];	// Indices of the 9 control points of the buckling patch:
							//	6----7----8 
							//	|',  |  ,'| 
							//	|  ',|,'  | 
							//	3----4----5
							//	|  ,'|',  |
							//	|,'  |  ',|
							//	0----1----2
		
		Vec2d texCoords[9];	// texture coordinates linked to control points (same indices order as cpIndices)
		int textureID;		// texture ID that will be used in glBindTexture for this patch
							
		double length0[6];	// initial lengths 
							// (0=bottom, 1=top, 2=left, 3=right, 4=diag topleft/bottomright, 5=diag bottomleft/topright)
		
        int cpAtlasIndices[9]; // holds the index for the point from the corresponding atlas
        std::string atlasName; // in case we are dealing with patches in a Section or Garment we need to get back to the correct atlas to lookup the atlas index							
        Vec3d interpPatternPoints[9]; // for storing interpolated control points in 2d for a smoother rendering of 2d control mesh										
		Patch();
		
	private:
		friend class ControlMesh;
		int _adjacencies[4];			// indices of its 4 adjacent patchs	(0=bottom, 1=top, 2=left, 3=right)
		std::vector<int> _vtxIndices[4];// vertex indices of the refined mesh (refers to _refinedVertices and _refinedNormals)
										// _vtxIndices[n].size()=(patchSubdiv+1)*(patchSubdiv+1)
	};
		 
	// 	Arrays to be filled to populate the control mesh
	std::vector<Vec3d>	controlPoints;
	std::vector<int>	cpSeamIDs;	// id==0 means that it is not a seam point
	std::vector<Patch>	patchs;

private:
	int addControlPoint(const Vec3d& cp, int seamID, double eps, int *currentPatchCPIndices=NULL);
	std::vector<Vec3d>	_cpNormals;	// control points' normals
	void computeNormals();
	void computeAdjacencies();
	Vec3d _bboxMin, _bboxMax;
	double _bboxRadius;
	void computeBBox();
	MoveOutsideCallback _moveOutsideCB;
	void *_moveOutsideClientData;
	void moveControlPointsOutside();
	
	void updateRefinedCurvedPatch(int patchIndex);
	void updateRefinedBucklingPatch(int patchIndex);
	void drawRefinedPatch(int patchIndex) const;	
	void updateRefinedNormals();
	void moveRefinedVerticesOutside();
	void smoothRefinedMesh();
	
	struct VtxAdj	// stores the indices of adjacent vertices of a given refined vertex.
	{
		enum Max		{ MAX = 12 };	// maximum number of adjacent vertices
		int				count[2];
		int				indices[2][MAX];
		inline			VtxAdj() { count[0] = 0; count[1] = 0; }
		void			append(int table, int i);
	};
	
	void initRefinedMesh();	// initialize Patch::_vtxIndices, and _refinedVertices and _refinedNormals
	int addRefinedVertex(const Vec3d& vertex, const Vec3d& normal);	// add vertex to _refinedVertices and _refinedNormals if it has not already been done.
	void computeRefinedAdjacencies();
	std::vector<Vec3d>	_refinedVertices;
	std::vector<Vec3d>	_refinedNormals;		// refined vertices' normals  (_refinedNormals.size()=_refinedVertices.size())
	std::vector<VtxAdj>	_refinedAdjacencies;	// refined vertices' adjacencies (_refinedNormals.size()=_refinedVertices.size())
};


// for test puprose
void bucklingTest(const Vec3d& bbMin, const Vec3d& bbMax, class QGLViewer *qglViewer);


#endif // BUCKLING_H
