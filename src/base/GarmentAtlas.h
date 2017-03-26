#ifndef  GARMENT_ATLAS_H
#define GARMENT_ATLAS_H

/* 

Holds a reference to the previous model representation (LCModel) and also
a description of the pattern (atlas) which allows us to find the external
edges.
*/

#include <string>
#include <map>
#include <set>
#include <vector>
//#include "Face.h"
//#include "Edge.h"
//#include "Vertex.h"
#include <utility>
#include "Parabola.h"
#include "buckling.h"
#include "vectypes.h"

class LCModel;
class EdgeSection;
class GarmentSection;
class Garment;
class pointDistance;
namespace polymesh {
class Face;
class Edge;
class Vertex;
}



class GarmentAtlas
{
  typedef std::vector<polymesh::Vertex*> VertexPList;
  typedef std::vector<polymesh::Face*> FacePList;
  
  public:
    
    // MEMBERS
    bool _edgeSectionsCategorised;
    bool _dartPointsDetermined;
    bool _atlas_populated;
    bool _controlMeshAssigned;
    bool _coincidentPointSearchVisited;
    bool _boundaryEdgesLinked;
    bool _seamEdgesDetermined;
    
    ControlMesh _controlMesh;
    std::string name;
    LCModel* _lcmodel;
    GarmentSection* parent;
    std::vector<Parabola*> parabolas; // parabolas to fit bottom [middle] and top of pattern. Middle only when required (Tshirts with sleeves etc)
    Vec3d _axisPointA;
    Vec3d _axisPointB;
    
    std::vector<EdgeSection*> _leftEdgeSections, _rightEdgeSections,_topEdgeSections,_bottomEdgeSections; 
    std::vector<long> _dartPoints; // Store indexes to darts for later comparison
    std::set<long> _dartPointsSet; // Stores indexes for quick lookup
    
    std::vector<polymesh::Edge*> _unsortedBoundaryEdges;
    std::vector<EdgeSection*> _listOfSeamEdgeSections;
    
    // METHODS
    
    GarmentAtlas() :  _edgeSectionsCategorised(false),
                      _dartPointsDetermined(false), 
                      _atlas_populated(false),
                      _controlMeshAssigned(false), 
                      _coincidentPointSearchVisited(false), 
                      _boundaryEdgesLinked(false), 
                      _seamEdgesDetermined(false) 
    {
      _vertexList.reserve(3000);
      _edgeList.reserve(1000);
      _faceList.reserve(2000);
      _unsortedBoundaryEdges.reserve(500);
      parabolas.reserve(3);
      leftEdgesReq = 1;
      rightEdgesReq = 1;
      topEdgesReq = 1;
      bottomEdgesReq =1;
      _lcmodel = 0;
      parent = 0;
      _averageMeshVectorLength = 0.0;

    };
    
    ~GarmentAtlas();

    

    bool addVertex(polymesh::Vertex*);
    bool addEdge(polymesh::Edge*);
    bool addFace(polymesh::Face*);
    
    // Lookup an Edge by it's two endpoints. Order doesn't matter.
    polymesh::Edge* lookup(long,long);
    int leftEdgesReq, rightEdgesReq, topEdgesReq, bottomEdgesReq;
    
    // for rendering a smoother pattern controlmesh version
    std::vector<Vec3d> interpolatedPatternPoints; //JDW
    
    VertexPList vertices();
    FacePList faces();
    
    std::vector<polymesh::Edge*> edges() { return _edgeList; };
    
    // once model has been assigned generate Vertex,Face and Edge descriptions
    void populateAtlasFromModel();
    
    void centreModelOnOrigin(); // needed for rotations to work
    // Returns the index of the closest texture co-ordinate 
    long findClosestTCTo(double x, double y);
    long findClosestTCToIgnoringDarts(double x, double y);
    // returns index of closest plus populates npArray with closest three
    long findClosest3TCTo(double x, double y, pointDistance *npArray);
    void determineCoincidentPoints(Garment*);
    void determineEdgesUsingSeams(); // uses coincident point counts to categorise boundary edges.
    void linkBoundaryEdges();
    void reduceEdgeSectionsToRequiredNumber(int);
    void fitParabola();
    void categoriseSeamEdgeSections();
    void updateEdgeSectionInformation(); // call after co-ordinate changes.
    
    void findBarycentricCoefficients(Vec3d, Vec3d, Vec3d, Vec3d, double*, double*);
    
    void populateControlMesh(int rows, int cols);

    double lengthOfInterpolatedParabola(Parabola *ppa, Parabola *ppb, double alpha, int samples) const;
    // Call after atlas is populated. Compare all boundary points to each other and mark darts
    // (where 3d mesh points are co-incident)

    
    
  private:
    
    VertexPList _vertexList;
    std::vector<polymesh::Edge*> _edgeList;
    FacePList _faceList;
    
    
    std::map<std::pair<long, long>, polymesh::Edge*> _edgeLookupMap; 
    
    double _averageMeshVectorLength;
    //void categoriseEdgeSections(); // obsolete but here in case required again
    // Rescale the Atlas to it's real dimensions (previously in 0-1 range for both co-ords)
    void scaleAtlasFromMesh();
    
    void determineDartPoints();


    
    
};

#endif // GARMENT_ATLAS_H
