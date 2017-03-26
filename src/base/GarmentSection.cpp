//
// C++ Implementation: GarmentSection
//
// Description: 
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "Garment.h"
#include "GarmentSection.h"
#include "GarmentAtlas.h"
#include <vector>
#include <cassert>

using namespace std;

GarmentSection::~GarmentSection(void)
{
  std::vector<GarmentAtlas*>::iterator gai;
  for(gai = atlasList.begin(); gai != atlasList.end(); gai++)
  {
    delete *gai;
  }
}

void GarmentSection::resetControlMesh()
{
  controlMesh = _origControlMesh;
}

void GarmentSection::addAtlas(std::string name, GarmentAtlas* pAtlas)
{
  assert(pAtlas);
  assert(name != "");
  atlasMap[name] = pAtlas;
  atlasList.push_back(pAtlas);
  pAtlas->parent = this;
}

void GarmentSection::populateControlMesh()
{
  std::cout << "Populating control mesh for section: " << name << std::endl;
  // Can only be done once all atlases patterns have been correctly populated (edges categorised etc).
  std::vector<GarmentAtlas*>::const_iterator gai;
  for(gai = atlasList.begin(); gai != atlasList.end(); gai++)
  {
    if(!(*gai)->_controlMeshAssigned) {
      std::cout << "Can't populate Section until all atlases complete" << std::endl;
      return;
    }
  }
  
  // loop over atlases, pushing control points onto aggregate controlmesh while skipping seam points
  // Assumes seams at either vertical edge of each control mesh.
  // Assumes all GarmentAtlases have the same number of rows.
  int iRows = 0;
  int iCols = 0; // populated below
  
  int totalAtlases = atlasList.size();
  
  vector<int> viCols;
  vector<int> viSubCols;
  
  cout << "Initialising counts" << endl;
  // Assume same rows and parabolas for all atlases in a section
  iRows = atlasList.at(0)->_controlMesh.rows;
  for(int i=0;i<totalAtlases;i++)
  {
    viCols.push_back(atlasList[i]->_controlMesh.cols);
    viSubCols.push_back(viCols[i]*2+1);
    iCols += viCols[i];
  }


  // Phil: add atlases patchs to this controlMesh    
  controlMesh.reset();
  // FIXME Note this order is important, forced from file load
  for(int a=0; a<totalAtlases;a++)
  { 
    const ControlMesh& acm = atlasList[a]->_controlMesh;
    for(int i=0; i<acm.getNumPatchs(); ++i)
    {
      Vec3d  patchCP[9];
      Vec3d interpPatternPoints[9];
      int   patchIsSeamCP[9];
      double patchLength0[6];
      Vec2d  patchTexCoord[9];
      int    patchTexID;
      int    patchAtlasIndices[9];
      std::string atlasName;
      acm.getPatch(i, patchCP, NULL, patchIsSeamCP, patchLength0, patchTexCoord,
                   &patchTexID, patchAtlasIndices, &atlasName, interpPatternPoints);
                
                // estimate an espilon that will be used to weld inside control points
      double cpInsideEps =
            ((patchLength0[0]+patchLength0[1]+patchLength0[2]+patchLength0[3])/4)/16;

                // estimate an espilon that will be used to weld seam control points
      double cpSeamEps;
      if(parent->name == "skirtdarts")
        
      {
        cpSeamEps = ((patchLength0[0]+patchLength0[1]+patchLength0[2]+patchLength0[3])/8);
        std::cout << "Using skirt WELD RADIUS for seam welding sections" << std::endl;
      }
      else
      { // tshirt with unique seam welds - can afford greater range
        cpSeamEps = ((patchLength0[0]+patchLength0[1]+patchLength0[2]+patchLength0[3])*2);
        std::cout << "Using tshirt WELD RADIUS for seam welding sections" << std::endl;
      }
                                
      controlMesh.addPatch(patchCP, cpInsideEps, patchIsSeamCP, cpSeamEps,
                           patchLength0, patchTexCoord, patchTexID, patchAtlasIndices, atlasName, interpPatternPoints);
    }
  }
  
    
  
  assert( controlMesh.getNumPatchs()==iCols*iRows );
  
  // FORCE length update as atlas ones are invalid.
  controlMesh.resetPatchsLength0();
  controlMesh.update();
  _origControlMesh = controlMesh;
  isValid = true;
  controlMeshAssigned=true;
 

} // End populateControlMesh()

