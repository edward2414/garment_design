//
// C++ Implementation: Garment
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
#include "GarmentAtlas.h"
#include "GarmentSection.h"
#include "lcmodel.h"

Garment::~Garment(void)
{
  std::vector<GarmentSection*>::iterator gsi, gsi_end;
  gsi= sectionList.begin();
  gsi_end= sectionList.end();
  for(;gsi!=gsi_end;gsi++)
  {
    delete *gsi;
  }
}

GarmentAtlas* Garment::getAtlas(std::string name)
{
  GarmentAtlas* pFoundGarmentAtlas=0;
      
  std::vector<GarmentSection*>::const_iterator gsi, gsi_end;
  gsi= sectionList.begin();
  gsi_end= sectionList.end();
  for(;gsi!=gsi_end;gsi++)
  {
    GarmentSection* pGarmentSection = *gsi;
    // returns null pointer if not found.
    pFoundGarmentAtlas = pGarmentSection->getAtlas( name );
    if(!pFoundGarmentAtlas)
      continue;
    return pFoundGarmentAtlas;
  }
  return pFoundGarmentAtlas;
}

void Garment::determineCoincidentPoints() // FIXME rename
{
  // Can only be done once all atlases patterns have been correctly populated
  
  std::vector<GarmentSection*>::const_iterator gsi;
  for(gsi = sectionList.begin(); gsi != sectionList.end(); gsi++)
  {
    
    std::vector<GarmentAtlas*>::const_iterator gai;
    for(gai = (*gsi)->atlasList.begin(); gai != (*gsi)->atlasList.end(); gai++)
    {
      if((*gai)->_unsortedBoundaryEdges.size() <= 0) {
        std::cout << "Can't determine coincident points until all atlases complete" << std::endl;
        return;
      }
    }
  }
  
  
  
  
  for(gsi = sectionList.begin(); gsi != sectionList.end(); gsi++)
  {
    std::vector<GarmentAtlas*>::const_iterator gai;
    for(gai = (*gsi)->atlasList.begin(); gai != (*gsi)->atlasList.end(); gai++)
    {
      (*gai)->determineCoincidentPoints(this);
    }
  }
  
  for(gsi = sectionList.begin(); gsi != sectionList.end(); gsi++)
  {
    
    std::vector<GarmentAtlas*>::const_iterator gai;
    for(gai = (*gsi)->atlasList.begin(); gai != (*gsi)->atlasList.end(); gai++)
    {
      (*gai)->determineEdgesUsingSeams();
      
    }
  }
  
}


void Garment::populateControlMesh()
{
  std::cout << "Populating control mesh for whole garment" << std::endl;
  // Can only be done once all atlases patterns have been correctly populated (edges categorised etc).
  // FIXME add assertion
  
  
  controlMesh.reset();
  
  std::vector<GarmentSection*>::const_iterator gsi, gsi_end;
  
gsi=sectionList.begin();
gsi_end=sectionList.end();
for(;gsi!=gsi_end;gsi++)
{
  GarmentSection* pGS=*gsi;
  if(!pGS->controlMeshAssigned)
  {
    std::cout << "Can't populate garment controlmesh until all sections assigned" << std::endl;
    return;
  }
}
  
  
  gsi=sectionList.begin();
  gsi_end=sectionList.end();
  for(;gsi!=gsi_end;gsi++)
  {
    int totalAtlases = (*gsi)->atlasList.size();
    for(int a=0; a<totalAtlases;a++)
    { 
      const ControlMesh& acm = (*gsi)->atlasList[a]->_controlMesh;
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
        if(name == "skirtdarts")
        
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
  }

  controlMesh.resetPatchsLength0();
  controlMesh.update();
  _origControlMesh = controlMesh;
//  isValid = true;
 

} // End populateControlMesh()

void Garment::addSection(std::string name, GarmentSection* pSection)
{

  if(!pSection) return;
  sectionMap[name] = pSection;
  sectionList.push_back(pSection);
  pSection->parent = this;
}

void Garment::resetControlMesh()
{
  controlMesh = _origControlMesh;
}

void Garment::printHierarchy()
{
  std::cout << "Garment: " << name << std::endl;
  std::vector<GarmentSection*>::const_iterator gsi, gsi_end;
  gsi= sectionList.begin();
  gsi_end= sectionList.end();
  for(;gsi!=gsi_end;gsi++)
  {
    GarmentSection* pGarmentSection = *gsi;
    std::cout << " GarmentSection: " << pGarmentSection->name << std::endl;
    std::map<std::string, GarmentAtlas*>::const_iterator gai, gai_end;
    gai= pGarmentSection->atlasMap.begin();
    gai_end= pGarmentSection->atlasMap.end();
    for(;gai!=gai_end;gai++)
    {
      GarmentAtlas* pGarmentAtlas = pGarmentSection->getAtlas ( (*gai).first );
      std::cout << "  GarmentAtlas: " << pGarmentAtlas->name << std::endl;
      std::cout << "   VMeshCount:" << pGarmentAtlas->_lcmodel->points().size() << std::endl;
      std::cout << "   TCCount:" << pGarmentAtlas->_lcmodel->tex_coords().size() << std::endl;
    }
  }
}
