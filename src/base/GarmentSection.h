#ifndef  GARMENT_SECTION_H
#define GARMENT_SECTION_H

//
// C++ Interface: GarmentSection
//
// Description: 
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <string>

#include <map>
#include <vector>
#include "buckling.h"

class Garment;
class GarmentAtlas;

class GarmentSection {

  public:
    
    bool isValid;
    bool controlMeshAssigned;
        
    GarmentSection(void) : isValid(false), controlMeshAssigned(false) {};
    ~GarmentSection(void);

    
    std::string name;
    std::map<std::string, GarmentAtlas*> atlasMap;
    std::vector<GarmentAtlas*> atlasList;
    ControlMesh controlMesh;
    Garment* parent;

    
    
    // Create a single control mesh from the atlases in this section
    void populateControlMesh();
    void resetControlMesh();
    
    
    
    // FIXME Add prime axis for this section as a vector, to be supplied via user sketching.
    
    void addAtlas(std::string name, GarmentAtlas* pAtlas);
    
    // Don't call with a name which doesn't exist or a blank entry is created
    // returns null pointer if not found
    GarmentAtlas* getAtlas(std::string name) {
      
      std::map<std::string, GarmentAtlas*>::iterator p = atlasMap.find( name );
      if(p == atlasMap.end()) {
        return (GarmentAtlas*) 0;
      }
      return (*p).second;
    }
    
    GarmentAtlas* getAtlas(int num) {
      return atlasList.at(num);
    }
    
  private:
    ControlMesh _origControlMesh;

    
    

};

#endif // GARMENT_SECTION_H
