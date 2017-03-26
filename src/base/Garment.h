//
// C++ Interface: Garment
//
// Description: 
// An object describing the whole garment. Consisting of multiple GarmentSections, which contain Garment Atlases.
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef  GARMENT_H
# define GARMENT_H

#include <cassert>
#include <string>
#include <map>
#include <vector>
#include "PushBackParams.h"

#include "buckling.h"

class GarmentSection;
class GarmentAtlas;


class Garment {


  public:
    
    Garment(void) {
      _push_back_params_left.valid = false;
      _push_back_params_right.valid = false;
    
    };
    ~Garment(void);
    
    std::string name;
    std::map<std::string, GarmentSection*> sectionMap;
    std::vector<GarmentSection*> sectionList;
    ControlMesh controlMesh;
    
    // Create a single control mesh from the GarmentSections in this Garment
    void populateControlMesh();
    void resetControlMesh();
    void determineCoincidentPoints(); // only call once all GarmentAtlases are populated
    

    void addSection(std::string name, GarmentSection* pSection);
    
    GarmentSection* getSection(std::string name) {
      std::map<std::string, GarmentSection*>::iterator p = sectionMap.find( name );
      if(p == sectionMap.end()) {
        return (GarmentSection*) 0;
      }
      return (*p).second;
    }
    
    // As atlases have unique names we can search through hierarchy for correct one
    GarmentAtlas* getAtlas(std::string name);
    
    GarmentSection* getSection(int num) {
      assert(num <= sectionList.size());
      return sectionList.at(num);
    }
    
    void printHierarchy();
    PushBackParams _push_back_params_right;
    PushBackParams _push_back_params_left;
    
  private:
    ControlMesh _origControlMesh;
};


#endif // GARMENT_H
