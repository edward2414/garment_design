//
//  Filename         : pattern.cpp
//  Author           : Jamie Wither
//  Purpose          : A viewer for the developable pattern
//  Date of creation : 08 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <vector>
#include <lib3ds/mesh.h>
#include <lib3ds/vector.h>
#include <lib3ds/quat.h>
#include <iostream>
#include <fstream>
#include <qfiledialog.h>
#include <qcursor.h>
#include <qstring.h>
#include "mdi_window.h"
#include "config.h"
#include "lcmodel.h"
#include "repository_pretreatment_result.h"
//#include "repository_canvas_result.h"
//#include "pretreatment.h"
//#include "canvas.h"
//#include "layer.h"
//#include "garment_map.h"
#include "point_utils.h"
#include "utils.h"
#include "config_pattern.h"
#include "pattern.h"
#include "EdgeSection.h"
#include "Parabola.h"
#include "dock_pattern.h"
#include "dock_pattern_widget.h"
#include "qlistview.h"
#include "pretreatment.h"
#include "math.h"
#include "Garment.h"
#include "GarmentSection.h"
#include "buckling.h"
#include "distance_field.h"
#include "stdlib.h"
#include "PushBackParams.h"
#include <qcolor.h>
 

using namespace qglviewer;

extern int other_model_dl;

int  Pattern::_clothes_horse_dl;
bool Pattern::_clothes_horse_available;

Pattern::Pattern(MDIWindow* parent, const char* name,const QGLWidget* shareWidget, int wflags, DockPattern* _dp)
  : BaseQGLViewer(parent, name, shareWidget, wflags)
{
  _mdi = parent;
  _model_dl = 0;
  _clothes_horse_dl = 0;
  _clothes_horse_available = false;
  _lcmodel = 0;
  _garment = 0;
  _hcompress_bottom = 1.0;
  _hcompress_top = 1.0;
  _twist_degree = 0.0;
  _twist_factor = 4.0;
  _twist_enabled = false;
  _collision_enabled = false;
  _axis_select_mode = false;
  _disp_model = false;
  _garmentInitialised = false;
  _dock_pattern = _dp; // needed so we can alter the list box in the dock area
  
  
  
  updateConfig();
}

Pattern::~Pattern()
{
}

void Pattern::setClothingModelDisplayList(int dl) 
{
  _clothes_horse_dl = dl;
  _clothes_horse_available = true;
}

void Pattern::updateConfig()
{
  makeCurrent();

  bool type_changed = false;
  bool atlas_changed = false;
  bool redraw = false;
  bool axis_align = false;
  
  
  int atlas_cols = 2;
  int atlas_rows = 2;
  
  int left, right, top, bottom;
  
  _disp_pattern = config::PATTERN_DISP_PATTERN;
  QString mName;
  


  
  
  // FIXME make this conditional on toggle box having changed - otherwise recompiling everytime
  
  
  if (_mdi && _mdi->config()) {
    _mdi->config()->io()->getValue("display/pattern", _disp_pattern);
    _mdi->config()->io()->getValue("display/type_changed", type_changed);
    _mdi->config()->io()->getValue("display/atlas_changed", atlas_changed);
    _mdi->config()->io()->getValue("display/control_mesh", _disp_control_mesh);
    _mdi->config()->io()->getValue("display/model_mesh", _disp_model_mesh);
    _mdi->config()->io()->getValue("display/redraw", redraw);
    _mdi->config()->io()->getValue("display/twist_factor_val", _twist_factor);
    _mdi->config()->io()->getValue("display/twist_degree_val", _twist_degree);
    _mdi->config()->io()->getValue("display/atlas_cols", atlas_cols);
    _mdi->config()->io()->getValue("display/atlas_rows", atlas_rows);
    _mdi->config()->io()->getValue("display/buckling_surface", _disp_buckling_surface);
    _mdi->config()->io()->getValue("display/twist_enabled", _twist_enabled);
    _mdi->config()->io()->getValue("display/ptf_val", _ptf_val);
    _mdi->config()->io()->getValue("display/hcompress_bottom_val", _hcompress_bottom);
    _mdi->config()->io()->getValue("display/hcompress_top_val", _hcompress_top);
    _mdi->config()->io()->getValue("display/collision_detection_enabled", _collision_enabled);
    _mdi->config()->io()->getValue("display/model", _disp_model);
    _mdi->config()->io()->getValue("display/axis_align_clicked", axis_align);
    _mdi->config()->io()->getValue("display/left_edges", left);
    _mdi->config()->io()->getValue("display/right_edges", right);
    _mdi->config()->io()->getValue("display/top_edges", top);
    _mdi->config()->io()->getValue("display/bottom_edges", bottom);
    
    
    if(_collision_enabled)
    {
      if(!RepositoryPretreatmentCanvas::distanceField()) {
        emit message(MsgHandler::MSG_ERROR, name(), "Cannot use collision detection until distance field is initialised");
        _mdi->config()->io()->setValue("display/collision_detection_enabled", false);
        _collision_enabled = false;
      }
    }
    
    // FIXME Improve the case logic for recalculation and redrawing
    if(axis_align)
    {
      // The edges required have to have been set by the user...
      _mdi->config()->io()->setValue("display/axis_align_clicked", false);
      
      if(!_garment) return;
      if(!_garment->getAtlas( currentAtlasName )) return;
      
      _garment->getAtlas( currentAtlasName )->leftEdgesReq=left;
      _garment->getAtlas( currentAtlasName )->rightEdgesReq=right;
      _garment->getAtlas( currentAtlasName )->topEdgesReq=top;
      _garment->getAtlas( currentAtlasName )->bottomEdgesReq=bottom;
      _garment->getAtlas( currentAtlasName )->_lcmodel->correctAtlasRotationUser(
          _garment->getAtlas( currentAtlasName )->_axisPointA,
          _garment->getAtlas( currentAtlasName )->_axisPointB
                                                                               );
      _garment->getAtlas( currentAtlasName )->updateEdgeSectionInformation(); // orientations need updating
      
      // FIXME Need a separate event for recalculating edges from requirements
      _garment->getAtlas( currentAtlasName )->reduceEdgeSectionsToRequiredNumber( left + right + top + bottom );
      _garment->getAtlas( currentAtlasName )->categoriseSeamEdgeSections();
      _garment->getAtlas( currentAtlasName )->fitParabola();
      
      // hardcode the row and column counts because we don't have time for a nice GUI interface
      if (currentAtlasName == "front")
      {
        // SKIRT
        _garment->getAtlas( currentAtlasName )->populateControlMesh(SKIRT_ROWS,SKIRT_COLS);
      }
      else if(currentAtlasName == "back")
      {
        _garment->getAtlas( currentAtlasName )->populateControlMesh(SKIRT_ROWS,SKIRT_COLS);
      }
      else if(currentAtlasName =="torso_front")
      {
        // TSHIRT
        _garment->getAtlas( currentAtlasName )->populateControlMesh(TSHIRT_ARM_COLS,TSHIRT_TORSO_COLS);
      }
      else if(currentAtlasName == "torso_back")
      {
        _garment->getAtlas( currentAtlasName )->populateControlMesh(TSHIRT_ARM_COLS,TSHIRT_TORSO_COLS);
      }
      else if (currentAtlasName == "right_arm_front")
      {
        _garment->getAtlas( currentAtlasName )->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
      }
      else if (currentAtlasName == "right_arm_back")
      {
        _garment->getAtlas( currentAtlasName )->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
      }
      else if (currentAtlasName ==  "left_arm_front")
      {
        _garment->getAtlas( currentAtlasName )->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
      }
      else if (currentAtlasName == "left_arm_back")
      {
        _garment->getAtlas( currentAtlasName )->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
      }
      else
      {
        _garment->getAtlas( currentAtlasName )->populateControlMesh(5,5);
      }
      
      //_garment->getAtlas( currentAtlasName )->populateControlMesh(5,5);
      _garment->getSection( currentSectionName )->populateControlMesh(); // If possible
      _garment->populateControlMesh();
      
      
      redraw = true;
    }
    
    
    if(redraw) { // only execute this if user clicks redraw
      _mdi->config()->io()->getValue("display/atlas_name", mName);
      if(!_garment) return;
      if(!_garment->getAtlas( mName ) ) return;
      if(_garment->getAtlas( mName )->parent->isValid) {
        // recalculate Control mesh if row or col count has changed.
        /* DISPLAY UNTIL WE HAVE DECENT USER INTERACTION MODEL
        if(atlas_cols != _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_controlMesh.cols
          ||
          atlas_rows != _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_controlMesh.rows
          )
        {
          if(atlas_rows != _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_controlMesh.rows)
          { // If rows have changed, change all atlases in this GarmentSection to keep row count the same
            
            emit message(MsgHandler::MSG_NORMAL, name(), "Changing row count in control mesh for all atlases to "
                + QString::number(atlas_rows));
            
            std::vector<GarmentAtlas*>::const_iterator ai, ai_end;
            ai = _garment->getSection( currentSectionName )->atlasList.begin();
            ai_end = _garment->getSection( currentSectionName )->atlasList.end();
            for(;ai!=ai_end;ai++) 
            {
              (*ai)->populateControlMesh(atlas_rows,(*ai)->_controlMesh.cols);
            }
            _garment->getSection( currentSectionName )->populateControlMesh();
          }
          else
          { // Only cols have changed so rebuild only current atlas
            emit message(MsgHandler::MSG_NORMAL, name(), "Changing column count in control mesh for atlas: " + currentAtlasName + " to " + QString::number(atlas_cols));
            _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->populateControlMesh( 
                _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_controlMesh.rows 
                ,atlas_cols
                                                                              );
            _garment->getSection( currentSectionName )->populateControlMesh();
          }
          
          // promote the changes
          //_garment->getSection( currentSectionName )->populateControlMesh();
          // FIXME promote to whole Garment
          
          QListView* qList = ((DockPatternWidget*) dock_pattern->widget())->listView1;
          
          // FIXME populate with section names also
          std::vector<GarmentAtlas*>::const_iterator ai, ai_end;
          ai = _garment->getSection( currentSectionName )->atlasList.begin();
          ai_end = _garment->getSection( currentSectionName )->atlasList.end();
          for(;ai!=ai_end;ai++)
          {
            QString qsCurrentName = QString( (*ai)->name.c_str() );
            QListViewItem *currentItem = qList->findItem( qsCurrentName , 0 );
            currentItem->setText( 1, QString::number( (*ai)->_controlMesh.rows ));
            currentItem->setText( 2, QString::number( (*ai)->_controlMesh.cols ));
          }
        } 
        */
      }
    }
    
    if(atlas_changed)
    {
      _mdi->config()->io()->setValue("display/atlas_changed", false);
      _mdi->config()->io()->getValue("display/atlas_name", mName);
      
      if(!_garment) return;
      
      currentAtlasName = std::string( mName.ascii() ) ; // FIXME must be a better conversion method
      currentSectionName = _garment->getAtlas( currentAtlasName )->parent->name;
      
      std::cout << "Switching to: " << currentSectionName << ", " << currentAtlasName << std::endl;
    
    
      
      _lcmodel = _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_lcmodel;
      
      int current_cols = _garment->getSection( currentSectionName )->getAtlas(currentAtlasName)->_controlMesh.cols;
      int current_rows = _garment->getSection( currentSectionName )->getAtlas(currentAtlasName)->_controlMesh.rows;
      _mdi->config()->io()->setValue("display/atlas_cols", current_cols);
      _mdi->config()->io()->setValue("display/atlas_rows", current_rows);
      _dock_pattern->updateConfig();
      
      redraw = true;
    }

    
    if(type_changed || redraw) {
      if(!_lcmodel) return;
      
      _mdi->config()->io()->setValue("display/redraw", false);
      
      _model_dl = prepareModel();
      
      if(type_changed || atlas_changed || axis_align) {
        _mdi->config()->io()->setValue("display/type_changed", false);
        if(_disp_pattern) {
          _bb_size = _lcmodel->tcbbMax() - _lcmodel->tcbbMin();
          setSceneBoundingBox(_lcmodel->tcbbMin().address(), _lcmodel->tcbbMax().address());
        }
        else {
          if(!_disp_model)
          {
            // Focus on the section of the garment
            if(_lcmodel) 
            {
              _bb_size = _lcmodel->bbMax() - _lcmodel->bbMin();
              setSceneBoundingBox(_lcmodel->bbMin().address(), _lcmodel->bbMax().address());
            }
          }
          else 
          {
            // unless we are displaying the whole doll also
          
            _bb_size = RepositoryPretreatmentResult::modelBBSize();
            setSceneBoundingBox(RepositoryPretreatmentResult::modelBBMin().address(), RepositoryPretreatmentResult::modelBBMax().address());
          }
        }
        camera()->showEntireScene();
      }
      updateGL();
    }
  }
}

void Pattern::init()
{
  glEnable(GL_BLEND);
#ifdef GL_FUNC_ADD    // Phil
  glBlendEquation(GL_FUNC_ADD);
#endif
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_NORMALIZE);
  glCullFace(GL_BACK);
  
  // FIXME Removed so we can share the same GL context as pretreatment
//JDW SharedContext//  setManipulatedFrame(new ManipulatedFrame());
  
}


void Pattern::draw()
{
  //if (_disp_pattern)
    renderModel();
    /*
    if(_garmentInitialised && _disp_pattern)
    {
      // Render the central axis lines supplied by the user
      glDisable(GL_DEPTH_TEST);
      glPointSize(15.0f);
      glColor3f(0.5, 0.9, 0.5);
      glBegin(GL_POINTS);
      GarmentAtlas * pGA = _garment->getAtlas( currentAtlasName );
      glVertex3dv( pGA->_axisPointA.address() );
      glColor3f(0.5, 0.5, 0.9);
      glVertex3dv( pGA->_axisPointB.address() );
      glEnd();
      glColor3f(1.0, 1.0, 1.0);
      glBegin(GL_LINES);
      glVertex3dv( pGA->_axisPointA.address() );
      glVertex3dv( pGA->_axisPointB.address() );
      glEnd();
      glEnable(GL_DEPTH_TEST);
    }
    */

}

void Pattern::renderModel()
{
  setBackgroundColor(QColor("white"));
  //glMultMatrixd(RepositoryPretreatmentResult::matrix().address());
  
  glDisable(GL_LIGHTING);
  
  if(!_disp_pattern && _disp_model)
  {
    if(_clothes_horse_available)
    {
	  glEnable(GL_LIGHTING);
	  glEnable(GL_COLOR_MATERIAL);
	  glColor3d(0.4,0.4,0.4);
      glCallList(_clothes_horse_dl);
	  glDisable(GL_COLOR_MATERIAL);
    }
  }

  glPushMatrix();
  if(!_disp_pattern) {
    glMultMatrixd(RepositoryPretreatmentResult::matrix().address());
  }
  
  if (_model_dl)
  {
    glCallList(_model_dl);
    
  }

  if(other_model_dl)
  {
    if(_garment)
    {
      if(_garment->name.find("skirt")==std::string::npos)
        // This is tshirt. Draw the skirt in the tshirt pattern if available
      {
        glPushMatrix();
        // shift back in Z
        //glTranslatef(0.0,0.0,+0.4);
        //glScalef(.9,.9,.9);
        glCallList(other_model_dl);
        glPopMatrix();
      }
    }
  }
  
  glPopMatrix();
}

void Pattern::freeModel()
{
  delete _lcmodel;
  _lcmodel = 0;
  if (_model_dl) {
    //RepositoryPretreatmentResult::setModelDisplayList(0);
    //RepositoryPretreatmentResult::pretreatmentUpdated();
    glDeleteLists(_model_dl, 1); // FIXME need to free display lists
    _model_dl = 0;
  }
}

// We are limited by the resolution of the distance field. The smaller the voxel sizes
// the closer are moved point will be to the model surface. We walk along the gradient
// of the distance field in small steps.
Vec3d Pattern::movePointOutsideModel(const Vec3d& point, void *clientData)
{
  DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
  if( df )
  {
    double step_sizeX = df->size()[0]/((double)df->sizeX());
    double step_sizeY = df->size()[1]/((double)df->sizeY());
    double step_sizeZ = df->size()[2]/((double)df->sizeZ());
    double step_size = (step_sizeY < step_sizeY) ? step_sizeY : step_sizeZ;
    step_size = (step_sizeX < step_size) ? step_sizeX : step_size;
                  
    Vec3d p = point;  
    Vec3d prevP = p;
    Vec3d grad;
    double dist, prevDist;
    point_utils::distFromZ(p[0], p[1], p[2], dist);
    prevDist = dist-1;
                
    const int MAX_ITERATIONS = 10;
    const double MIN_DIST = 1.5*step_size;
                
    int iterations=0;
    while( dist<MIN_DIST && (-prevDist)>(-dist) && iterations<MAX_ITERATIONS )
    {
      prevDist = dist;
      prevP = p;
      grad = (df->gradTrilinear(p[0],p[1],p[2])).normalizeSafe();
      p += step_size * grad;
      point_utils::distFromZ(p[0], p[1], p[2], dist);
      ++iterations;
                        //std::cout<<" "<<dist;
    }
    if( iterations>0 )
    {
                        // interpolate between the final and the previous steps
      p = prevP + (-prevDist+MIN_DIST)*(p-prevP);
                        //std::cout<<" COL_END "<<iterations<<" iter"<<std::endl;
      if( iterations>=MAX_ITERATIONS )
        std::cout<<"MAX_ITERATIONS reached"<<std::endl;
    }
    return p;
  }
  else
    return point;
}

 
int Pattern::prepareModel()
{
  // pattern selects whether to draw the pattern or the 3d version
  int model_dl = glGenLists(1);
  if (!model_dl)
    return 0;
  

  
  glNewList(model_dl, GL_COMPILE);
//  for (Lib3dsNode* p=_file->nodes; p!=0; p=p->next)
//    prepareModelRec(p);
  
  glPushMatrix();
  
  glDisable(GL_CULL_FACE);
  float no_shininess[] = {0.0, 0.0, 0.0, 0.0}; // FIXME
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, no_shininess);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glColor3f(0.9, 0.3, 0.3);
  glEnable(GL_COLOR_MATERIAL);
  

    
#ifdef BUCKLING_TEST	//Phil: defined in buckling.h
 	std::cout<<"BUCKLING TEST"<<std::endl;
	bucklingTest(_lcmodel->bbMin(), _lcmodel->bbMax());             
#else // BUCKLING_TEST

  
  if(!_disp_pattern) { // DISPLAY 3D MESH
    
    if(_disp_model_mesh)
    {
      std::vector<const Triangle*>::const_iterator it = _lcmodel->triangles().begin();
      std::vector<const Triangle*>::const_iterator itEnd = _lcmodel->triangles().end();
      
      glPointSize(3.0f);
      
      
      
      for(;it != itEnd; it++) {
        
        
        
        glBegin(GL_LINE_LOOP);
        glVertex3dv( (*it)->pointA()->address() );
        glVertex3dv( (*it)->pointB()->address() );
        glVertex3dv( (*it)->pointC()->address() );
        glEnd();
        
        
        glBegin(GL_POINTS);
        glVertex3dv( (*it)->pointA()->address() );
        glVertex3dv( (*it)->pointB()->address() );
        glVertex3dv( (*it)->pointC()->address() );
        glEnd();
        
      
      } // end of iteration over triangles in model
      

    }
    
    _garment->resetControlMesh();
    _garment->getSection( currentSectionName )->resetControlMesh();
    
    _garment->controlMesh.collisionEnabled = _collision_enabled;
    _garment->controlMesh.patchTangentFactor = _ptf_val;
    
    //_garment->getSection( currentSectionName )->controlMesh.collisionEnabled = _collision_enabled;
    //_garment->getSection( currentSectionName )->controlMesh.patchTangentFactor = _ptf_val;
    
    
    if(_twist_enabled)
    {
      Garment* pG = _garment;
      
      pG->controlMesh.twistFactor = _twist_factor;
      pG->controlMesh.twistBucklingEnabled = true;
      pG->controlMesh.twistVert(
          6.28318531*_twist_degree/360.0, 0.0 ); // FIXME add another angle option to GUI
      
      /*
      _garment->getSection( currentSectionName )->controlMesh.twistFactor = _twist_factor;
      _garment->getSection( currentSectionName )->controlMesh.twistBucklingEnabled = true;
      _garment->getSection( currentSectionName )->controlMesh.twistVert(
          6.28318531*_twist_degree/360.0, 0.0 ); // FIXME add another angle option to GUI
      */
    }
    else
    {
      Garment* pG = _garment;
      pG->controlMesh.twistBucklingEnabled = false;
      //_garment->getSection( currentSectionName )->controlMesh.twistBucklingEnabled = false;
    }
      
    _garment->getSection( currentSectionName )->controlMesh.compressHoriz(_hcompress_bottom, _hcompress_bottom, _hcompress_top, _hcompress_top);
    
    _garment->controlMesh.compressHoriz(_hcompress_bottom, _hcompress_bottom, _hcompress_top, _hcompress_top);
    
    //if(_garment->_push_back_params_right.valid)
    
    int e1,e2,e3,e4;
    e1 = e2 = e3 = e4 = false;
    int s1,s2;
    s1 =s2 =0;
    
    if (_mdi && _mdi->config()) {
      _mdi->config()->io()->getValue("display/effect1", e1);
      _mdi->config()->io()->getValue("display/effect2", e2);
      _mdi->config()->io()->getValue("display/effect3", e3);
      _mdi->config()->io()->getValue("display/effect4", e4);
      _mdi->config()->io()->getValue("display/slider1", s1);
      _mdi->config()->io()->getValue("display/slider2", s2);
    }
    
    
    if(e1)
    {
      
      PushBackParams *pright = &_garment->_push_back_params_right;
      if(pright)
      {
        std::cout << "APPLYING PUSH BACK RIGHT ARM" << std::endl;
        std::cout << "Slide1: " << s1 << std::endl;
        double len = (pright->startPoint - pright->endPoint).norm();
  
        pright->shrunkLength = (((double)s1/99.0)) * len;
        pright->a = 0.3 * len;
        pright->b = 0.4 * len;
        
        _garment->controlMesh.pushBackAlongAxis(
            pright->startPoint,
        pright->endPoint,
        pright->shrunkLength,
        pright->a,
        pright->b,
        pright->influenceRadius
                                              ); // use params
        
        // VALIDATE PARAMS
        glPointSize(7.0f);
        glBegin(GL_POINTS);
        glColor3f(1.0,1.0,1.0);
        glVertex3dv( _garment->_push_back_params_right.startPoint.address() );
        glColor3f(1.0,0.5,0.7);
        glVertex3dv( _garment->_push_back_params_right.endPoint.address() );
        glPointSize(3.0f);
        glEnd();
      }
      
    }
    
    //if(_garment->_push_back_params_left.valid)
      if(e2)
    {
      
      PushBackParams *pleft = &(_garment->_push_back_params_left);
      if(pleft)
      {
        std::cout << "APPLYING PUSH BACK" << std::endl;
        _garment->controlMesh.pushBackAlongAxis(
            pleft->startPoint,
        pleft->endPoint,
        pleft->shrunkLength,
        pleft->a,
        pleft->b,
        pleft->influenceRadius
                                              
                                              ); // use params
        
        // VALIDATE PARAMS
        glPointSize(7.0f);
        glBegin(GL_POINTS);
        glColor3f(1.0,1.0,1.0);
        glVertex3dv( _garment->_push_back_params_left.startPoint.address() );
        glColor3f(1.0,0.5,0.7);
        glVertex3dv( _garment->_push_back_params_left.endPoint.address() );
        glPointSize(3.0f);
        glEnd();
      }
    }
    
    
        // COMPRESS & TWIST TORSO
    //if(_garment->name.find("tshirt")!=std::string::npos)
    if(e3)
    {
      std::cout << "Slide2: " << s2 << std::endl;
      Vec3d twistStart( 0.5*(_lcmodel->bbMin()+_lcmodel->bbMax()) );
      twistStart[1] = 0.2*_lcmodel->bbMin()[1] + 0.8*_lcmodel->bbMax()[1];
      Vec3d twistEnd( twistStart[0], _lcmodel->bbMin()[1], twistStart[2]);
      //double ratio = (double)s2/99.0; // slider range is int [0-99] // doesn't help the effect
      double twistInfluence = 0.6*(_lcmodel->bbMax()[0]-_lcmodel->bbMin()[0]);
          //_garment->controlMesh.scaleAlongAxis(twistStart, twistEnd, 1, 0.6,    twistInfluence);
      _garment->controlMesh.twistAlongAxis(twistStart, twistEnd, 0, 0.1*M_PI,
                                          twistInfluence);
    }
    
    
    _garment->getSection( currentSectionName )->controlMesh.update( Pattern::movePointOutsideModel );
    _garment->controlMesh.update( Pattern::movePointOutsideModel );
    
    if(_disp_control_mesh)
    {
      //_garment->getAtlas( currentAtlasName )->_controlMesh.drawControlMesh();
//      std::cout << "Twist is: " << _twist << std::endl;
      //_garmentSection->getAtlas( currentAtlasName )->_controlMesh.twistVert( _twist );
      //_garment->getSection(currentSectionName)->getAtlas( currentAtlasName )->_controlMesh.drawControlMesh();
      //_garment->getSection( currentSectionName )->controlMesh.drawControlMesh();
      _garment->controlMesh.drawControlMesh();
    }
    
    if(_disp_buckling_surface)
    {
      //_garment->getSection(currentSectionName)->getAtlas( currentAtlasName )->_controlMesh.drawBucklingSurface();
      // FIXME Bizarre behaviour. Can't display buckling surface for garment in second window, but can for section!!??
      if(_garment->name.find("skirt")!=std::string::npos)
      {
        //_garment->getSection( currentSectionName )->controlMesh.drawBucklingSurface();
        _garment->controlMesh.drawBucklingSurface();
      }
      else
      {
      _garment->controlMesh.drawBucklingSurface();
      }
    }
    
    

    
    
  }
  else {  // DISPLAY THE PATTERN
    
    bool displayPatternFaces = true;
    if(displayPatternFaces)
    {
      renderPatternFaces();
    }
    
    // HIGHLIGHT OUTER EDGES
    
    
    if(false) // off while I try new method
    {
      std::vector< polymesh::Edge* >::const_iterator li =
          _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_unsortedBoundaryEdges.begin();
      std::vector< polymesh::Edge* >::const_iterator li_end =
          _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_unsortedBoundaryEdges.end();
      
      glDisable(GL_DEPTH_TEST);
      int edgeCount=0;
      for(;li!=li_end;li++)  
      {
  
        polymesh::Edge *pEdge = (*li);
        glColor3f(1.0,1.0,1.0);

        
            glBegin(GL_LINES);
         
            glVertex3dv( _lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]->address() );
            glVertex3dv( _lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]->address() );
            glEnd();
            
            // draw normal
            bool drawNormal = false;
            if(drawNormal)
            {
              glColor3f(1.0, 0.5, 1.0);
              
              glBegin(GL_LINES);
              Vec3d vectorAlongEdge( _lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]->x()
                  - _lcmodel->tex_coords()[ pEdge->_vertexIndex[0]]->x(),
              _lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]->y()
                  - _lcmodel->tex_coords()[ pEdge->_vertexIndex[0]]->y(), 0 );
              Vec3d halfWayPoint = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]) + (0.5* vectorAlongEdge);
              
              glVertex3dv( halfWayPoint.address() );
              Vec3d n = pEdge->normal();
              //std::cout << "Draw Normal: " << n.x() << ", " << n.y() << ", " << n.z() << std::endl;
              
              Vec3d normalEndPoint = halfWayPoint + n;
              glVertex3dv( normalEndPoint.address() );
              glEnd();
            }
            
        
        edgeCount++;
      }
      
    }
    
    // Courtesy of colorbrewer (on the web)
    GLubyte COLOURS[][3] = { 
      {103,0,31},
      {178,24,43},
      {214,96,77},
      {244,165,130},
      {253,219,199},
      {247,247,247},
      {209,229,240},
      {146,197,222},
      {67,147,195},
      {33,102,172},
      {5,48,97} 
    };
      
    // Draw categorized edges of pattern
    bool renderPatternEdges = true;
    if(renderPatternEdges) {
      if(_garment->getAtlas( currentAtlasName )->_seamEdgesDetermined)
      //if(false)
      {
        std::vector< EdgeSection* >::const_iterator li =
            _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_listOfSeamEdgeSections.begin();
        std::vector< EdgeSection* >::const_iterator li_end =
            _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_listOfSeamEdgeSections.end();
        
        glDisable(GL_DEPTH_TEST);
        int edgeSectionCount=0;
        int numberOfEdgeSections = _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_listOfSeamEdgeSections.size();
        //std::cout << "Drawing: " << numberOfEdgeSections << " edge sections" << std::endl;
        for(;li!=li_end;li++)  // Boundary list iterator
        {
    
          EdgeSection *pEdgeSection = (*li);
    
          std::vector<polymesh::Edge*>::const_iterator ei =     pEdgeSection->_edges.begin();
    
          std::vector<polymesh::Edge*>::const_iterator ei_end = pEdgeSection->_edges.end();
    
  
          for(;ei!=ei_end;ei++) { // Edge iterator
  
            /*
            switch(pEdgeSection->_orientation)
            {
              case EdgeSection::UP:
                glColor3ubv(COLOURS[0]);
                break;
              case EdgeSection::DOWN:
                glColor3ubv(COLOURS[2]);
                break;
              case EdgeSection::LEFT:
                glColor3ubv(COLOURS[4]);
                break;
              case EdgeSection::RIGHT:
                glColor3ubv(COLOURS[6]);
                break;
            }
            //glColor3ubv(COLOURS[(edgeSectionCount)%(sizeof(COLOURS)/sizeof(COLOURS[0]))]); 
            */
            glColor3f(0.7,0.3,0.3);
    
            glLineWidth(4.0);
            glBegin(GL_LINES);
              //glVertex3dv( _lcmodel->tex_coords()[ (*ei)->_vertexIndex[0] ]->address() );
            glVertex3dv( _lcmodel->tex_coords()[ (*ei)->_vertexIndex[0] ]->address() );
            glVertex3dv( _lcmodel->tex_coords()[ (*ei)->_vertexIndex[1] ]->address() );
            glEnd();
            
            glLineWidth(1.0);
            
            bool markEdgeSectionBeginnings = false;
            if(markEdgeSectionBeginnings)
            {
              if(ei == pEdgeSection->_edges.begin()){
                glPointSize(5.0);
                glColor3f(1.0,1.0,1.0);
                glBegin(GL_POINTS);
                glVertex3dv( _lcmodel->tex_coords()[ (*ei)->_vertexIndex[0] ]->address() );
                //std::cout << "Begin point for [" << edgeSectionCount << "] is index: " << (*ei)->_vertexIndex[0] << std::endl;
                glPointSize(3.0);
                glEnd();
              }
            }
              
              
            // Draw normals
            /*
            glColor3f(1.0, 1.0, 1.0);
            glBegin(GL_LINES);
            Vec3d vectorAlongEdge( _lcmodel->tex_coords()[ (*ei)->_vertexIndex[1] ]->x()
            - _lcmodel->tex_coords()[ (*ei)->_vertexIndex[0] ]->x(),
            _lcmodel->tex_coords()[ (*ei)->_vertexIndex[1] ]->y()
            - _lcmodel->tex_coords()[ (*ei)->_vertexIndex[0] ]->y(), 0 );
            Vec3d halfWayPoint = *(_lcmodel->tex_coords()[ (*ei)->_vertexIndex[0] ]) + (0.5* vectorAlongEdge);
              
            glVertex3dv( halfWayPoint.address() );
            Vec3d n = (*ei)->normal();
              //std::cout << "Draw Normal: " << n.x() << ", " << n.y() << ", " << n.z() << std::endl;
              
            Vec3d normalEndPoint = halfWayPoint + n;
            glVertex3dv( normalEndPoint.address() );
            glEnd();
            */
              
              
          } // End of edge iterator
          edgeSectionCount++;
        } // End of boundary list iterator
        
        
        std::vector<polymesh::Edge*>::const_iterator dei = 
            _garment->getAtlas( currentAtlasName )->_unsortedBoundaryEdges.begin();
        for(;dei!= _garment->getAtlas( currentAtlasName )->_unsortedBoundaryEdges.end();
             dei++)
        {
          polymesh::Edge* pEdge = *dei;
          if(pEdge->isDartBoundary || pEdge->isDartEdge)
          {
            glColor3f(0.7,0.3,0.3);
      
            glLineWidth(4.0);
            glBegin(GL_LINES);
                
            glVertex3dv( _lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]->address() );
            glVertex3dv( _lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]->address() );
            glEnd();
          }
        }
        
        // Draw a simple legend to represent the orientation colour scheme
        bool drawLegend = false;
        if(drawLegend)
        {
          glPointSize(10.0);
          glBegin(GL_POINTS);
          // UP
          glColor3ubv(COLOURS[0]);
          glVertex3dv( (_lcmodel->tcbbMin() + Vec3d(0.0,0.3,0.0)).address() );
          //DOWN
          glColor3ubv(COLOURS[2]);
          glVertex3dv( (_lcmodel->tcbbMin() - Vec3d(0.0,0.3,0.0)).address() );
          // LEFT
          glColor3ubv(COLOURS[4]);
          glVertex3dv( (_lcmodel->tcbbMin() - Vec3d(0.3,0.0,0.0)).address() );
          // RIGHT
          glColor3ubv(COLOURS[6]);
          glVertex3dv( (_lcmodel->tcbbMin() + Vec3d(0.3,0.0,0.0)).address() );
          
          glPointSize(3.0);
          glEnd();
        }
        
        
      }
    }
    

    
    // HIGHLIGHT PARABOLA POINTS
    glPointSize(15.0f);
    glBegin(GL_POINTS);
    
    // FIXME don't assume number of parabola - this just for development
    for(int i = 0; i< _garment->getAtlas( currentAtlasName )->parabolas.size();i++)
    {
      glColor3ubv( COLOURS[ i%(sizeof(COLOURS)/sizeof(COLOURS[0])) ] );
      Parabola* p = _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->parabolas.at(i);
      glVertex3dv( p->fit_points[0].address() );
      glVertex3dv( p->fit_points[1].address() );
      glVertex3dv( p->fit_points[2].address() );
    }

    glEnd();
    
    /*
    // plot  parabolas
    
    Parabola* y0 = _garmentSection->getAtlas( currentAtlasName )->parabolas.at(0);
    Parabola* y1 = _garmentSection->getAtlas( currentAtlasName )->parabolas.at(1);
    double xRange = y1->fit_points[2].x() - y1->fit_points[0].x();
    glColor3f(0.9, 0.9, 0.9);
    for(int i = 0; i < 10; i++) 
    {
      double x = y1->fit_points[0].x() + i*(xRange /9.0);
      Vec3d p = Vec3d( x,
                       y1->y(x),
                       0
                     );
      glBegin(GL_POINTS);
      glVertex3dv(p.address());
      glEnd();
         
    }
    */
    bool renderDartPoints = false;
    if(renderDartPoints)
    {
      glPointSize(8.0f);
      
      // Plot dart points
      glDisable(GL_DEPTH_TEST);
      glColor3f(0.0, 1.0, 1.0);
      std::vector<long>::const_iterator vi = _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_dartPoints.begin();
      std::vector<long>::const_iterator vi_end = _garment->getSection( currentSectionName )->getAtlas( currentAtlasName )->_dartPoints.end();
      glBegin(GL_POINTS);
      
      
      for(;vi!=vi_end;vi++)
      {
        glVertex3dv(_lcmodel->tex_coords()[ *vi ]->address() );  
      }
      glEnd();
    }
    
    
    // Show the control mesh if it exists
    if(_garment->getAtlas( currentAtlasName )->_controlMeshAssigned )
    {
      const ControlMesh& acm = _garment->getAtlas( currentAtlasName )->_controlMesh;
      for(int i=0; i<acm.getNumPatchs(); ++i)
      {
        Vec3d  patchCP[9];
        Vec3d  interpPatternPoints[9];
        int   patchIsSeamCP[9];
        double patchLength0[6];
        Vec2d  patchTexCoord[9];
        int    patchTexID;
        int    patchAtlasIndices[9];
        std::string atlasName;
        acm.getPatch(i, patchCP, NULL, patchIsSeamCP, patchLength0, patchTexCoord,
                     &patchTexID, patchAtlasIndices, NULL, interpPatternPoints);
        /*
                  // estimate an espilon that will be used to weld inside control points
        double cpInsideEps =
              ((patchLength0[0]+patchLength0[1]+patchLength0[2]+patchLength0[3])/4)/16;
  
                  // estimate an espilon that will be used to weld seam control points
        double cpSeamEps =
              ((patchLength0[0]+patchLength0[1]+patchLength0[2]+patchLength0[3])/4)/4;
                                  
        controlMesh.addPatch(patchCP, cpInsideEps, patchIsSeamCP, cpSeamEps,
                            patchLength0, patchTexCoord, patchTexID, patchAtlasIndices, atlasName);
        */
        //Mark the control points.
        glPointSize(10.0f);
        glBegin(GL_POINTS);
    
    
        for(int ii = 0; ii< 9;ii++)
        {
          switch(ii)
          {
            case 0:
            case 1:
            case 2:
            case 6:
            case 7:
            case 8:
              glColor3ubv( COLOURS[ (3)%(sizeof(COLOURS)/sizeof(COLOURS[0])) ] );
              break;
            case 3:
            case 4:
            case 5:
              glColor3ubv( COLOURS[ (6)%(sizeof(COLOURS)/sizeof(COLOURS[0])) ] );
              break;
              
          }
          
          
          
            //glVertex3dv( _lcmodel->tex_coords()[ patchAtlasIndices[ii] ]->address() );
          
        }

        glEnd();
        
                // Draw edges
        glDisable(GL_DEPTH_TEST);
        const int EDGES[8][3] = { {0,1,2}, {3,4,5}, {6,7,8}, {0,3,6}, {1,4,7}, {2,5,8}, {2,4,6}, {0,4,8} }; 
        std::vector<int> edgesToRender;
        edgesToRender.push_back(6);
        edgesToRender.push_back(7);
        std::vector<int>::const_iterator p = edgesToRender.begin();
        
        glLineWidth(3.0);
        glColor3d(0.0,0.7,0.0); // green diagonals
        glBegin(GL_LINES);
        for( ;p!=edgesToRender.end();p++)
        {
          int j = *p;
          //Vec3d a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][0]]];
          Vec3d a = interpPatternPoints[EDGES[j][0]];
          glVertex3dv( a.address() );
          //a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][1]]];
          a = interpPatternPoints[EDGES[j][1]];
          glVertex3dv( a.address() );

          glVertex3dv( a.address() );
          //a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][2]]];
          a = interpPatternPoints[EDGES[j][2]];
          glVertex3dv( a.address() );
        }
        glEnd();
        
        edgesToRender.clear();
        edgesToRender.push_back(3);
        edgesToRender.push_back(5);
        
        glColor3d(0.0,0.0,0.7); // blue verticals
        glLineWidth(4.0);
        p = edgesToRender.begin();
        glBegin(GL_LINES);

        for( ;p!=edgesToRender.end();p++)
        {
          int j = *p;
          //Vec3d a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][0]]];
          Vec3d a = interpPatternPoints[EDGES[j][0]];
          glVertex3dv( a.address() );
          //a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][1]]];
          a = interpPatternPoints[EDGES[j][1]];
          glVertex3dv( a.address() );

          glVertex3dv( a.address() );
          //a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][2]]];
          a = interpPatternPoints[EDGES[j][2]];
          glVertex3dv( a.address() );
        }
        glEnd();
        
        edgesToRender.clear();
        edgesToRender.push_back(0);
        edgesToRender.push_back(1);
        edgesToRender.push_back(2);
        
        // grey dotted horizontals
        glLineStipple(3, 0xAAAA);
        glEnable(GL_LINE_STIPPLE);
        glBegin(GL_LINES);
        glColor3d(0.4,0.4,0.4);
        glLineWidth(1.0);
        
        p = edgesToRender.begin();
        for( ;p!=edgesToRender.end();p++)
        {
          int j = *p;
          //Vec3d a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][0]]];
          Vec3d a = interpPatternPoints[EDGES[j][0]];
          glVertex3dv( a.address() );
          //a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][1]]];
          a = interpPatternPoints[EDGES[j][1]];
          glVertex3dv( a.address() );

          glVertex3dv( a.address() );
          //a = *_lcmodel->tex_coords()[patchAtlasIndices[EDGES[j][2]]];
          a = interpPatternPoints[EDGES[j][2]];
          glVertex3dv( a.address() );
        }
        glEnd();
        glDisable(GL_LINE_STIPPLE);
      }
    } // end if garment atlas has control mesh
    
    
    glEnable(GL_DEPTH_TEST);
    glPointSize(3.0f);
    glLineWidth(1.0);
    
    
  }
   
#endif // BUCKLING_TEST     
  
  glPopMatrix();
  glEndList();
  
    // Put the skirt in the "other model" for displaying in the other pattern window
  if(_garment->name.find("skirt") != std::string::npos)
  {
    other_model_dl = model_dl;
  }
  
  return model_dl;
}

void Pattern::loadFile()
{
  QString path;
  if (_mdi && _mdi->config())
    _mdi->config()->io()->getValue("paths/file_obj", path);
  QString filename = QFileDialog::getOpenFileName(path, "OBJ files (*.obj);;All files (*)", this);
  
  // In case of Cancel.
  if (filename.isEmpty())
    return;

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

  // If a scene is already loaded, it has to be freed.
  
  if (_lcmodel) {
    freeModel();
    // FIXME - Free all models and GarmentAtlases in GarmentSection
//    freeBBox();
  }
  
  
  
  if(_garment) {
    delete _garment;
  }
  
  
  
  //_lcmodel = new LCModel(); // holds model for pattern or mesh we are looking at currently.
  _garment = new Garment(); // holds all the patterns/meshes
  
  // Load an OBJ model.
//  _file = lib3ds_file_load(filename.latin1());
  
 // if( !utils::importOBJAtlas(_lcmodel,filename.latin1()) ) {
  if( !utils::importOBJGarment(_garment,filename.latin1(),_dock_pattern) ) {
    emit statusBarMessage("Unable to open file \"" + filename + "\".", 2000);
    emit message(MsgHandler::MSG_ERROR, name(), "Unable to open file \"" + filename + "\".");
    QApplication::restoreOverrideCursor();
    return;
  }
  
  _garmentInitialised = true;
  
  // garmentSections consist of multiple lcModels, but we need one to display, so pick the first
  _lcmodel = _garment->getSection( 0 )->getAtlas(0)->_lcmodel;
  // FIXME temporary until we add correct section attributes and loading
  currentSectionName = _garment->getSection( 0 )->name;
  currentAtlasName = _garment->getSection( 0 )->getAtlas( 0 )->name;
  
  // FIXME HARDCODED!! Fix orientations
  // Original skirt
  GarmentAtlas* pGA = _garment->getAtlas( "skirt_front" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(-10.9715, 2.21425,0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 6 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    _garment->getAtlas( currentAtlasName )->populateControlMesh(SKIRT_ROWS,SKIRT_COLS);
  }

  pGA = _garment->getAtlas( "skirt_back" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(-12.6393, -6.14883,0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 6 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(SKIRT_ROWS,SKIRT_COLS);
  }
  // NEW SKIRT
  pGA = _garment->getAtlas( "skirt_front1" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(-10.164, 1.35229, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 6 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    _garment->getAtlas( currentAtlasName )->populateControlMesh(SKIRT_ROWS,SKIRT_COLS);
  }
  
  pGA = _garment->getAtlas( "skirt_back1" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(-11.9422, -1.418, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 6 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(SKIRT_ROWS,SKIRT_COLS);
  }
  
  // Populate torso_back first as it has no darts and so has better seam point spacing when we force welds
  pGA = _garment->getAtlas( "torso_back" );
  if(pGA)
  {
    std::cout << "POP TORSO BACK" << std::endl;
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(1.55953, 8.57744, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 6 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_COLS,TSHIRT_TORSO_COLS);
  }
  
  pGA = _garment->getAtlas( "torso_front" );
  if(pGA)
  {
    std::cout << "POP TORSO FRONT" << std::endl;
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(-0.724259, 9.44618, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 6 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_COLS,TSHIRT_TORSO_COLS);
  }

  pGA = _garment->getAtlas( "right_arm_front" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(-0.228738, 7.29674, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 4 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
  }
  pGA = _garment->getAtlas( "right_arm_back" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(0.347412, 7.3883, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 4 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
  }
  pGA = _garment->getAtlas( "left_arm_front" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(3.89566, 7.30123, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 4 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
  }
  pGA = _garment->getAtlas( "left_arm_back" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(-0.439643, 7.45078, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 4 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
  }
  
  // NEW TSHIRT tshirt_dev2.obj
  // Populate torso_back first as it has no darts and so has better seam point spacing when we force welds
  pGA = _garment->getAtlas( "torso_back1" );
  if(pGA)
  {
    std::cout << "POP TORSO BACK" << std::endl;
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(1.35759, 7.70918,0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 6 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_COLS,TSHIRT_TORSO_COLS);
  }
  
  
  pGA = _garment->getAtlas( "torso_front1" );
  if(pGA)
  {
    std::cout << "POP TORSO FRONT" << std::endl;
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(0.191, 9.3097, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 6 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_COLS,TSHIRT_TORSO_COLS);
  }

  pGA = _garment->getAtlas( "right_arm_front1" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(0.964, 7.336,0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 4 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
  }
  pGA = _garment->getAtlas( "right_arm_back1" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(-0.7017, 7.349, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 4 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
  }
  pGA = _garment->getAtlas( "left_arm_front1" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(0.013, 7.307, 0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 4 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
  }
  pGA = _garment->getAtlas( "left_arm_back1" );
  if(pGA)
  {
    pGA->_lcmodel->correctAtlasRotationUser( Vec3d(0.0383, 7.170,0) );
    pGA->updateEdgeSectionInformation(); // orientations need updating
    pGA->reduceEdgeSectionsToRequiredNumber( 4 );
    pGA->categoriseSeamEdgeSections();
    pGA->fitParabola();
    pGA->populateControlMesh(TSHIRT_ARM_ROWS,TSHIRT_ARM_COLS);
  }
  
  GarmentSection* pGS = _garment->getSection( "skirt" );
    if(pGS) pGS->populateControlMesh(); 
    
    pGS = _garment->getSection( "torso" );
    if(pGS) pGS->populateControlMesh();  
    pGS = _garment->getSection( "right_arm" );
      if(pGS) pGS->populateControlMesh(); 
      pGS = _garment->getSection( "left_arm" );
      if(pGS) pGS->populateControlMesh(); 
  
  
  _garment->populateControlMesh();
  
  
  
  if(_disp_pattern) {
    
    _bb_size = _lcmodel->tcbbMax() - _lcmodel->tcbbMin();
    setSceneBoundingBox(_lcmodel->tcbbMin().address(), _lcmodel->tcbbMax().address());
    
    
  }
  else {
    _bb_size = _lcmodel->bbMax() - _lcmodel->bbMin();
    setSceneBoundingBox(_lcmodel->bbMin().address(), _lcmodel->bbMax().address());
    
    
  }
  camera()->showEntireScene();
//JDW SharedContext//  (_bb_size.norm()); // Voluntarily too big (factor 2).
  //std::cout << "BB:" << _bb_size << std::endl;
  
  

  
//JDW SharedContext//  manipulatedFrame()->setTranslationAndRotation(sceneCenter(), Quaternion());
  
  
  _model_dl = prepareModel();
  
  if (_mdi && _mdi->config())
  {
    _mdi->config()->io()->setValue("paths/file_obj", filename);
    _mdi->config()->io()->setValue("display/atlas_changed", true);
    _mdi->config()->io()->setValue("display/atlas_name", currentAtlasName);
    
  }


  QApplication::restoreOverrideCursor();
  if (_mdi && _mdi->config())
  {
    updateConfig();
  }
}

void Pattern::keyPressEvent(QKeyEvent* e)
{
  switch (e->key())
  {
    case Qt::Key_A :
      _axis_select_mode = !_axis_select_mode;
      std::cout << "Axis mode: " << _axis_select_mode << std::endl;
      break;
    
    default:
      QGLViewer::keyPressEvent(e);
  }
}

void Pattern::mousePressEvent(QMouseEvent *e)
{
  
  if (    (e->button() == Qt::LeftButton) 
       && (e->state() == Qt::NoButton)
       && _axis_select_mode
     ) {
    

    double coord[2];
    double ratio;
    convertCoordinates(e->x(), e->y(), coord[0], coord[1], ratio);
    std::cout << "MX: " << e->x() << ", MY: " << e->y() << std::endl;
    std::cout << "X: " << coord[0] << ", Y: " << coord[1] << std::endl;
    
    
    qglviewer::Vec screen=qglviewer::Vec(e->x(), e->y(), 0.5);
    qglviewer::Vec world = camera()->unprojectedCoordinatesOf(screen);
    
    // Find intersection with XY plane
    qglviewer::Vec dir = camera()->viewDirection();
    
    double denom = dir[2];
    double t = -(world[2]) / denom;
    Vec3d origin(world[0], world[1], world[2]);
    Vec3d direction(dir[0], dir[1], dir[2]);
    _garment->getAtlas( currentAtlasName )->_axisPointA = origin + t * direction;
    
    //_clickPoint= Vec3d(world[0], world[1], world[2]);
    updateGL();
  }
  else
    BaseQGLViewer::mousePressEvent(e);
}

void Pattern::mouseMoveEvent(QMouseEvent *e)
{

  if (_axis_select_mode) {
    // Don't move camera
//    updateGL();
    return;
  }

  BaseQGLViewer::mouseMoveEvent(e);
}

void Pattern::mouseReleaseEvent(QMouseEvent *e)
{
  
  if (    (e->button() == Qt::RightButton)
           && _axis_select_mode
     ) {
    

    double coord[2];
    double ratio;
    convertCoordinates(e->x(), e->y(), coord[0], coord[1], ratio);
    std::cout << "MX: " << e->x() << ", MY: " << e->y() << std::endl;
    std::cout << "X: " << coord[0] << ", Y: " << coord[1] << std::endl;

    
    qglviewer::Vec screen=qglviewer::Vec(e->x(), e->y(), 0.5);
    qglviewer::Vec world = camera()->unprojectedCoordinatesOf(screen);
    
    // Find intersection with XY plane
    qglviewer::Vec dir = camera()->viewDirection();
    
    double denom = dir[2];
    double t = -(world[2]) / denom;
    Vec3d origin(world[0], world[1], world[2]);
    Vec3d direction(dir[0], dir[1], dir[2]);
    _garment->getAtlas( currentAtlasName )->_axisPointB = origin + t * direction;
    
    //_clickPoint= Vec3d(world[0], world[1], world[2]);
    updateGL();
     }
     else
       BaseQGLViewer::mouseReleaseEvent(e);
}

void Pattern::convertCoordinates(double x, double y, double& new_x, double& new_y, double& ratio)
{
  float camera_pos[3];
  double w, h;
  //JDW - method depreciated in later versions of QGLViewer. FIXME needs more tidying up. Should use Vec everywhere
  //camera()->getPosition(camera_pos[0], camera_pos[1], camera_pos[2]);
  qglviewer::Vec pos=qglviewer::Vec();
  pos = camera()->position();
  camera_pos[0] = pos[0];
  camera_pos[1] = pos[1];
  camera_pos[2] = pos[2];

  

  camera()->getOrthoWidthHeight(w, h);
  ratio = 2.0 * w / width();
  new_x = camera_pos[0] + w * (2.0 * x / width() - 1);
  new_y = camera_pos[1] - h * (2.0 * y / height() - 1);
  
  
}

void Pattern::renderPatternFaces()
{
        // use the texture co-ordinates to generate the display list
  LCModel::IndexedTriangleList::iterator it = _lcmodel->indexedTriangles().begin();
  LCModel::IndexedTriangleList::iterator itEnd = _lcmodel->indexedTriangles().end();
      
  glPointSize(3.0f);
      
  glColor3f(1.0,0.8,0.8);
  for(;it != itEnd; it++) {
      
  
        //std::cout << "TriIndex: " << (*it)->indexA() << std::endl;
        //std::cout << "Array size: " << _lcmodel->tex_coords().size() << std::endl;
    glBegin(GL_LINE_LOOP);
    glVertex3dv( _lcmodel->tex_coords()[ (*it)->indexA() ]->address() );
    glVertex3dv( _lcmodel->tex_coords()[ (*it)->indexB() ]->address() );
    glVertex3dv( _lcmodel->tex_coords()[ (*it)->indexC() ]->address() );
        
    glEnd();
        
  
        
        // end draw edges
        
    glBegin(GL_POINTS);
    glVertex3dv( _lcmodel->tex_coords()[ (*it)->indexA() ]->address() );
    glVertex3dv( _lcmodel->tex_coords()[ (*it)->indexB() ]->address() );
    glVertex3dv( _lcmodel->tex_coords()[ (*it)->indexC() ]->address() );
    glEnd();
  }
}

