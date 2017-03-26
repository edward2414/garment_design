//
//  Filename         : distance_field.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Distance Field.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <limits>
#include <qstring.h>
#include "lcmodel.h"
#include "octree.h"
#include "utils.h"
#include "point_utils.h"
#include "distance_field.h"
#include "chrono.h"

#include "glutils.h"





DistanceField::DistanceField(unsigned x, unsigned y, unsigned z) :
  VolumeData<valuetype>(x, y, z)
{
  _model = 0;
  _max_dist = 1;
  if (x && y)
    _type.resize(x * y);
}

DistanceField::DistanceField(DistanceField *df) :
    VolumeData<valuetype>(df->_x, df->_y, df->_z)
{
  //DistanceField *newDF = new DistanceField(df->sizeX(),df->sizeY(),df->sizeZ());
  
  // Copy distance field fields
  _max = df->_max;
  _max_dist = df->_max_dist;
  _min = df->_min;
  _model= df->_model;
  _type.assign(df->_type.begin( ), df->_type.end( ));
  
  // copy VolumeData fields
  _max_val = df->_max_val;
  _val.assign(df->_val.begin(), df->_val.end());
  _x = df->_x;
  _y = df->_y;
  _z = df->_z;
  
  
}
    

DistanceField::~DistanceField()
{

}
// JDW Added for gradient vector calculation
// Uses central differences so we don't have a valid value at the boundaries.
// As we don't care about values far from the surface we just repeat the nearest
// valid value
// Values outside the volume are clipped
Vec3d DistanceField::gradCentralDifference(unsigned i, unsigned j, unsigned k) const
{
  double xRight, xLeft, yRight, yLeft, zRight, zLeft;
  
  // handle edge effects 
  if(i<=0) i = 1; 
  if(i>=(sizeX()-1)) i = sizeX() - 2;
  if(j<=0) j = 1; 
  if(j>=(sizeY()-1)) j = sizeY() - 2;
  if(k<=0) k = 1; 
  if(k>=(sizeZ()-1)) k = sizeZ() - 2;
 
  
  xRight = dist(i + 1, j, k);
  xLeft = dist(i - 1, j, k);
  
  yRight = dist(i, j + 1, k);
  yLeft = dist(i, j - 1, k);
  
  zRight = dist(i, j, k + 1);
  zLeft  = dist(i, j, k - 1);
      
  double x = (xRight - xLeft) / 2. * (size()[0]/sizeX());
  double y = (yRight - yLeft) / 2. * (size()[1]/sizeY());
  double z = (zRight - zLeft) / 2. * (size()[2]/sizeZ());
  return Vec3d(x,y,z);

}


Vec3d DistanceField::gradTrilinear(double x, double y, double z) const
{
    x = sizeX() * (x - min()[0]) / (size()[0]) - 0.5;
    y = sizeY() * (y - min()[1]) / (size()[1]) - 0.5;
    z = sizeZ() * (z - min()[2]) / (size()[2]) - 0.5;
  
    int i = static_cast<int>(x), j = static_cast<int>(y), k = static_cast<int>(z);

    if (i < 0)
      i = 0;
    else if (i >= static_cast<int>(_x) - 1)
      i = _x - 2;
    if (j < 0)
      j = 0;
    else if (j >= static_cast<int>(_y) - 1)
      j = _y - 2;
    if (k < 0)
      k = 0;
    else if (k >= static_cast<int>(_z) - 1)
      k = _z - 2;
      
    Vec3d a11, a12, a21, a22, b11, b12, b21, b22;
    Vec3d z1, a1, a2, z2, b1, b2, final;

    a11 = gradCentralDifference(i, j, k);
    a12 = gradCentralDifference(i + 1, j, k);
    a21 = gradCentralDifference(i, j + 1, k);
    a22 = gradCentralDifference(i + 1, j + 1, k);
    b11 = gradCentralDifference(i, j, k + 1);
    b12 = gradCentralDifference(i + 1, j, k + 1);
    b21 = gradCentralDifference(i, j + 1, k + 1);
    b22 = gradCentralDifference(i + 1, j + 1, k + 1);
    
//    if (!a11 && !a12 && !a21 && !a22 && !b11 && !b12 && !b21 && !b22)
 //     return 0;
      
    // interpolate along X between Y slices in front Z slice 
      // bottom Y
    a1 = bilinearVectorInterpolation(a11, a12, x - i);
    
      // top Y
    a2 = bilinearVectorInterpolation(a21, a22, x - i);
    
      // interpolated vector for front Z slice
    z1 = bilinearVectorInterpolation(a1, a2, y - j);
    
    // interpolate along X between Y slices in back Z slice 
      // bottom Y
    b1 = bilinearVectorInterpolation(b11, b12, x - i);
    
      // top Y
    b2 = bilinearVectorInterpolation(b21, b22, x - i);
    
    // interpolated vector for back Z slice
    z2 = bilinearVectorInterpolation(b1, b2, y - j);
    
    // finally trilinearly interpolate 2 bilinear Vectors between the two Z slices
    final = bilinearVectorInterpolation(z1, z2, z - k);

    return final;
  
}

Vec3d DistanceField::bilinearVectorInterpolation(const Vec3d &v1, const Vec3d &v2, const double &alpha) const
{
  return ( (1.0-alpha)*v1 + alpha*v2 );
}

double DistanceField::actualMaxDist()
{
  return (_max - _min).norm();
}

void DistanceField::compute(const Octree *octree)
{
  if (!octree)
    return;

  _model = octree->model();
  _min = octree->min();
  _max = octree->max();

  // Voluntarily too small (factor 2), to increase precision.
  _max_dist = (_max - _min).norm() / 2;

  double dx = (_max[0] - _min[0]) / _x;
  double dy = (_max[1] - _min[1]) / _y;
  double dz = (_max[2] - _min[2]) / _z;
  double dist;
  const Triangle *tr;
  Vec3d point;
  point[0] = _min[0] + dx / 2;
  for (unsigned i = 0; i < _x; ++i) {
    point[1] = _min[1] + dy / 2;
    for (unsigned j = 0; j < _y; ++j) {
      point[2] = _min[2] + dz / 2;
      _type[i * _y + j] = point_utils::type(point[0], point[1]);
      for (unsigned k = 0; k < _z; ++k) {
	dist = _max_dist * _max_dist;
	tr = 0;
	computeRec(dist, tr, point, octree->root(), _min, _max, octree->level());
	dist = sqrt(dist);
	if (_type[i * _y + j] != Point::OUT && tr && (point - *tr->pointA()) * tr->normal() < 0)
	  dist = -dist;
	_val[k * _y * _x + i * _y + j] =
	  static_cast<valuetype>(std::numeric_limits<valuetype>::max() * dist / _max_dist);
	point[2] += dz;
      }
      point[1] += dy;
    }
    point[0] += dx;
  }
}

void DistanceField::computeRec(double& dist, const Triangle*& tr, const Vec3d& point, const OctreeNode *node, const Vec3d& min, const Vec3d& max, unsigned level)
{
  if (node->level() == level) {
    double dist_tmp;
    OctreeNode::TriangleList::const_iterator it = node->triangles().begin();
    OctreeNode::TriangleList::const_iterator it_end = node->triangles().end();
    for (; it != it_end; ++it) {
      dist_tmp = utils::squareDistPointTriangle(point, **it);
      if (dist_tmp < dist) {
	dist = dist_tmp;
	tr = *it;
      }
    }
    return;
  }

  Vec3d center((min + max) / 2);
 
  unsigned best_guess = 0;
  if (point [0] > center[0])
    best_guess += 1;
  if (point [1] > center[1])
    best_guess += 2;
  if (point [2] > center[2])
    best_guess += 4;

  unsigned i;
  Vec3d new_min, new_max;
  for (i = best_guess; i < 8; ++i) {
    if (!node->child(i))
      continue;
    new_min[0] = i & 1 ? center[0] : min[0];
    new_min[1] = i & 2 ? center[1] : min[1];
    new_min[2] = i & 4 ? center[2] : min[2];
    new_max[0] = i & 1 ? max[0] : center[0];
    new_max[1] = i & 2 ? max[1] : center[1];
    new_max[2] = i & 4 ? max[2] : center[2];
    if (utils::squareDistPointBox(point, new_min, new_max) >= dist)
      continue;
    computeRec(dist, tr, point, node->child(i), new_min, new_max, level);
  }
  for (i = 0; i < best_guess; ++i) {
    if (!node->child(i))
      continue;
    new_min[0] = i & 1 ? center[0] : min[0];
    new_min[1] = i & 2 ? center[1] : min[1];
    new_min[2] = i & 4 ? center[2] : min[2];
    new_max[0] = i & 1 ? max[0] : center[0];
    new_max[1] = i & 2 ? max[1] : center[1];
    new_max[2] = i & 4 ? max[2] : center[2];
    if (utils::squareDistPointBox(point, new_min, new_max) >= dist)
      continue;
    computeRec(dist, tr, point, node->child(i), new_min, new_max, level);
  }
}

const int jdwWIDTH=400;
const int jdwHEIGHT=400;
float jdwimage[jdwWIDTH*jdwHEIGHT*4]; // RGB

//std::vector<float> image;

bool dfsliceready=false;
double dMin = 9999999.0;
double dMax = -9999999.0;

// write a image to the pixel buffer
void DistanceField::drawSlice() {
 if(dfsliceready)
 {
    std::cout << "Drawing slice" << std::endl;
    //glDrawPixels(jdwWIDTH,jdwHEIGHT,GL_GREEN,GL_FLOAT,jdwimage);
    glDrawPixels(jdwWIDTH,jdwHEIGHT,GL_RGBA,GL_FLOAT,jdwimage);
 }
}
	
void DistanceField::initSlice(DistanceField *df, float percentage) {
	
    
	double dx = df->size()[0] / jdwWIDTH;
	double dz = df->size()[2] / jdwHEIGHT;
    double y = df->min()[1] + (df->size()[1] * percentage);
    
    std::cout << "Per: " << percentage << std::endl;
    
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
    
    // setup colours
    const int nbCol = 20;
    float rgbaColourArray[nbCol*4];
    for(int i=0; i<nbCol*4;i+=4)
    {
      //rgbaColourArray[i] = static_cast<float>( (nbCol-1) - (i/4.0) ) / static_cast<float>(nbCol-1);
      
      
      //float green = static_cast<float>( (nbCol-1) - (i/8.0) ) / static_cast<float>(nbCol-1);
      // clamp
      //if(green > 1.0) green = 1.0;
      //rgbaColourArray[i+1] = green;
      
      //rgbaColourArray[i+2] = static_cast<float>(i/4.0) / static_cast<float>(nbCol-1);
      
      rgbaColourArray[i] = COLOURS[(i/4)%12][0]/255.0;
      rgbaColourArray[i+1] = COLOURS[(i/4)%12][1]/255.0;
      rgbaColourArray[i+2] = COLOURS[(i/4)%12][2]/255.0;
      
      // alpha
      rgbaColourArray[i+3] = 1.0;
    }
      
    
    
    //image.clear();
    //image.resize(jdwHEIGHT * jdwWIDTH * 4);
    
	
    // first determine range of values
    for (int j=0;j<jdwHEIGHT;j++)
	{
        for(int i=0;i<jdwWIDTH;i++) 
		{
            jdwimage[j*jdwWIDTH + i] =0.0;
            //image[j*jdwWIDTH + i] =0.0;
            double d;
            point_utils::distFromZ(df->min()[0] + i*dx, y, df->max()[2] - j*dz, d);
            //point_utils::distNormalisedTo1stOrder(df->min()[0] + i*dx, y, df->max()[2] - j*dz, d);
			if(d<dMin) dMin = d;
			if(d>dMax) dMax = d;
		}
	}
 
    std::cout << "dMin: " << dMin << std::endl;
    std::cout << "dMax: " << dMax << std::endl;
    std::cout << "DFMax: " << df->maxDist() << std::endl;
    std::cout << "dx: " << dx << std::endl;
    std::cout << "dz: " << dz << std::endl;
    std::cout << "dfSize: " << df->size() << std::endl;
    
    std::cout << "dfmin: " << df->min() << std::endl;
    std::cout << "dfmax: " << df->max() << std::endl;
    
    // now render slice
    for (int j=0;j<jdwHEIGHT;j++)
    {
        for(int i=0;i<jdwWIDTH;i++) 
        {
          double d;
          point_utils::distFromZ(df->min()[0] + i*dx, y, df->max()[2] - j*dz, d);
          //point_utils::distNormalisedTo1stOrder(df->min()[0] + i*dx, y, df->max()[2] - j*dz, d);
          //jdwimage[i*jdwHEIGHT + j] = 1.0 - (d / dMax);
          

          
          // group as bands of values
          //double lum = (d-dMin) / dMax; // 0 to 1
          double lum = d;
          if(d <= 0) d = 0; // not interested in interior field for now
          
          lum = d / dMax; // 0 to 1 for euclidean dist - brighter at larger dist
          
          //double lum = d; // 0 to 1 - for normalised - brighter at larger dist
          lum = lum * (nbCol-1);
          //lum = floor(lum)/20.0;
          
          int index = static_cast<int>(lum);
          //int surfaceIndex = static_cast<int>( (0-dMin / dMax) * (nbCol-1) );
          int surfaceIndex = 0;

          // Make surface and interior black
          if(index == surfaceIndex && d <= (0.001 * dMax)) 
          {
            // closest to surface in black
            jdwimage[j*jdwWIDTH*4 + 4*i] = 0.0;
            jdwimage[j*jdwWIDTH*4 + 4*i +1] = 0.0;
            jdwimage[j*jdwWIDTH*4 + 4*i +2] = 0.0;
            jdwimage[j*jdwWIDTH*4 + 4*i +3] = 1.0;
            continue; 
          }
          
          //index = i%nbCol;
          
          //jdwimage[j*jdwWIDTH + i] = lum; 
         
          jdwimage[j*jdwWIDTH*4 + 4*i] = rgbaColourArray[index*4]; 
          jdwimage[j*jdwWIDTH*4 + 4*i +1] = rgbaColourArray[index*4+1]; 
          jdwimage[j*jdwWIDTH*4 + 4*i +2] = rgbaColourArray[index*4+2]; 
          jdwimage[j*jdwWIDTH*4 + 4*i +3] = rgbaColourArray[index*4+3]; 
          
          
          
        }
    }
	dfsliceready = true;
}

// can be optimized further
// produces new normalised distance field given a distance field.
DistanceField* giveNormalisedField(DistanceField *df, VecData<Vec3d> *gf) 
{
  // indexing outside the field is clamped back to the range
  DistanceField *pNormDF = new DistanceField(df);
  
  for(unsigned i=0; i<df->sizeX(); i++)
    for(unsigned j=0; j<df->sizeY(); j++)
      for(unsigned k=0; k<df->sizeZ(); k++)
  {

    double d = df->dist(i,j,k);
    
    //grad_d = df->gradTrilinear(i,j,k);
    Vec3d grad_d = gf->getVal(i,j,k);
    //grad_d.normalizeSafe();
    
    double grad_d_sq = grad_d[0]*grad_d[0]+
          grad_d[1]*grad_d[1]+
          grad_d[2]*grad_d[2];
    
    double d_sq = d * d;

    // this value between -1 and 1
    double normalised = d / sqrt(d_sq + grad_d_sq);
    
    // scale back into distance field numeric range.
    pNormDF->_val[k * pNormDF->sizeY() * pNormDF->sizeX() + i * pNormDF->sizeY() + j] = static_cast<DistanceField::valuetype>(
          std::numeric_limits<DistanceField::valuetype>::max() * normalised
         );
  }
  
  return pNormDF;
}

VecData<Vec3d>* DistanceField::giveGradientField(DistanceField *df)
{
  VecData<Vec3d> *pGradField = new VecData<Vec3d>(df->sizeX(), df->sizeY(), df->sizeZ());
  for(unsigned i=0; i<df->sizeX(); i++)
    for(unsigned j=0; j<df->sizeY(); j++)
      for(unsigned k=0; k<df->sizeZ(); k++)
  {

    pGradField->setVal(i,j,k, df->gradTrilinear(i,j,k).normalizeSafe() );
    //pGradField->setVal(i,j,k, df->gradCentralDifference(i,j,k)); //.normalizeSafe() );

  }
  return pGradField;
}

DistanceField *DistanceField::give2NormalField(DistanceField *df, VecData<Vec3d> *gf)
{
  Chronometer chrono;
  chrono.start();
  DistanceField *norm2df = new DistanceField(df);
  std::cout << "DF copy in: " << QString::number(chrono.stop()) << " secs" << std::endl;
  chrono.start();
      
  DistanceField *norm1df = giveNormalisedField(df,gf);
  std::cout << "norm1 in: " << QString::number(chrono.stop()) << " secs" << std::endl; chrono.start();
  DistanceField *dirDiff1df = giveDirectionallyDifferentiatedField(norm1df,gf);
  
  //DistanceField *norm1df = df; // don't normalise the field
  //DistanceField *dirDiff1df = giveDirectionallyDifferentiatedField(df,gf);
  
  std::cout << "dirdiff1 in: " << chrono.stop() << " secs" << std::endl; chrono.start();
  DistanceField *dirDiff2df = giveDirectionallyDifferentiatedField(dirDiff1df,gf);
  std::cout << "dirdiff2 in: " << chrono.stop() << " secs" << std::endl; chrono.start();
  for(unsigned i=0; i<df->sizeX(); i++)
    for(unsigned j=0; j<df->sizeY(); j++)
      for(unsigned k=0; k<df->sizeZ(); k++)
  {
    norm2df->val(i,j,k) =
        norm1df->val(i,j,k) -
        (0.5 * norm1df->val(i,j,k) * norm1df->val(i,j,k)) *
        dirDiff2df->val(i,j,k);
  }
  std::cout << "norm2 in: " << chrono.stop() << " secs" << std::endl; 
  return norm2df;
}

DistanceField *DistanceField::giveXNormalField(DistanceField *df, VecData<Vec3d> *gf, int order)
{
  if(order <= 0)
  {
    std::cerr << "Error: shouldn't call giveXNormalField for order less than 2" << std::endl;
    return df; // shouldn't be called for order less than 2
  }
 
  DistanceField *firstOrderField;

  if(order == 1) {
    firstOrderField = giveNormalisedField(df,gf);
    return firstOrderField;
  }
  
  DistanceField *orderMinus1Field = giveXNormalField( df,  gf, order-1);
  DistanceField *toDifferentiate = orderMinus1Field;
  DistanceField *result = new DistanceField(df);
  
  int count = 2;
  while(count<=order)
  {
    
    for(int m=1; m<=order;m++) 
    {
      toDifferentiate = giveDirectionallyDifferentiatedField(toDifferentiate,gf);
    }
    
    // order must be available at compile time
    double denom = utils::factorial(order);
    double coeff = 1.0 / static_cast<double>(denom);
  
    for(unsigned i=0; i<df->sizeX(); i++)
        for(unsigned j=0; j<df->sizeY(); j++)
        for(unsigned k=0; k<df->sizeZ(); k++)
    {
      double power_term = orderMinus1Field->val(i,j,k);
      for(int m=2; m<=order;m++) {
        power_term *= power_term;
      }
      
      result->val(i,j,k) =
          orderMinus1Field->val(i,j,k) -
          (coeff * power_term) *
          toDifferentiate->val(i,j,k);
    }
    count++;
  }
  delete orderMinus1Field;orderMinus1Field=0;toDifferentiate=0;
  return result;
}

DistanceField *DistanceField::giveDirectionallyDifferentiatedField
    (DistanceField *df, VecData<Vec3d> *gf)
{
  // ( f(x+ha, y+hb, z+hc) - f(x-ha, y-hb, z-hc) ) / 2h
  // where
  //  h is unit distance 
  //  a, b, c are the co-efficients for the supplied direction vector ai + bj + ck
  
  double left, right;
  DistanceField* newdf = new DistanceField(df);
  
  for(unsigned i=0; i<df->sizeX(); i++)
  {
        for(unsigned j=0; j<df->sizeY(); j++)
        {
            for(unsigned k=0; k<df->sizeZ(); k++)
            {

                
                Vec3d normal = gf->getVal(i,j,k); // unit normal
                
                
                double dx = i + normal.x();
                double dy = j + normal.y();
                double dz = k + normal.z();
                
                // handle edges
                
                if(dx < 0.0) { dx = 0; }
                else if(dx > df->sizeX()-1) { dx = df->sizeX() - 1; };
                
                if(dy < 0.) { dy = 0; }
                else if(dy > df->sizeY()-1) { dy = df->sizeY() - 1; };
                
                if(dz < 0.) { dz = 0; }
                else if(dz > df->sizeZ()-1) { dz = df->sizeZ() - 1; };
                
                left = df->dist(dx,dy,dz);
                
                dx = i - normal.x();
                dy = j - normal.y();
                dz = k - normal.z();
                if(dx < 0.) { dx = 0; }
                else if(dx > df->sizeX()-1) { dx = df->sizeX() - 1; };
                if(dy < 0.) { dy = 0; }
                else if(dy > df->sizeY()-1) { dy = df->sizeY() - 1; };
                if(dz < 0.) { dz = 0; }
                else if(dz > df->sizeZ()-1) { dz = df->sizeZ() - 1; };
                
                right = df->dist(dx,dy,dz);
                
                newdf->val(i,j,k) = (left - right) / 2.0; // divide by twice h
            
            }
        }
  }
  return newdf;
}

// print out some field statistics
void DistanceField::dumpStats() 
{
  double actualMax=numeric_limits<double>::min();
  double actualMin=numeric_limits<double>::max();
  double d=0.0;
  for(unsigned i=0; i<sizeX(); i++)
   for(unsigned j=0; j<sizeY(); j++)
    for(unsigned k=0; k<sizeZ(); k++)
      {
         d = dist(i,j,k);
         if(d<actualMin) actualMin = d;
         if(d>actualMax) actualMax = d;
      }

  
  std::cout << "Field Actual Max: " << actualMax << std::endl;
  std::cout << "Field Actual Min: " << actualMin << std::endl;
  

}
