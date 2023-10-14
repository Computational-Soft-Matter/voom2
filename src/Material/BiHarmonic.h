#if !defined(__LB_h__)
#define __LB_h__

#include "VoomMath.h"
#include "ShellGeometry.h"

namespace voom
{
  
  class BiHarmonic
  {

  public:
    struct Scalarresults {
      Scalarresults(){
	W = 0.0;
	n1 = 0;;
	n2 << 0, 0;
	request = 0;
      };
      
      double W; //Energy Density
      double n1; //coefficient of variation in conc
      Vector2d n2; //coefficient of variation in conc_alpha
      double n3; //coefficient of variation in Lap(conc)
      int request;
    };

  BiHarmonic(){;};
  ~BiHarmonic(){;};
  BiHarmonic(int MatID): _matID(MatID) {;}
  void compute(Scalarresults & R,  
	       Vector2d uv, Vector2d uPartials, Vector2d vPartials,
	       Matrix2d u2Partials, Matrix2d v2Partials, ShellGeometry& geom);
  protected:

    int _matID;


  };
  
} //namespace voom

#endif 
