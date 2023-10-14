#include "BiHarmonic.h"
#include <iostream>

namespace voom 
{
  void BiHarmonic::compute(Scalarresults & R, 
			   Vector2d uv, Vector2d uPartials, Vector2d vPartials,
			   Matrix2d u2Partials, Matrix2d v2Partials, ShellGeometry & geom)
  {
    Matrix2d gab = geom.metricTensorInverse();
    Real LapU = (gab*u2Partials).trace() - uPartials.dot(geom.gklGammaI_kl());
    Real LapV = (gab*v2Partials).trace() - vPartials.dot(geom.gklGammaI_kl());
    
    //if (abs(LapPsi/psi+20) > 3) cout <<"LapPsi/psi : " << LapPsi/psi << endl << "---------------------" << endl;
    if( R.request & ENERGY ) {
      
      R.W = 0;// 0.5*(LapU*LapU +LapV*LapV);
      
    }
		
    if( R.request & FORCE ) {
      
      R.n1 = 0;; //LapU; //Coefficient of Lap(delta_u)

      R.n2 << 0,0;//LapV;
      
    return;
  }
}

}
