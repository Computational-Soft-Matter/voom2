#include <iostream>
#include "GaugeLipid.h"

namespace voom 
{
  void GaugeLipid::compute(Shellresults & R, const ShellGeometry & defGeom)
  {

    typedef vector<Vector3d> BasisVectors;

    const BasisVectors & basis = defGeom.a();
    const BasisVectors & dual = defGeom.aDual();
    const BasisVectors & dPartials = defGeom.dPartials();
    _H = - 0.5 * ( dual[0].dot( dPartials[0] ) + \
		   dual[1].dot( dPartials[1]) );
    
    const Matrix2d & aInv = defGeom.metricTensorInverse();
    const Matrix2d & G_ab = _referenceGeometry.metricTensor();
    Matrix2d mixedCurvatureTensor;
    // b^alpha_beta = mixedCurvatureTensor(alpha,beta)
    // first index is upper, second is lower
    mixedCurvatureTensor(0,0) = -dual[0].dot( dPartials[0] );
    mixedCurvatureTensor(0,1) = -dual[0].dot( dPartials[1] );
    mixedCurvatureTensor(1,0) = -dual[1].dot( dPartials[0] );
    mixedCurvatureTensor(1,1) = -dual[1].dot( dPartials[1] );

    _K =  (  mixedCurvatureTensor(0,0) * mixedCurvatureTensor(1,1) 
	     - mixedCurvatureTensor(0,1) * mixedCurvatureTensor(1,0) );

    //double twoHminusC0 = 2.0*_H;
    const Matrix2d G_abgag = (defGeom.metricTensorInverse()*_referenceGeometry.metricTensor());
    double gabG_ab = G_abgag.trace();
    double lambda = _lambda;

    //cout << twoHminusC0 << endl;

    R.Wg = lambda*gabG_ab;
    if( R.request & ENERGY ) {
      // compute strain energy
      R.W =   _kC * _H * _H;
      //R.W += lambda*gabG_ab;
    }
    if( R.request & FORCE ) {
      // -------------- stress resultants --------------
     
      R.n[0] = _kC * _H * 
	( aInv(0,0) * dPartials[0] + 
	  aInv(0,1) * dPartials[1] )
	+ (_kC * _H * _H )* dual[0];
      R.ng[0] = - 2*lambda*(G_abgag(0,0)*dual[0]+G_abgag(0,1)*dual[1]) + lambda*gabG_ab * dual[0];

      R.n[1] =  _kC * _H * 
	( aInv(1,0) * dPartials[0] + 
	  aInv(1,1) * dPartials[1] )
	+ (_kC * _H * _H )* dual[1];
      R.ng[1] = - 2*lambda*(G_abgag(1,0)*dual[0]+G_abgag(1,1)*dual[1]) + lambda*gabG_ab * dual[1];

      // ------------------- moment resultants ------------------
      R.m[0] = - _kC * _H * dual[0];
      R.m[1] = - _kC * _H * dual[1];
      //cout << "HERE\n";
    
    }
    return;
  }
}


