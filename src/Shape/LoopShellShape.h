//-*-C++-*_
/*!
  \file LoopShellShape.h
  
  \brief Shape function implementation for Loop Shell Element
*/

#if !defined(__LoopShellShape_h__)
#define __LoopShellShape_h__

#include "Shape.h"
#include <stdio.h>
#include <stdlib.h>

#define IX(i,j,n) ((i)+(n)*(j))

typedef
struct
{
  double * val;
  double * vecI;
  double ** Phi;
} EigValStruct;

namespace voom{

  class LoopShellShape: public Shape {

  public:
    typedef Eigen::Matrix<double,Dynamic,Dynamic> SubdivisionMatrix;
    typedef vector<int> CornerValences; //Valences of nodes at three corners

    //! LoopShellShape constructor fills in N, DN and DDN
    LoopShellShape(int nodes, CornerValences & V, const VectorXd & Point) {
      read_eval(); // read Stam's file to deal with valences<=50
      _nodes = nodes;
      _coords = Point;
      _Valences = V;
      _N.resize(nodes, 0.0); //12 shape functions
      _DN.resize(nodes, Vector2d::Zero()); //partial derivatives in 'v' and 'w' directions
      _DDN.resize(nodes, Matrix2d::Zero()); //partial deriatives in ['vv', 'vw'; 'wv' 'ww'] directions
      update(Point);
      //Vector2d temp;
      //temp << 0.6, 0.2; 
      //checkConsistency(temp, 1e-7,1e-8);
      //checkConsistency2(temp, 1e-7,1e-8);
    }

    ~LoopShellShape(){
      _N.resize(0);
      _DN.resize(0);
      _DDN.resize(0);
      _N.shrink_to_fit();
      _DN.shrink_to_fit();
      _DDN.shrink_to_fit();
      free_memory();
    }


      void free_memory ()
  {

    FILE * f;
    int Nmax, N, K;

    if ( !(f=fopen("/Users/sanjay/Google Drive/bucknell/research/projects/FEM_shell/code/voom2/src/Shape/lpdata50NT.dat","rb")) ) {cout<< "Error. lpdata50NT.dat not found" << endl; exit(-1); };
  
    fread ( &Nmax, sizeof(int), 1, f );
    //Nmax = 6;
  
    for (int i=0 ; i<Nmax-2 ; i++ )
      {
	N = i+3;
	K = N+6;

	free(_eigen[i]->val);
	free(_eigen[i]->vecI);
	free(_eigen[i]->Phi[0]);
	free(_eigen[i]->Phi[1]);
	free(_eigen[i]->Phi[2]);

	free(_eigen[i]->Phi);
	free(_eigen[i]);

      }

    free(_eigen);
    fclose ( f );

    //*pNmax = Nmax;

  }
    //! Update recomputes N, DN and DDN at a Point assuming it's in a regular patch
    void update(const VectorXd & Point);

    //! Get number of shape functions
    uint getShapeFunctionNum() {return _nodes; };

    //! getN returns shape function values at a given node a.
    Real getN(const uint a) {
      return _N[a];
    };
    
    //! GetDN returns shape function derivatives at node a, in direction 'v'/'w'.
    Real getDN(const uint a, const uint i) {
      return _DN[a](i);
    };
   
    //! GetDDN returns shape function second derivatives at node a
    //! DDN is a vector of Matrix2d 
    Real getDDN(const uint a, const uint i, const uint j) {
      return _DDN[a](i,j);
    };
	
    bool checkConsistency2(Vector2d Point, const Real eps, const Real tol)
  {
    // Test dimension (1D, 2D, 3D)
    const uint dim = Point.size();
    
    const uint NumF = this->getShapeFunctionNum();
    vector<Vector2d > DN(NumF, Vector2d::Zero());
    vector<Matrix2d > DDNnumerical(NumF, Matrix2d::Zero());
   
    Real error = 0.0, norm = 0.0;

    // Compute numerical derivatives
    for(uint j = 0; j < 2; j++){
      for(uint i = 0; i < 2; i++)
	{
	  Point(i) += eps;
	  this->update(Point);
	  for(uint a = 0; a < NumF; a++) 
	    DDNnumerical[a](i,j) = this->getDN(a,j);
	  
	  Point(i) -= 2*eps;
	  this->update(Point);
	  for(uint a = 0; a < NumF; a++)
	    DDNnumerical[a](i,j)-= this->getDN(a,j);
	  
	  Point(i) += eps;
	  for(uint a = 0; a < NumF; a++) 
	    DDNnumerical[a](i,j) /= (2.0*eps);    
	}
    }
    // Compute error norm
    for(uint a = 0; a < NumF; a++) {
      for(uint j = 0; j < dim; j++) {
	for(uint i = 0; i < dim; i++) {
	  error += pow(this->getDDN(a,i,j) - DDNnumerical[a](i,j), 2);
	  norm  += pow(this->getDDN(a,i,j), 2);
	}
      }
    }
    norm = sqrt(norm);
    error = sqrt(error);
    
    cout << "Error = " << error << " Norm = " << norm << endl;
    if ( abs(error) < norm * tol) {
      cout << "Shape consistency check passed" << endl;
      return true;
    }
    else {
      cout << "Shape consistency check failed" << endl
	   << setw(10) << "n"
	   << setw(24) << "analytical"
	   << setw(24) << "numerical" << endl;   
      for(uint a = 0; a < NumF; a++)
	for(uint i = 0; i < dim; i++){
	  for(uint j = 0; j < dim; j++){
	    cout << setw(10) << a << setw(24) << this->getDDN(a,i,j)
		 << setw(24) << DDNnumerical[a](i,j) ;
	  }
	  cout << endl;
	}
    }

    return false;
  } // consistency check

  public:
    vector<Real> _N;        //require: vector dim = _nodes
    vector<Vector2d > _DN;  //---------- do ------------
    vector<Matrix2d> _DDN;  //---------- do ------------
    Vector2d _coords;
    CornerValences _Valences;
    EigValStruct ** _eigen;

    int _nodes; //Regular patch: 12, Irregular patch: variable

    void _computeFunctions( const SubdivisionMatrix & S, bool needSubdivide);

    void _computeDerivatives( const SubdivisionMatrix & S, bool needSubdivide);

    void _computeSecondDerivatives( const SubdivisionMatrix & S, bool needSubdivide);

    //! compute subdivision Matrix
    void _computeSubdivisionMatrix(SubdivisionMatrix & S, 
				   bool needSubdivide, int N );
    
    void _convertParaCoords();

    void _initialize(const int nodes, const VectorXd & paraCoords);
    
    void _initialize(const int nodes) {
      Vector2d paraCoords(1.0/3.0, 1.0/3.0);
      _initialize(nodes, paraCoords);
    };

    // Methods for evaluation shape functions for arbitrary patch

    void read_eval();
    double EvalBasis(double* c, int I, double v,double w, int K);
    double Eval_D_Basis(double* c, int I, double v,double w, int K, Vector2d& ret);
    void Eval_D2_Basis(double* c, int I, double v,double w, int K, Matrix2d& ret);

  };
} // namespace voom

#endif
