#include "LoopShellShape.h"
namespace voom{

  
  //! constructor function with parrametric coords
  void LoopShellShape::update( const VectorXd & Point)
  {
    _coords = Point;

    //need subdivision
    bool needSubdivide = false;
    int nds =  _nodes;
    
    // sub-division matrix
    SubdivisionMatrix  S(nds,_nodes);
    S = SubdivisionMatrix::Zero(nds, _nodes);
    _computeSubdivisionMatrix(S, needSubdivide, nds-6);
    
    // compute shape functions
    _computeFunctions(S, needSubdivide);

    // compute the 1st derivatives of shape functions
    _computeDerivatives(S, needSubdivide);

    // compute the 2nd derivatives of shape functions
    _computeSecondDerivatives(S, needSubdivide);
	
  }


  void LoopShellShape::_computeFunctions( const SubdivisionMatrix & SDMatrix, bool needSubdivide )
  {
    // set the parametric coords
    double v = _coords(0);
    double w = _coords(1);
    double u = 1 - v - w;
    
      // Using Jos Stam's superior method for arbitrary (N>=3) valences
      int K=_nodes;
      int N = K-6;

      // determine in which domaon Omega_k^n the parameter lies
      Real n = floor(1-log2(v+w));
      Real pow2 = pow(2,n-1);
      v *= pow2;
      w *= pow2;

      int k;
      if (v> 0.5){
	k=0; v=2*v-1; w=2*w;
      }
      else if (w>0.5){
	k=2; v=2*v; w=2*w-1;
      }
      else{
	k=1; v=1-2*v; w=1-2*w;
      }
    
      for (int i=0; i < K; i++){
	double temp = 0;
	for (int j=0; j < K; j++){
	  temp += _eigen[N-3]->vecI[IX(j,i,K)] *pow(_eigen[N-3]->val[j],n-1) * EvalBasis(_eigen[N-3]->Phi[k], j, v, w, K);
	}
	_N[i] = temp;
      
      }
      SubdivisionMatrix  S(K,K);
      S = SubdivisionMatrix::Zero(K, K);
      _computeSubdivisionMatrix(S, false, N);
      Matrix<double,1,Dynamic> prod(1,K);
      for (int i=0; i<K; i++){
	prod[i] = 0;
	for (int j=0; j<K; j++)
	  prod[i] += S(i,j)*_N[j];
      }
      for (int i=0; i<K; i++)
	_N[i] = prod[i];
		
  }



  void LoopShellShape::_computeDerivatives( const SubdivisionMatrix & SDMatrix, bool needSubdivide )
  {
    // set the parametric coords
    double v = _coords(0);
    double w = _coords(1);
    double u = 1 - v - w;

          // Using Jos Stam's superior method for arbitrary (N>=3) valences
      int K=_nodes;
      int N = K-6;

      // determine in which domaon Omega_k^n the parameter lies
      Real n = floor(1-log2(v+w));
      Real pow2 = pow(2,n-1);
      v *= pow2;
      w *= pow2;
      int k;
      if (v> 0.5){
	k=0; v=2*v-1; w=2*w;
      }
      else if (w>0.5){
	k=2; v=2*v; w=2*w-1;
      }
      else{
	k=1; v=1-2*v; w=1-2*w;
      }
      
      double factr=1.0;
      factr = (k==1 ? -pow(2,n) : pow(2,n)); // because we are taking derivative
      //cout << "k=" << k << ", factr= " << factr << ", n = " << n <<  endl;
      for (int i=0; i < K; i++){
	double temp_0 = 0;
	double temp_1 = 0;
	for (int j=0; j < K; j++){
	  Vector2d ret(Vector2d::Zero());
	  double temp = Eval_D_Basis(_eigen[N-3]->Phi[k], j, v, w, K, ret);
	  double prefact = _eigen[N-3]->vecI[IX(j,i,K)] *pow(_eigen[N-3]->val[j],n-1);
	  temp_0 +=  prefact*ret(0);
	  temp_1 +=  prefact*ret(1);
	}
	_DN[i](0) = factr*temp_0; _DN[i](1) = factr*temp_1;     
      }

      SubdivisionMatrix  S(K,K);
      S = SubdivisionMatrix::Zero(K, K);
      _computeSubdivisionMatrix(S, false, N);
      Matrix<double,2,Dynamic> prod(2,K);
      for (int i=0; i<K; i++){
	prod(0,i) = 0; prod(1,i)=0;
	for (int j=0; j<K; j++){
	  prod(0,i) += S(i,j)*_DN[j](0);
	  prod(1,i) += S(i,j)*_DN[j](1);
	}
	
      }
      for (int i=0; i<K; i++){
	_DN[i](0) = prod(0,i);
	_DN[i](1) = prod(1,i);
      }
    
    return;
		
  }



  void LoopShellShape::_computeSecondDerivatives( const  SubdivisionMatrix & SDMatrix, bool needSubdivide )
  {
    /* second order derivatives of the box spline shape functions */
    /* der(0, *) derivative with respect to vv                     */
    /* der(1, *) derivative with respect to ww                     */
    /* der(2, *) derivative with respect to vw                     */
    double v = _coords(0);
    double w = _coords(1);
    double u = 1 - v - w;

      // Using Jos Stam's superior method for arbitrary (N>=3) valences
      int K=_nodes;
      int N = K-6;
      //cout << "here\n\n" << endl;
      // determine in which domaon Omega_k^n the parameter lies
      Real n = floor(1-log2(v+w));
      Real pow2 = pow(2,n-1);
      v *= pow2;
      w *= pow2;
      int k;
      if (v> 0.5){
	k=0; v=2*v-1; w=2*w;
      }
      else if (w>0.5){
	k=2; v=2*v; w=2*w-1;
      }
      else{
	k=1; v=1-2*v; w=1-2*w;
      }
      double factr=1.0;
      factr = k==1 ? -pow(2,2*n)*factr : pow(2,2*n)*factr; // because we are taking derivative
      for (int i=0; i < K; i++){
	double temp_00 = 0;       double temp_01 = 0; double temp_11= 0;
	for (int j=0; j < K; j++){
	  Matrix2d ret(Matrix2d::Zero());
	  Eval_D2_Basis(_eigen[N-3]->Phi[k], j, v, w, K, ret);
	  double prefact = _eigen[N-3]->vecI[IX(j,i,K)] *pow(_eigen[N-3]->val[j],n-1);
	  temp_00 +=  prefact*ret(0,0);
	  temp_01 +=  prefact * ret(0,1);
	  temp_11 += prefact * ret(1,1);
	
	}
	_DDN[i](0,0) = factr*temp_00; _DDN[i](1,1) = factr*temp_11;
	_DDN[i](1,0) = factr*temp_01; _DDN[i](0,1) = factr*temp_01;
      }
      SubdivisionMatrix  S(K,K);
      S = SubdivisionMatrix::Zero(K, K);
      _computeSubdivisionMatrix(S, false, N);
      Matrix<double,3,Dynamic> prod(3,K);
      for (int i=0; i<K; i++){
	prod(0,i) = 0; prod(1,i)=0; prod(2,i)=0;
	for (int j=0; j<K; j++){
	  prod(0,i) += S(i,j)*_DDN[j](0,0);
	  prod(1,i) += S(i,j)*_DDN[j](1,1);
	  prod(2,i) += S(i,j)*_DDN[j](1,0);
	}
	
      }
      for (int i=0; i<K; i++){
	_DDN[i](0,0) = prod(0,i);
	_DDN[i](1,1) = prod(1,i);
	_DDN[i](0,1) = prod(2,i);
	_DDN[i](1,0) = prod(2,i);
      }
    
    return;	
  }


  void LoopShellShape::_computeSubdivisionMatrix(SubdivisionMatrix & S, 
						 bool needSubdivide, int N)
  {
  	N = _nodes-6;
	//cout << "N=" << _nodes-6 << endl ;
	
	//cout << _Valences[0] << ", " << _Valences[1] << ", " << _Valences[2] << endl;
	//if (_Valences[0]==N) {
	  S(0 ,0 ) = 1.0;
	  S(1 ,1 ) = 1.0;
	  S(2, N) = 1.0;
	  S(3,2) = 1.0;
	  S(4,N+3) = 1.0;
	  S(5,N+2) = 1.0;
	  S(6,N+1) = 1.0;
	  S(7,N+4) = 1.0;
	  S(8,N+5) = 1.0;
	  for (int i=0;i<=N-4; i++)
	    S(9+i,N-i-1) = 1.0;
	  //}
	  /*
	else if (_Valences[2]==N) {
	  S(0 ,N) = 1.0;
	  S(1 ,0 ) = 1.0;
	  S(2, 1) = 1.0;
	  S(N+1,N+3) = 1.0;
	  S(N+2,N+2) = 1.0;
	  S(N+3,N+3) = 1.0;
	  S(N+4,N+4) = 1.0;
	  S(N+5,N+5) = 1.0;
	  for (int i=1;i<=N-2; i++)
	    S(2+i,N-(i-1)-1) = 1.0;
	}
	else if (_Valences[1]==N) {
	  S(0 ,1 ) = 1.0;
	  S(1 ,N ) = 1.0;
	  S(2, 0) = 1.0;
	  S(3,N+1) = 1.0;
	  S(4,N+4) = 1.0;
	  S(5,N+5) = 1.0;
	  S(N+3,2) = 1.0;
	  S(N+4,N+3) = 1.0;
	  S(N+5,N+2) = 1.0;
	  for (int i=0;i<=N-4; i++)
	    S(6+i,N-i-1) = 1.0;
	}
	else {
	  cout << "One of _Valences must be in while computing subdivision matrix!" << endl;
	  exit(-1);

	}
	  */  

    const bool Output_Flag = false;
		
    if( Output_Flag ) {
      //       std::cout.precision(4);
      //       std::cout.setf(std::ios_base::scientific,std::ios_base::floatfield);
      std::cout << "***** Begin Subdivision Matrix *****"<< std::endl;
      //	      << S << std::endl
      for(int i=0; i<S.rows(); i++) {
	std::cout << "    ";
	for(int j=0; j<S.cols(); j++) {
	  std::cout << S(i,j) << "    "; 
	}
	std::cout << "}" << std::endl;
      }
      std::cout << "*****  End Subdivision Matrix  *****"<< std::endl;
    }	
  }



  void LoopShellShape::_convertParaCoords()
  {
    _coords(0) = 1.0 - 2.0 * _coords(0);
    _coords(1) = 1.0 - 2.0 * _coords(1);	
    //! now, the parametric coords are related to the sub-patch
  }



  void LoopShellShape::read_eval ()
  {

    FILE * f;
    int Nmax, N, K;

    if ( !(f=fopen("/Users/sanjay/Google Drive/bucknell/research/projects/FEM_shell/code/voom2/src/Shape/lpdata50NT.dat","rb")) ) {cout<< "Error. lpdata50NT.dat not found" << endl; exit(-1); };
  
    fread ( &Nmax, sizeof(int), 1, f );
    //Nmax = 6;
    _eigen = (EigValStruct **) malloc ( (Nmax-2)*sizeof(EigValStruct *) );
  
    for (int i=0 ; i<Nmax-2 ; i++ )
      {
	N = i+3;
	K = N+6;

	_eigen[i] = (EigValStruct *) malloc ( sizeof(EigValStruct) );
	_eigen[i]->val = (double *) malloc ( K*sizeof(double) );
	_eigen[i]->vecI = (double *) malloc ( K*K*sizeof(double) );
	_eigen[i]->Phi = (double **) malloc ( 3*sizeof(double) );
	_eigen[i]->Phi[0] = (double *) malloc ( K*12*sizeof(double) );
	_eigen[i]->Phi[1] = (double *) malloc ( K*12*sizeof(double) );
	_eigen[i]->Phi[2] = (double *) malloc ( K*12*sizeof(double) );

	fread ( _eigen[i]->val, sizeof(double), K, f );
	fread ( _eigen[i]->vecI, sizeof(double), K*K, f );
	fread ( _eigen[i]->Phi[0], sizeof(double), K*12, f );
	fread ( _eigen[i]->Phi[1], sizeof(double), K*12, f );
	fread ( _eigen[i]->Phi[2], sizeof(double), K*12, f );
      }

    fclose ( f );

    //*pNmax = Nmax;

  }

  double LoopShellShape::EvalBasis(double* c, int I, double v,double w, int K){
    // set the parametric coords
    //const double v = _coords(0);
    //const double w = _coords(1);
    const double u = 1 - v - w;

    Matrix<double,1,12> boxSplines;
    // computing ...
    // 12 shape functions for the regular patch element
    //
    boxSplines(0)=(u*u*u*u + 2.0*u*u*u*v)/12.0;

    boxSplines(1)=(u*u*u*u + 2.0*u*u*u*w)/12.0;

    boxSplines(2)=(u*u*u*u + 2.0*u*u*u*w + 6.0*u*u*u*v + 6.0*u*u*v*w +
		   12.0*u*u*v*v + 6.0*u*v*v*w + 6.0*u*v*v*v + 2.0*v*v*v*w +
		   v*v*v*v)/12.0;

    boxSplines(3)=(6.0*u*u*u*u + 24.0*u*u*u*w + 24.0*u*u*w*w + 8.0*u*w*w*w + 
		   w*w*w*w + 24.0*u*u*u*v + 60.0*u*u*v*w + 36.0*u*v*w*w + 
		   6.0*v*w*w*w + 24.0*u*u*v*v + 36.0*u*v*v*w + 12.0*v*v*w*w + 
		   8.0*u*v*v*v + 6.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(4)=(u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 6.0*u*w*w*w + 
		   w*w*w*w + 2.0*u*u*u*v + 6.0*u*u*v*w + 6.0*u*v*w*w + 
		   2.0*v*w*w*w)/12.0;

    boxSplines(5)=(2.0*u*v*v*v + v*v*v*v)/12.0;

    boxSplines(6)=(u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 
		   6.0*u*w*w*w + w*w*w*w + 8.0*u*u*u*v + 36.0*u*u*v*w + 
		   36.0*u*v*w*w + 8.0*v*w*w*w + 24.0*u*u*v*v + 60.0*u*v*v*w + 
		   24.0*v*v*w*w + 24.0*u*v*v*v + 24.0*v*v*v*w + 
		   6.0*v*v*v*v)/12.0;

    boxSplines(7)=(u*u*u*u + 8.0*u*u*u*w + 24.0*u*u*w*w + 24.0*u*w*w*w + 
		   6.0*w*w*w*w + 6.0*u*u*u*v + 36.0*u*u*v*w + 60.0*u*v*w*w + 
		   24.0*v*w*w*w + 12.0*u*u*v*v + 36.0*u*v*v*w + 
		   24.0*v*v*w*w + 6.0*u*v*v*v + 8.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(8) =(2.0*u*w*w*w + w*w*w*w)/12.0; 
    boxSplines(9)=(2.0*v*v*v*w + v*v*v*v)/12.0;


    boxSplines(10)=(2.0*u*w*w*w + w*w*w*w + 6.0*u*v*w*w + 6.0*v*w*w*w + 
		    6.0*u*v*v*w + 12.0*v*v*w*w + 2.0*u*v*v*v + 
		    6.0*v*v*v*w + v*v*v*v)/12.0;

    boxSplines(11)=(w*w*w*w + 2.0*v*w*w*w)/12.0;



          
    double phi = 0;
    for (int l=0; l<12; l++){
      phi += boxSplines(l)*c[IX(I,l,K)];
      
    }
    return phi;
  }


  double LoopShellShape::Eval_D_Basis(double* c, int I, double v,double w, int K, Vector2d& ret){
    // set the parametric coords
    //const double v = _coords(0);
    //const double w = _coords(1);
    double u = 1 - v - w;

    MatrixXd bsDerivatives(2,12); //rows: directions v,w; cols: shape fns

    // 12 * 2 components of the 1st derivatives of the shape function for the
    // regular element
    bsDerivatives(0,0) = (-6.0*v*u*u - 2.0*u*u*u)/12.0;
    bsDerivatives(1,0) = (-6.0*v*u*u - 4.0*u*u*u)/12.0;

    bsDerivatives(0,1) = (-4.0*u*u*u-6.0*u*u*w)/12.0;
    bsDerivatives(1,1) = (-2.0*u*u*u-6.0*u*u*w)/12.0;

    bsDerivatives(0,2) = (-2.0*v*v*v-6.0*v*v*u
			  + 6.0*v*u*u+2.0*u*u*u)/12.0;
    bsDerivatives(1,2) = (-4.0*v*v*v-18.0*v*v*u
			  - 12.0*v*u*u-2.0*u*u*u
			  - 6.0*v*v*w-12.0*v*u*w
			  - 6.0*u*u*w)/12.0;

    bsDerivatives(0,3) = (-4.0*v*v*v-24.0*v*v*u
			  - 24.0*v*u*u-18.0*v*v*w 
			  - 48.0*v*u*w-12.0*u*u*w
			  - 12.0*v*w*w - 12.0*u*w*w
			  - 2.0*w*w*w)/12.0;

    bsDerivatives(1,3) = (-2.0*v*v*v-12.0*v*v*u
			  - 12.0*v*u*u-12.0*v*v*w
			  - 48.0*v*u*w-24.0*u*u*w
			  - 18.0*v*w*w-24.0*u*w*w
			  - 4.0*w*w*w)/12.0;

    bsDerivatives(0,4) = (-6.0*v*u*u-2.0*u*u*u
			  - 12.0*v*u*w-12.0*u*u*w
			  - 6.0*v*w*w-18.0*u*w*w
			  - 4.0*w*w*w)/12.0;

    bsDerivatives(1,4) = (2.0*u*u*u+6.0*u*u*w
			  - 6.0*u*w*w-2.0*w*w*w)/12.0;

    bsDerivatives(0,5) = (2.0*v*v*v+6.0*v*v*u)/12.0;
    bsDerivatives(1,5) = -v*v*v/6.0;

    bsDerivatives(0,6) = (24.0*v*v*u+24.0*v*u*u
			  + 4.0*u*u*u+12.0*v*v*w
			  + 48.0*v*u*w+18.0*u*u*w
			  + 12.0*v*w*w+12.0*u*w*w
			  + 2.0*w*w*w)/12.0;
  
    bsDerivatives(1,6) = (12.0*v*v*u+12.0*v*u*u
			  + 2.0*u*u*u-12.0*v*v*w
			  + 6.0*u*u*w-12.0*v*w*w
			  - 6.0*u*w*w-2.0*w*w*w)/12.0;

    bsDerivatives(0,7) = (-2.0*v*v*v-6.0*v*v*u
			  + 6.0*v*u*u+2.0*u*u*u
			  - 12.0*v*v*w+12.0*u*u*w
			  - 12.0*v*w*w+12.0*u*w*w)/12.0;

    bsDerivatives(1,7) = (2.0*v*v*v+12.0*v*v*u
			  + 18.0*v*u*u+4.0*u*u*u
			  + 12.0*v*v*w+48.0*v*u*w
			  + 24.0*u*u*w+12.0*v*w*w
			  + 24.0*u*w*w)/12.0;

    bsDerivatives(0,8) = -w*w*w/6.0;
    bsDerivatives(1,8) = (6.0*u*w*w+2.0*w*w*w)/12.0;

    bsDerivatives(0,9) = (4.0*v*v*v+6.0*v*v*w)/12.0;
    bsDerivatives(1,9) = v*v*v/6.0;

    bsDerivatives(0,10) = (2.0*v*v*v+6.0*v*v*u
			   + 12.0*v*v*w+12.0*v*u*w
			   + 18.0*v*w*w+6.0*u*w*w
			   + 4.0*w*w*w)/12.0;

    bsDerivatives(1,10)= (4.0*v*v*v+6.0*v*v*u
			  + 18.0*v*v*w+12.0*v*u*w
			  + 12.0*v*w*w+6.0*u*w*w
			  + 2.0*w*w*w)/12.0;

    bsDerivatives(0,11) = w*w*w/6.0;
    bsDerivatives(1,11) = (6.0*v*w*w+4.0*w*w*w)/12.0;

    //double fact = 1.0; //need to multiply by -2 if subdivision

    double dphi_0 = 0; double dphi_1 = 0;
    for (int l=0; l<12; l++){
      dphi_0 += bsDerivatives(0,l)*c[IX(I,l,K)];
      dphi_1 += bsDerivatives(1,l)*c[IX(I,l,K)];
    }

    
    ret(0) = dphi_0; ret(1) = dphi_1;

    return dphi_0;
  }

  void LoopShellShape::Eval_D2_Basis(double* c, int I, double v,double w, int K, Matrix2d& ret){
    // set the parametric coords
    //const double v = _coords(0);
    //const double w = _coords(1);
    double u = 1 - v - w;

    Matrix<double,3,12> bs2ndDerivatives;

    bs2ndDerivatives(0,0) = v*u;
    bs2ndDerivatives(1,0) = v*u+u*u;
    bs2ndDerivatives(2,0) = (12.0*v*u+6.0*u*u)/12.0;

    bs2ndDerivatives(0,1) = u*u+u*w;
    bs2ndDerivatives(1,1) = u*w;
    bs2ndDerivatives(2,1) = (6.0*u*u+12.0*u*w)/12.0;
             
    bs2ndDerivatives(0,2) = -2.0*v*u;
    bs2ndDerivatives(1,2) = v*v+v*u+v*w+u*w;
    bs2ndDerivatives(2,2) = (6.0*v*v-12.0*v*u
			     -6.0*u*u)/12.0;
             
    bs2ndDerivatives(0,3) = v*v-2.0*u*u
      + v*w-2.0*u*w;
    bs2ndDerivatives(1,3) = -2.0*v*u-2.0*u*u
      + v*w+w*w;
    bs2ndDerivatives(2,3) = (6.0*v*v-12.0*u*u
			     + 24.0*v*w+6.0*w*w)/12.0;
             
    bs2ndDerivatives(0,4) = v*u+v*w+u*w+ w*w;
    bs2ndDerivatives(1,4) = -2.0*u*w;
    bs2ndDerivatives(2,4) = (-6.0*u*u-12.0*u*w 
			     + 6.0*w*w)/12.0;
             
    bs2ndDerivatives(0,5) = v*u;
    bs2ndDerivatives(1,5) = 0.0;
    bs2ndDerivatives(2,5) = -v*v/2.0;
             
    bs2ndDerivatives(0,6) = (-24.0*v*v+12.0*u*u-24.0*v*w
			     + 12.0*u*w)/12.0;
    bs2ndDerivatives(1,6) = (-24.0*v*v-24.0*v*u-24.0*v*w
			     - 24.0*u*w)/12.0;
    bs2ndDerivatives(2,6) = (-12.0*v*v+6.0*u*u-24.0*v*w
			     - 12.0*u*w-6.0*w*w)/12.0;
             
    bs2ndDerivatives(0,7) = -2.0*v*u-2.0*v*w-2.0*u*w- 2.0*w*w;
    bs2ndDerivatives(1,7) = v*u+u*u-2.0*v*w - 2.0*w*w;
    bs2ndDerivatives(2,7) = (-6.0*v*v-12.0*v*u+6.0*u*u 
			     - 24.0*v*w-12.0*w*w)/12.0;
             
    bs2ndDerivatives(0,8) = 0.0;
    bs2ndDerivatives(1,8) = u*w;
    bs2ndDerivatives(2,8) = -w*w/2.0; 
             
    bs2ndDerivatives(0,9) = (12.0*v*v+12.0*v*w)/12.0;
    bs2ndDerivatives(1,9) = 0.0;
    bs2ndDerivatives(2,9) = v*v/2.0;
             
    bs2ndDerivatives(0,10)= (12.0*v*u+12.0*v*w+12.0*u*w
			     + 12.0*w*w)/12.0;
    bs2ndDerivatives(1,10)= v*v+v*u+v*w+u*w;
    bs2ndDerivatives(2,10)= (6.0*v*v+12.0*v*u+24.0*v*w 
			     + 12.0*u*w+6.0*w*w)/12.0;
             
    bs2ndDerivatives(0,11)= 0.0;
    bs2ndDerivatives(1,11)= v*w+w*w;
    bs2ndDerivatives(2,11)= w*w/2.0;

    //double fact = 1.0; //need to multiply by -2 if subdivision

    double d2phi_00 = 0; double d2phi_01 = 0;
    double d2phi_11 = 0;
    for (int l=0; l<12; l++){
      d2phi_00 += bs2ndDerivatives(0,l)*c[IX(I,l,K)];
      d2phi_01 += bs2ndDerivatives(2,l)*c[IX(I,l,K)];
      d2phi_11 += bs2ndDerivatives(1,l)*c[IX(I,l,K)];

    }
    
    ret(0,0) = d2phi_00; ret(0,1) = d2phi_01;
    ret(1,0) = d2phi_01; ret(1,1) = d2phi_11;

    return;
  }
  //double LoopShellShape::

} // namespace voom

/*
  77| 0.0616329, 0.0544738:4       14| 0.116077, 0.0606434:3

0::0.111326||0.106872              0::-0.189859||-0.154509
1::-0.362142||-0.139661            1::0.41051||0.164463
2::-0.133769||-0.360068            2::0.124912||0.391136
3::-0.0715439||0.234421            3::0.0757783||-0.250523
4::-0.0017178||3.90199e-05         4::0.00580701||-0.000260665
5::-0.000181502||-3.90199e-05      5::0.000929878||0.000260665
6::-0.00653194||-0.00655604        6::0.0146479||0.0150665
7::-2.69409e-05||-0.000145327      7::3.71707e-05||0.000287785
8::2.69409e-05||-0.00133837        8::-3.71707e-05||0.00155103
9::0.232861||-0.0678055            9::-0.230352||0.057445
10::0.2317||0.234281               10::-0.212374||-0.224917
--------












*/
