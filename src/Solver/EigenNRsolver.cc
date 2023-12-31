#include "EigenNRsolver.h"

namespace voom
{
  void EigenNRsolver::solve(SolveFor UKN)
  { 
    // Identify number of DoF and unique material models
    uint PbDoF = ( (_myModel->getMesh())->getNumberOfNodes() )*( _myModel->getDoFperNode() );

    // Initialize field vector
    for (int i = 0; i < _DoFid.size(); i++) {
      _myModel->setField(_DoFid[i], _DoFvalues[i]); // Set known displacements
    }

    // Find unique material parameters
    vector<MechanicsMaterial *> materials = _myModel->getMaterials();
    set<MechanicsMaterial *> UNIQUEmaterials;
    for (uint i = 0; i < materials.size(); i++) 
      UNIQUEmaterials.insert(materials[i]);
    uint TotNumMatProp =  UNIQUEmaterials.size()*(materials[0]->getMaterialParameters()).size();

    // Create Eigen Elliptic results
    EigenEllipticResult myResults( PbDoF, TotNumMatProp );

    switch (UKN)
    {
    case 0: // Displacements
      {	
	myResults.setRequest(1); 
	_myModel->compute(myResults);
	cout << "Model energy before solving " << myResults.getEnergy() << endl << endl;
	
	// Initialize solver parameters
	Real error = 1.0;
	uint iter = 0;
	
	// NR loop
	while (iter < _NRmaxIter && error > _NRtol)
	{
	  // Compute stiffness and residual
	  myResults.setRequest(6);
	  _myModel->compute(myResults);
	  
	  // Apply essential boundary conditions
	  this->applyEBC(myResults);
	  //Modify residual where x is known
	  for (int i = 0; i < _DoFid.size(); i++) {
	    myResults.setResidual(_DoFid[i], _DoFvalues[i]); 
	  }
	  // cout << "BC applied" << endl;
	  
	  VectorXd Deltax;
	  // Solve
	  switch (_linSolType)
	  {
	  case 0: 
	    {
	      SimplicialCholesky<SparseMatrix<Real > > chol(*(myResults._stiffness));  // performs a Cholesky factorization of _stiffness
	      Deltax = chol.solve(*(myResults._residual));                             // use the factorization to solve for the given right hand side
	      break;
	    }
	  case 1:
	    {
	      ConjugateGradient<SparseMatrix<Real> > cg;
	      // cg.setTolerance(1.0e-8);
	      // cg.setMaxIterations(100000);
	      cg.compute(*(myResults._stiffness));
	      Deltax; Deltax = cg.solve(*(myResults._residual));
	      // std::cout << "#iterations:     " << cg.iterations() << std::endl;
	      // std::cout << "estimated error: " << cg.error()      << std::endl;
	      break;
	    }
	  case 2:
	    {
	      SparseLU<SparseMatrix<Real, ColMajor> > solver;
	      solver.analyzePattern(*(myResults._stiffness)); 
	      solver.factorize(*(myResults._stiffness));
	      Deltax = solver.solve(*(myResults._residual)); 
	      break;
	    }
	  default: 
	    {
	      cout << "Error - Linear solver type not implemented" << endl;
	    }
	  } // Solve switch 
	  
	  // Change solver solution so that EBC are not added to field in Model multiple times
	  for (int i = 0; i < _DoFid.size(); i++) {
	    Deltax(_DoFid[i]) = 0.0; // So we do not add the known EBC mutiple times
	  }
	  // Update field in the body
	  _myModel->linearizedUpdate(Deltax.data(),-1.0);
	  
	  // Update iter and error
	  iter++;
	  error = Deltax.norm();
	  myResults.setRequest(3);
	  _myModel->compute(myResults);
	  cout << "Energy = " << myResults.getEnergy() << "   - NR iter = " << iter << "   -  NR error = " << error << endl;
	
	} // while loop
	// After finding current field, update prev field
	// _myModel->setPrevField();
	
	break;  
      } // Solve for displacement
    case 1: // Material parameters
      {	
	// Initialize solver parameters
	Real error = 1.0;
	uint iter = 0;
       
	// NR loop
	while (iter < _NRmaxIter && error > _NRtol)
	{
	  // Compute stiffness and residual
	  myResults.setRequest(8);
	  _myModel->compute(myResults);

	  VectorXd DeltaAlpha;
	  // No choice of solver for now - Only cholesky
	  SimplicialCholesky<SparseMatrix<Real > > chol(*(myResults._Hg));  // performs a Cholesky factorization of _stiffness
	  DeltaAlpha = chol.solve(*(myResults._Gradg));                             // use the factorization to solve for the given right hand side

	  // Update material properties in material
	  for (set<MechanicsMaterial *>::iterator MatIt = UNIQUEmaterials.begin(); MatIt != UNIQUEmaterials.end(); MatIt++) {
	    int MatID = (*MatIt)->getMatID();
	    vector<Real > MatProp = (*MatIt)->getMaterialParameters();
	    int NumPropPerMat = MatProp.size();
	    for (int m = 0; m < NumPropPerMat; m++) {
	      MatProp[m] -= DeltaAlpha(MatID*NumPropPerMat + m);
	    }
	    (*MatIt)->setMaterialParameters(MatProp);
	  }
	  
	  // Update iter and error
	  iter++;
	  error = DeltaAlpha.norm();
	  myResults.setRequest(1);
	  _myModel->compute(myResults);
	  cout << "Energy = " << myResults.getEnergy() << "   - NR iter = " << iter << "   -  NR error = " << error << endl;
	
	} // while loop
	break;  
      } // Solve for alphas
    default: 
      {
	cout << "Unknown to solve for not found" << endl;
      }
    } // switch for selecting unknown variables
    
  } // solve function



  void EigenNRsolver::applyEBC(EigenEllipticResult & myResults) {
    
// | A    B |  | x       |   | f |
// |        |  |         | = |   |
// | B^T  C |  | \bar{x} |   | r |

    // Build auxiliary set with EBC dof id
    set<int > EBC;
    for (int i = 0; i < _DoFid.size(); i++) {
      EBC.insert(_DoFid[i]);
    }

    // Fill in matrix B and change matrix A
    for (int k = 0; k < (myResults._stiffness)->outerSize(); k++) {
      for (SparseMatrix<Real>::InnerIterator it(*(myResults._stiffness), k); it; ++it) {
	 
	  // int col = it.col(); == k here because we are using the default column major storing
	  int row = it.row();
	  if ( EBC.find(k) != EBC.end() ) {
	    if ( row == k ) { // setting C = I
	      it.valueRef() = 1; 
	    }
	    else {
	      it.valueRef() = 0;
	    }
	  }
	  else if ( EBC.find(row) != EBC.end() ) { // setting B^T = 0
	    it.valueRef() = 0;
	  }
	    
	} // end of loop over rows
      } // end of loop over columns

    // Clean A
    (myResults._stiffness)->prune(0.0);
      
  } // applyEBC function



} // namespace voom
