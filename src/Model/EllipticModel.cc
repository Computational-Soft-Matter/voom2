#include "EllipticModel.h"

namespace voom {

  // Consistency Checks //
  void EllipticModel::checkConsistency(EllipticResult & R, Real perturbationFactor, int request,
				       Real h, Real tol)
  {
    // Check only for local nodes
    //int PARTICLEDOF = 4;

       const uint nodeNum   = _myMesh->getNumberOfNodes();
    const uint nLocalDoF = nodeNum*_nodeDoF + 4;

    // Perturb field randomly to change from reference configuration
    // Save perturbed field to set the configuration back to reference 
    // at the end of the test
    vector<Real > perturb(nLocalDoF, 0.0);
    // srand( time(NULL) );
    // for(uint a = 0; a < nodeNum; a++) {
    //   for(uint i = 0; i < _nodeDoF; i++) {
    // 	Real randomNum =  perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
    // 	perturb[a*_nodeDoF + i] = randomNum;
    // 	this->linearizedUpdate(a, i, randomNum);
    //   }
    // }
    
    


    // Force Check
    if ( request & FORCE ) {
      std::cout << "Checking consistency of forces..." << std::endl;
      Real error = 0.0, norm = 0.0;

      R.setRequest(FORCE); // First compute forces numerically
      this->compute(R);

      R.setRequest(ENERGY); // Reset result request so that only energy is computed 
      // this->compute(R);

      cout << "Model energy at test start = " <<  R.getEnergy() << endl;

      cout << "nodeNum  = " << nodeNum << endl;
      cout << "_nodeDoF = " << _nodeDoF << endl;
      for(int a = 0; a < nodeNum; a++) {		 
	for(int i = 0; i < _nodeDoF; i++) {
	  // if (i%500 == 0 && a%500 ==0)
	  //   cout << '(' << a << ',' << i << ')' << endl;
	  // Perturb +h
	  this->linearizedUpdate(a, i, _nodeDoF, h);
	  this->compute(R);
	  Real Wplus = R.getEnergy();
	  
	  // Perturb -2h
	  this->linearizedUpdate(a, i, _nodeDoF,  -2*h);
	  this->compute(R);
	  Real Wminus = R.getEnergy();
	  
	  // Bring back to original position
	  this->linearizedUpdate(a, i, _nodeDoF, h);

	  //if (a>nodeNum-1)
	    cout << i <<"..." << a<<"..."<< (Wplus-Wminus)/(2.*h) << "<-->" << R.getResidual(a*_nodeDoF + i) << " :::: " << fabs((Wplus-Wminus)/(2.*h) - R.getResidual(a*_nodeDoF + i)) << endl;
	  error += pow( (Wplus-Wminus)/(2.*h) - 
			R.getResidual(a*_nodeDoF + i), 2);
	  norm += pow(R.getResidual(a*_nodeDoF + i), 2);
	} // Loop over dimension
      } // Loop over nodes

      cout << "Particle's Equilibirum" << endl;
      double NP=10;
      double M = nodeNum*3;
      for(int a = 0; a < NP; a++) {		 
	for(int i = 0; i < 2; i++) {
	  // if (i%500 == 0 && a%500 ==0)
	  //   cout << '(' << a << ',' << i << ')' << endl;
	  // Perturb +h
	  this->linearizedUpdate(a, M + i, 2, h);
	  this->compute(R);
	  Real Wplus = R.getEnergy();
	  
	  // Perturb -2h
	  this->linearizedUpdate(a, M+i, 2,  -2*h);
	  this->compute(R);
	  Real Wminus = R.getEnergy();
	  
	  // Bring back to original position
	  this->linearizedUpdate(a, M+i, 2, h);

	  //if (a>nodeNum-1)
	    cout << i <<"..." << a<<"..."<< (Wplus-Wminus)/(2.*h) << "<-->" << R.getResidual(2*a + M + i) << " :::: " << fabs((Wplus-Wminus)/(2.*h) - R.getResidual(2*a + M + i)) << endl;
	  error += pow( (Wplus-Wminus)/(2.*h) - 
			R.getResidual(2*a + M + i), 2);
	  norm += pow(R.getResidual(2*a + M + i), 2);
	} // Loop over dimension
      } // Loop over nodes
      
      error = sqrt(error);
      norm  = sqrt(norm);

      if ( abs(error) < norm * tol) {
	cout << "** Elliptic Model Force consistency check PASSED" << endl;
	cout << "** Error: " << error << " Norm: " << norm  << " Norm*tol: " << norm*tol << endl;
      }
      else {
	cout << "** Elliptic Model Force consistency check FAILED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl;
      }
    } // Check Forces loop

    // Stiffness check
    if ( request & STIFFNESS ) {
      Real error = 0.0, norm = 0.0;

      R.setRequest(4); // First compute stiffness numerically
      this->compute(R);

      R.setRequest(2); // Reset result request so that only forces are computed 
      for(int a = 0; a < nodeNum; a++) {
	for(int i = 0; i < _nodeDoF; i++) {
	  for(int b = 0; b < nodeNum; b++) {
	    for(int j = 0; j < _nodeDoF; j++) {
	      // Perturb +h
	      this->linearizedUpdate(b, j, h);
	      this->compute(R);
	      Real Fplus = R.getResidual(a*_nodeDoF + i);

	      // Perturb -2h
	      this->linearizedUpdate(b, j, -2*h);
	      this->compute(R);
	      Real Fminus = R.getResidual(a*_nodeDoF + i);

	      // Bring back to original position
	      this->linearizedUpdate(b, j, h);

	      // Computing Error and Norm;
	      error += pow((Fplus - Fminus)/(2.*h) - R.getStiffness(a*_nodeDoF+i, b*_nodeDoF+j), 2.0);
	      norm += pow( R.getStiffness(a*_nodeDoF+i, b*_nodeDoF+j), 2.0); 
	    } // j loop
	  } // b loop
	} // i loop
      } // a loop
      
      error = sqrt(error);
      norm  = sqrt(norm);
      if ( abs(error) < norm * tol) {
	cout << "** Elliptic Model Hessian consistency check PASSED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl << endl;
      }
      else {
	cout << "** Elliptic Model Hessian consistency check FAILED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl << endl;
      }
    } // Check Stiffness    

    // Reset field to initial values
    for(uint a = 0; a < nodeNum; a++) {
      for(uint i = 0; i < _nodeDoF; i++) {
	this->linearizedUpdate(a, i, -perturb[a*_nodeDoF + i]);
      }
    } 
    
  } // Check consistency

} // namespace voom




  
