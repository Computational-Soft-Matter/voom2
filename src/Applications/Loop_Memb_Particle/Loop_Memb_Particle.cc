
#include "EigenEllipticResult.h"
#include "LoopShellMesh.h"
#include "GaugeLipid.h"
#include "MembParticleModel.h"
#include "LBFGSB.h"
#include <stdlib.h>
#define PI (3.141592653589793)

using namespace voom;

int main(int argc, char** argv) {
  srand (time(NULL));
  cout << " ------------------------------- " << endl;
  //LoopShellMesh icosa_mesh("sphere_42_nodes.dat","sphere_42_conn.dat");
  //LoopShellMesh icosa_mesh("sphere_212_nodes.dat","sphere_212_conn.dat");
  //LoopShellMesh icosa_mesh("sphere_642_nodes.dat","sphere_642_conn.dat"); // chage quadrature to 3
  LoopShellMesh icosa_mesh("sphere_2562_nodes.dat","sphere_2562_conn.dat");
  //LoopShellMesh icosa_mesh("sphere_10242_nodes.dat","sphere_10242_conn.dat");

  int NumParticles = 200;
  string filename = "N200_kappa2_lambda30";
  //string filename = "RefConfSchematic"; 
  uint NumMat = icosa_mesh.getNumberOfElements();
  uint NodeDoF = 3;
  vector<GaugeLipid *> materials;
  materials.reserve(NumMat);
  Real epsLJ = 0.1;
  Real sigLJ = 0.15;2*sqrt(4.0/NumParticles);//0.15;

  for(int k = 0; k < NumMat; k++)
    materials.push_back(new GaugeLipid(2,0,0,15));

  MembParticleModel model( &icosa_mesh, materials, NodeDoF, NumParticles, sigLJ, epsLJ);

  // Set Requests
  uint PbDoF = (icosa_mesh.getNumberOfNodes())*model.getDoFperNode() + 2*NumParticles;
  int TotNumMatProp = NumMat*2;
  EigenEllipticResult myResults(PbDoF, TotNumMatProp);

  /*
  for (int i=0; i<10000; i++){
    myResults.setRequest(ENERGY);
    model.compute( myResults );
  }
  return 0;

  */
  // model.checkConsistency(myResults, 0, FORCE, 1.0e-5, 1.0e-6); 
  // myResults.setRequest(ENERGY);
  // model.compute( myResults );
  // return 0;
  // cout << setprecision(15) << endl;
  // cout << "Energy : " << myResults.getEnergy()<<endl;
    
    
  // myResults.setRequest(FORCE|ENERGY);
  // model.compute( myResults );
        
  // cout << "Energy : " << myResults.getEnergy()<<endl;
  // for (int i=0; i< PbDoF; i++){
  //   cout << i << " ##### " << (*(myResults._residual))(i) << endl;
  // 	//cout << "Detected at " << i << endl;
  // }
  // return 0;

  // Set bounds
  vector<double > l;
  vector<double > u;
  vector<int > nbd;
  l.assign(PbDoF, 0.0);
  u.assign(PbDoF, 0.0);
  nbd.assign(PbDoF, 0);
  int Tot_len = (icosa_mesh.getNumberOfNodes())*model.getDoFperNode();
  for (int j=0; j< NumParticles; j++) {
    l[Tot_len + 2*j] = 0.00001;
    u[Tot_len + 2*j] = PI-0.00001;
    l[Tot_len + 2*j + 1] = -PI;
    u[Tot_len + 2*j + 1] = PI;
    nbd[Tot_len + 2*j] = 2;
    nbd[Tot_len + 2*j+1] = 2;        
  }
    
  ofstream energyVFile;
  energyVFile.open(filename.c_str());
  energyVFile << "Lambda = " << materials[0]->getLambda() << endl;
  energyVFile << "V " << "E " << "norm " << "\n";
  energyVFile << setprecision(15);
  energyVFile.close();

  for(int k = 0; k < NumMat; k++)
    materials[k]->setRedVol(1.0); //0.98

  LBFGSB mySolver( & model, & myResults, 5, 0, 1.0e-4, 98, 4000);

  mySolver.setBounds( nbd, l, u);
  double kappa = 3;
  for (int j=0;j<=1000;j++){
    energyVFile.open(filename.c_str(), ofstream::app);

    energyVFile << "kappa = " << kappa << ", sigma=" << model._sigma << endl;
    mySolver.solve();
    model.writeOutputVTK(filename,j);

    // cout << " ------------------------------ kappa = " << kappa << endl << endl;
    // if (kappa>0.15){
    //   kappa = kappa-0.1;
    // for(int k = 0; k < NumMat; k++){
    //   materials[k]->setKappa(kappa);
    // }
    // }
    // else if (kappa<=0.15 && kappa>0.06){
    //   kappa = kappa-0.05;
    // for(int k = 0; k < NumMat; k++){
    //   materials[k]->setKappa(kappa);
    // }
    // }
    // else{
    //   kappa = kappa-0.0025;
    //  for(int k = 0; k < NumMat; k++){
    //   materials[k]->setKappa(kappa);
    // }
    // }

    model._sigma = model._sigma+0.01;
    cout << "###########################\n";
    cout << materials[0]->bendingModulus();

    energyVFile.close();
  }

  return 0;
    

}
