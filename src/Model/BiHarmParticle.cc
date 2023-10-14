#include "BiHarmParticle.h"
#include "LoopShellElement.h"
namespace voom {
// Constructor
BiHarmParticle::BiHarmParticle(Mesh* aMesh, vector<BiHarmonic * > materials, 
                 const uint NodeDoF, string ghostBCfilename1,
                 string ghostBCfilename2, string boundaryFilename):
  EllipticModel(aMesh, NodeDoF), _materials(materials)
{
  _field.resize(  2*(_myMesh->getNumberOfNodes())  );
  this->initializeField();
  _implicitDynamics = false; //default
  //---- Start:  Boundary Conditions ------
  // Ghost BC1
  ifstream fid_G1; 
  //fid_G1.open("Undap4c2p5_ghostBoundary1Pair.txt", ios::in);
  fid_G1.open(ghostBCfilename1.c_str(), ios::in);
  //fid_G1.open("Unda1c2_ghostBoundary1Pair.txt", ios::in);
  int intNode, ghostNode;
  while(fid_G1 >> intNode >> ghostNode){
    ghostBC1[ghostNode] = intNode;
  }
  fid_G1.close();
  //Ghost BC2
  ifstream fid_G2;
  fid_G2.open(ghostBCfilename2.c_str(), ios::in);
  //fid_G2.open("Undap4c2p5_ghostBoundary2Pair.txt", ios::in);
  //fid_G2.open("Unda1c2_ghostBoundary2Pair.txt", ios::in);
  while(fid_G2 >> intNode >> ghostNode){
    ghostBC2[ghostNode] = intNode;
  }
  fid_G2.close();
  ifstream fid_B; 
  fid_B.open(boundaryFilename.c_str(), ios::in);
  //fid_B.open("Undap4c2p5_BoundaryPair.txt", ios::in);
  //fid_B.open("Unda1c2_BoundaryPair.txt", ios::in);
  int b1Node, b2Node;
  while(fid_B >> b1Node >> b2Node){
    Boundary[b2Node] = b1Node;
  }
  // 
  _NoOfGhostNodes = ghostBC1.size() + ghostBC2.size();
  // Each ghost node contributes two faces/elements
  _NoOfGhostFaces = _NoOfGhostNodes*2;
  _NoOfIrrelevantNodes = _NoOfGhostNodes + Boundary.size();
  fid_B.close();
  cout << "Boundary nodes=" << Boundary.size() << endl;
  cout << "Ghost nodes=" << _NoOfGhostNodes << endl;
  //---- Start:  Boundary Conditions ------
}
// Compute Function - Compute Energy, Force, Stiffness
  void BiHarmParticle::compute(EllipticResult & R)
  {
//     //----- START: This has to be fixed for appropriate BC
//     //Loop over ghostBC maps and enforce periodic BC
//     for (map<int,int>::iterator it=ghostBC1.begin(); it!=ghostBC1.end(); ++it){
//       _field[it->first] = _field[it->second];
//     }
//     for (map<int,int>::iterator it=ghostBC2.begin(); it!=ghostBC2.end(); ++it){
//       _field[it->first] = _field[it->second];
//     }
//     for (map<int,int>::iterator it=Boundary.begin(); it!=Boundary.end(); ++it){
//       _field[it->first] = _field[it->second];
//     }
//     //----- END: This has to be fixed for appropriate BC
//     const vector<GeomElement* > elements = _myMesh->getElements();
//     const int AvgNodePerEl = ((elements[0])->getNodesID()).size();
//     const int NumEl = elements.size() - _NoOfGhostFaces ; //6464; elements.size() - 8*2*2; //56*2*2
//     const int dim = _myMesh->getDimension();
//     int PbDoF = R.getPbDoF();
//     int TotNumMatProp = R.getNumMatProp();

//     // Reset values in result struct
//     if ( R.getRequest() & ENERGY ) {
//       R.setEnergy(0.0);
//     }
//     if ( R.getRequest() & FORCE || R.getRequest() & DMATPROP )  {
//       R.resetResidualToZero();
//     }
//     if ( R.getRequest() & STIFFNESS ) { 
      
//     }
   
//     // Loop through elements, also through material points array, which is unrolled
//     BiHarmonic::Scalarresults ShRes;
//     ShRes.request = R.getRequest();
    
// #ifdef _OPENMP  
// #pragma omp parallel for schedule(static) firstprivate(ShRes)
// #endif

//     for(int e = 0; e < NumEl; e++)
//       {
//         LoopShellElement* element = dynamic_cast<LoopShellElement*>(elements[e]);
//         const vector<int  >& NodesID = element->getNodesID();
//         //cout << "here "<< e << endl;;
//         const int numQP    = element->getNumberOfQuadPoints();
//         const int numNodes = NodesID.size();
//         // Loop over quadrature points
//         for(int q = 0; q < numQP; q++) {
//           // Initialize a and aPartials
//           vector<Vector3d> a(2, Vector3d::Zero());
//           vector<Vector3d> aPartials(3, Vector3d::Zero());
//           Vector3d temp;
//           Vector2d uv;
//           Vector2d uPartials = Vector2d::Zero();
//           Matrix2d u2Partials = Matrix2d::Zero();
// 	  Vector2d vPartials = Vector2d::Zero();
//           Matrix2d v2Partials = Matrix2d::Zero();
// 	  Vector2d uvPrev;
//           for(uint J = 0; J < 2; J++) { //J =0,1 represents two surface coordinates
//             temp = Vector3d::Zero();
//             uv = Vector2d::Zero();
//             uvPrev = Vector2d::Zero();
            
//             for(uint n = 0; n < numNodes; n++) {
//               for(uint i = 0; i < 3; i++) { //i =0,1,2 represents x,y,z
//                 temp(i) += _myMesh->getX(NodesID[n],i) * element->getN(q,n);
//                 a[J](i) += _myMesh->getX(NodesID[n],i) * element->getDN(q, n, J);
//                 aPartials[J](i) += _myMesh->getX(NodesID[n],i) * element->getDDN(q,n,J,J); //a_tt, a_pp
//               }
              
//               uv(0) += _field[2*NodesID[n]] * element->getN(q,n);
// 	      uv(1) += _field[2*NodesID[n]+1] * element->getN(q,n);
//               if (_implicitDynamics == true){
//                 // Implict time stepping ***************************************************** XXXXXXXXXXXXXXXXXXXXXXX
//                 uvPrev(0) += _prevField[2*NodesID[n]] * element->getN(q,n);
// 		uvPrev(1) += _prevField[2*NodesID[n]+1] * element->getN(q,n); 
//               }
//               uPartials(J) += _field[2*NodesID[n]] * element->getDN(q,n,J);
// 	      vPartials(J) += _field[2*NodesID[n]+1] * element->getDN(q,n,J);
//               u2Partials(J,J) += _field[2*NodesID[n]] * element->getDDN(q,n,J,J);
// 	      v2Partials(J,J) += _field[2*NodesID[n]+1] * element->getDDN(q,n,J,J);
//             }
//           }
          
//           for(uint n = 0; n < numNodes; n++) {
//             for(uint i = 0; i < 3; i++) { //i =0,1,2 represents three displacement indices
//               aPartials[2](i) += _myMesh->getX(NodesID[n],i) *element->getDDN(q,n,0,1); //a_tp
//             }
//             u2Partials(0,1) += _field[2*NodesID[n]] * element->getDDN(q,n,0,1);
// 	    v2Partials(0,1) += _field[2*NodesID[n]+1] * element->getDDN(q,n,0,1);
//           }
//           u2Partials(1,0) = u2Partials(0,1);
// 	  v2Partials(1,0) = v2Partials(0,1);
//           ShellGeometry geometry = ShellGeometry(a,aPartials);
//           Real metric = geometry.metric();
//           Matrix2d gab = geometry.metricTensorInverse();
//           Vector2d GammaI = geometry.gklGammaI_kl();

//           _materials[e]->compute(ShRes, uv, uPartials, vPartials, u2Partials, v2Partials, geometry);
//           // Volume associated with QP q
//           Real Vol = element->getQPweights(q) * metric;
          
//           Real dt = 0.01;
//           // Compute energy

//           if (R.getRequest() & ENERGY) {
//             if (_implicitDynamics == true){
//               ShRes.W = ShRes.W + dt/2.0*pow( (phi-phiPrev)/dt, 2.0 );
//             }
//             R.addEnergy(ShRes.W * Vol);
//           }
          
//           // Compute Residual
//           if ( (R.getRequest() & FORCE) ) {     
//             Real tempResidual = 0.0;
//             for(uint n = 0; n < numNodes; n++) {
                
//               tempResidual = ShRes.n1 * element->getN(q,n) +
//                 ShRes.n2(0)*element->getDN(q,n,0) + ShRes.n2(1)*element->getDN(q,n,1) + 
//                 ShRes.n3 * ( gab(0,0) * element->getDDN(q,n,0,0)  + 
//                              gab(0,1) * element->getDDN(q,n,0,1) +
//                              gab(1,0) * element->getDDN(q,n,1,0) +
//                              gab(1,1) * element->getDDN(q,n,1,1) +
//                              - GammaI(0)*element->getDN(q,n,0) 
//                              - GammaI(1)*element->getDN(q,n,1) );
              
              
//               if (_implicitDynamics == true){
//                 tempResidual += (phi-phiPrev)/(dt)*element->getN(q,n);
//               }
//               int mynode = NodesID[n];
//               //Loop over ghostBC maps and enforce periodic BC
//               for (map<int,int>::iterator it=ghostBC1.begin(); it!=ghostBC1.end(); ++it){
//                 if (it->first == NodesID[n])
//                   mynode = it->second;
//               }
//               for (map<int,int>::iterator it=ghostBC2.begin(); it!=ghostBC2.end(); ++it){
//                 if (it->first == NodesID[n])
//                   mynode = it->second;
//               }
//               for (map<int,int>::iterator it=Boundary.begin(); it!=Boundary.end(); ++it){
//                 if (it->first == NodesID[n])
//                   mynode = it->second;
//               }
//               R.addResidual(mynode, tempResidual*Vol);
//               //R.addResidual(mynode, 1.0);
//             } // n loop
//               //cout << endl << e << " " << numNodes << " " << output << endl;        
//           } // Internal force loop
//         }
//       }// Element loop
  } // Compute Mechanics Model
// Writing output
  void BiHarmParticle::writeOutputVTK(const string OutputFile, int step) 
  {
    // Create outputFile name
    stringstream FileNameStream;
    FileNameStream << OutputFile << step << ".vtk";
    ofstream out;
    out.open( (FileNameStream.str()).c_str() );
    int NumNodes = _myMesh->getNumberOfNodes()-_NoOfGhostNodes; //_myMesh->getNumberOfNodes()-8*2; //56*2
    // Header
    char spchar = '#';
    out << spchar << " vtk DataFile Version 3.1" << endl;
    out << "vtk output" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << NumNodes << " FLOAT" << endl;

    // For now we assumed dim == 3 !
    for (int i = 0; i < NumNodes; i++ ) {
      out << (_myMesh->getX(i)).transpose() << endl;
    }

    int noOfElements =  _myMesh->getNumberOfElements()-_NoOfGhostFaces; //6464; _myMesh->getNumberOfElements() - 8*4; //56*4
    out << "POLYGONS " << noOfElements << " " 
        << noOfElements * 4 << endl;
    vector<GeomElement*> elements = _myMesh->getElements();
    for(int e=0; e< noOfElements; e++){
      vector<int> nodesID = elements[e]->getNodesID();
      out << "3 "<<nodesID[0]<< " " << nodesID[1] <<" "<< nodesID[2] << endl;
    }

    out << endl << "POINT_DATA " << NumNodes << endl
        << "SCALARS Field double" << endl
        << "LOOKUP_TABLE default" << endl;
  
    int dim = 1;
    for (int i = 0; i < NumNodes; i++ ) {
      VectorXd x = VectorXd::Zero(dim);
      for (int j = 0; j < dim; j++) {
        x(j) = _field[i*dim + j];
      }
      out << x  << endl;
    }    

    // Close file
    out.close();
  } // writeOutput
} // namespace voom
