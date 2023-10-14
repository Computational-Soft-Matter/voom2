#include "MembParticleModel.h"
namespace voom {
  // Constructor
  MembParticleModel::MembParticleModel(Mesh* aMesh, vector<GaugeLipid * > materials, 
				       const uint NodeDoF, int NParticles, Real sigmaLJ, Real epsilonLJ):
    EllipticModel(aMesh, NodeDoF), _materials(materials), _NP(NParticles), _epsilon(epsilonLJ), _sigma(sigmaLJ) 
  {

    _NMeshDoF = (_myMesh->getNumberOfNodes() )*_nodeDoF;
    // Resize and initialize (default function) _field vector
    _field.resize( _NMeshDoF  + _NP*2); //2 dof per particle

    _particle3DPos.resize(_NP,Vector3d::Zero());
    _aParticle.resize(_NP);
    _particlesElementID.resize(_NP);
    _NPart.resize(_NP);

    this->initializeField(0.00);

    cout << "_NMeshDOF = " << _NMeshDoF << endl << endl;
    cout << "sigma_LJ = " << _sigma << endl;
    cout << "Number of particles = " << _NP << endl;

    _count = 0;

    ifstream inp("dict.txt");
    int i,j,k;

    string temp;
    cout << "Reading the hash table ...";
    inp >> _ax; inp >> _ay; inp >> _az;
    inp >> _bx; inp >> _by; inp >> _bz;

    cout << _ax << "," << _ay <<"," << _az << endl;
    cout << _bx << "," << _by <<"," << _bz << endl;
    
    while(!inp.eof()){
      getline(inp, temp);
      istringstream iss(temp, istringstream::in);
      int word;
						 
      iss >> i; iss >>j; iss>> k;
      vector<int> ijk{i,j,k};
      vector<int> elements;
      while( iss >> word ){
        elements.push_back(word);
      }
      _hashTable.insert(pair<vector<int>,vector<int>>(ijk, elements));

    }
    cout << "....done" << endl;
    //   while(inp >> i && inp >> j && inp >> k){
    //   vector<int> ijk{i,j,k};
    //   cout << i <<" " <<j << " " << k;
    //   vector<int> elementList;
    //   int temp;
    //   while(      inp >> temp ){
    // 	elementList.push_back(temp);
    // 	//inp>>temp;
    //   }
      
    // }



  }
  
  // Compute Function - Compute Energy, Force, Stiffness
  void MembParticleModel::compute(EllipticResult &Res)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();

    //_particle3DPos.clear();
    //_aParticle.clear();
    //_particlesElementID.clear();
    _NPart.clear();
    
    //_particle3DPos.resize(_NP,Vector3d::Zero());
    //_aParticle.resize(_NP);
    //_particlesElementID.resize(_NP);
    _NPart.resize(_NP);

    
    VectorXd Residual;
    // Reset values in result struct
    if (Res.getRequest() & ENERGY ) {
      Res.setEnergy(0.0);
    }
    if (Res.getRequest() & FORCE ||Res.getRequest() & DMATPROP )  {
      Res.resetResidualToZero();
    }
    if (Res.getRequest() & STIFFNESS ) { 
      Res.resetStiffnessToZero();
    }
    // Loop through elements to compute area
    double epsilon = _materials[0]->getEpsilon(); //area constraint penalty
    double area = 0;
    double RefArea = 0;
    double volume = 0;
    double gaugeConstX = 0; double gaugeConstY = 0; double gaugeConstZ = 0;
    double Vbar = _materials[0]->getRedVol()*4.18879; //volume of a unit sphere = 4.18879
    double lambda = _materials[0]->getLambda();
    double press = 0;
    double epsilonA = epsilon;
    double epsilonX = epsilon; double epsilonY = epsilon; double epsilonZ = epsilon;
    double epsilonV = epsilon*0; //No Vol constr 
    Real trX = 0.0; Real trY = 0.0; Real trZ = 0.0;
    
    double gaugeEnergy = 0.0;

    
    // Compute Particles' energies and forces
    for (int k=0; k < _NP; k++){
      // 1. Locate element corresponding to each particle
      _particle3DPos[k] = Vector3d::Zero();
      double theta_k = _field[_NMeshDoF + 2*k];   // global coordinate1 of the particle
      double phi_k   = _field[_NMeshDoF + 2*k+1]; // global coordinate2 of the particle

      // 2. Use (theta_i,phi_i) to locate element ID and compute Barycentric coordinates in the element
      Vector3d q; q<< sin(theta_k) * cos(phi_k), sin(theta_k) * sin(phi_k), cos(theta_k);
      Vector3d intersect;
      int eleID;
      find_q_in_tri(q, elements, intersect, eleID);
      _particlesElementID[k] = eleID;
      LoopShellElement* ele = dynamic_cast<LoopShellElement*>(elements[eleID]);
      int numNodes = ele->getNodesPerElement();
      const vector<int  >& NodesID = ele->getNodesID();
      Vector3d a; a << _myMesh->getX(NodesID[0],0), _myMesh->getX(NodesID[0],1), _myMesh->getX(NodesID[0],2);
      Vector3d b; b << _myMesh->getX(NodesID[1],0), _myMesh->getX(NodesID[1],1), _myMesh->getX(NodesID[1],2);
      Vector3d c; c << _myMesh->getX(NodesID[2],0), _myMesh->getX(NodesID[2],1), _myMesh->getX(NodesID[2],2);
      
      Vector2d vw = Vector2d::Zero();
      Vector2d vw_theta = Vector2d::Zero();
      Vector2d vw_phi = Vector2d::Zero();
      computeBarycentricCoords(theta_k,phi_k, a, b, c, vw, vw_theta, vw_phi);

      // 3. Use element ID and barycentric coordinates to compute 3d position, derivatives, etc
      //LoopShellShape* shapeAtP = new LoopShellShape(numNodes, ele->_cornerValences, vw); //Shape function and derivatives at particle position
      LoopShellShape shapeAtP(numNodes, ele->_cornerValences, vw); //Shape function and derivatives at particle position
      
      vector<Vector3d> a_particle(2, Vector3d::Zero());
      Vector3d R_particle;

      vector<Real> NAtP;
      
      for(uint J = 0; J < 2; J++) { //J =0,1 (represents two surface coordinates
	R_particle = Vector3d::Zero();	    
	for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	  for(uint n = 0; n < numNodes; n++) {
	    R_particle(i) += _field[NodesID[n]*dim + i] * shapeAtP.getN(n);
	    
	    if (J==0)
	      a_particle[J](i) += _field[NodesID[n]*dim + i] * (shapeAtP.getDN(n, 0)*vw_theta(0) + shapeAtP.getDN(n, 1)*vw_theta(1));
	    if (J==1)
	      a_particle[J](i) += _field[NodesID[n]*dim + i] * (shapeAtP.getDN(n, 0)*vw_phi(0) + shapeAtP.getDN(n, 1)*vw_phi(1));
	    
	    // Store shape functions and derivatives (to be used for residue computation)
	    if (i==0 && J==0){
	      double temp =  shapeAtP.getN(n);
	      NAtP.push_back(temp);
	    }
	  }
	  // Store shape functions and derivatives (to be used for residue computation)
	  if (i==0 && J==0) {
	    //_NPart.push_back(NAtP);
	    _NPart[k] = NAtP;
	    //copy(NAtP.begin(), NAtP.end(), back_inserter(_NPart[k]));
	  }
	}
      }
      //cout << "in model " << k << endl;
      // 4. Store these in _particle3DPos, etc
      _particle3DPos[k] << R_particle(0) , R_particle(1) , R_particle(2);
      _aParticle[k].resize(2,Vector3d::Zero());
      _aParticle[k][0] << a_particle[0](0) , a_particle[0](1) , a_particle[0](2);
      _aParticle[k][1] << a_particle[1](0) , a_particle[1](1) , a_particle[1](2);
      //delete shapeAtP;
      }

    // 5. Compute total interaction energy of all particles
    for (int i=0; i <_NP; i++){
      LoopShellElement* ele_i = dynamic_cast<LoopShellElement*>(elements[_particlesElementID[i]]);
      const vector<int  >& NodesID_i = ele_i->getNodesID();
      int numNodes_i = ele_i->getNodesPerElement();
      for (int j=i+1; j<_NP; j++){
	LoopShellElement* ele_j = dynamic_cast<LoopShellElement*>(elements[_particlesElementID[j]]);
	const vector<int  >& NodesID_j = ele_j->getNodesID();
	int numNodes_j = ele_j->getNodesPerElement();
	double energy_ij = particleEnergy(_particle3DPos[i],_particle3DPos[j]);
	
	if (Res.getRequest() & ENERGY){
	  Res.addEnergy(energy_ij);
	}

	if (Res.getRequest() & FORCE) {
	  Vector3d gij = gradParticleEnergy(_particle3DPos[i],_particle3DPos[j]);
	  for (int dd=0;dd<dim;dd++){
	    for (int n=0; n<numNodes_i;n++){
	      Res.addResidual(NodesID_i[n]*dim+dd, gij(dd)*(_NPart[i][n]));
	    }
	  }

	  for (int dd=0;dd<dim;dd++){
	    for (int n=0; n<numNodes_j;n++){
	      Res.addResidual(NodesID_j[n]*dim+dd, -gij(dd)* _NPart[j][n]);
	    }
	  }

	  Res.addResidual(_NMeshDoF+2*i,   gij(0)*_aParticle[i][0](0)+gij(1)*_aParticle[i][0](1)+gij(2)*_aParticle[i][0](2)); //gij.a_theta
	  Res.addResidual(_NMeshDoF+2*i+1, gij(0)*_aParticle[i][1](0)+gij(1)*_aParticle[i][1](1)+gij(2)*_aParticle[i][1](2)); //gij.a_phi

	  Res.addResidual(_NMeshDoF+2*j,   -(gij(0)*_aParticle[j][0](0)+gij(1)*_aParticle[j][0](1)+gij(2)*_aParticle[j][0](2)) ); 
	  Res.addResidual(_NMeshDoF+2*j+1, -(gij(0)*_aParticle[j][1](0)+gij(1)*_aParticle[j][1](1)+gij(2)*_aParticle[j][1](2))) ;
	}
	
      }
    }
    GaugeLipid::Shellresults ShRes;
    ShRes.request =Res.getRequest();

    for(int e = 0; e < NumEl; e++)
      {
	LoopShellElement* element = dynamic_cast<LoopShellElement*>(elements[e]);
	const vector<int  >& NodesID = element->getNodesID();
	const int numQP    = element->getNumberOfQuadPoints();
	const int numNodes = NodesID.size();
      
	// Loop over quadrature points
	for(int q = 0; q < numQP; q++) {
	  // Initialize a and aPartials
	  vector<Vector3d> a(2, Vector3d::Zero());
	  vector<Vector3d> aPartials(3, Vector3d::Zero());
	  Vector3d R;
	  vector<Vector3d> Refa(2, Vector3d::Zero());
	  vector<Vector3d> RefaPartials(3, Vector3d::Zero());
	  Vector3d RefR;
	  
	  for(uint J = 0; J < 2; J++) { //J =0,1 (represents two surface coordinates
	    R = Vector3d::Zero();
	    RefR = Vector3d::Zero();
	    for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	      for(uint n = 0; n < numNodes; n++) {
		R(i) += _field[NodesID[n]*dim + i] * element->getN(q,n);
		a[J](i) += _field[NodesID[n]*dim + i] * element->getDN(q, n, J);
		aPartials[J](i) += _field[NodesID[n]*dim + i] * element->getDDN(q,n,J,J); //a_tt, a_pp

		RefR(i) += _myMesh->getX(NodesID[n],i) * element->getN(q,n);
		Refa[J](i) += _myMesh->getX(NodesID[n],i) * element->getDN(q, n, J);
		RefaPartials[J](i) += _myMesh->getX(NodesID[n],i) * element->getDDN(q,n,J,J); //a_tt, a_pp
	      }
	    }
	  }
	  // a_tp component
	  for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	    for(uint n = 0; n < numNodes; n++) {
	      aPartials[2](i) += _field[NodesID[n]*dim + i]*element->getDDN(q,n,0,1); //a_tp
	      RefaPartials[2](i) += _myMesh->getX(NodesID[n],i) *element->getDDN(q,n,0,1); //a_tp
	    }
	  }
	  ShellGeometry geometry = ShellGeometry(a,aPartials);
	  ShellGeometry Refgeometry = ShellGeometry(Refa,RefaPartials);
	  _materials[e]->setRefGeometry(Refgeometry);
	  Real metric = geometry.metric();
	  // Volume associated with QP q
	  Real da = element->getQPweights(q) * metric;
	  Real dA = element->getQPweights(q) * Refgeometry.metric();
	  area += da;
	  RefArea += dA;
	  volume += da*(geometry.d()).dot(R)/3.;
	  gaugeConstX += (RefR(0)-trX)*da;
	  gaugeConstY += (RefR(1)-trY)*da;
	  gaugeConstZ += (RefR(2)-trZ)*da;

	  //---------------pre compute gauge energy
	  _materials[e]->compute(ShRes, geometry);
	 
	  gaugeEnergy += ShRes.Wg*da;
	  // ------------ pre compute gauge energy
	}
      }
    // End: Total area computation
    // Loop through elements, also through material points array, which is unrolled
    for(int e = 0; e < NumEl; e++)
      {
	LoopShellElement* element = dynamic_cast<LoopShellElement*>(elements[e]);
	const vector<int  >& NodesID = element->getNodesID();
	const int numQP    = element->getNumberOfQuadPoints();
	const int numNodes = NodesID.size();
      
	// Loop over quadrature points
	for(int q = 0; q < numQP; q++) {
	  // Initialize a and aPartials
	  vector<Vector3d> a(2, Vector3d::Zero());
	  vector<Vector3d> aPartials(3, Vector3d::Zero());
	  Vector3d R;
	  vector<Vector3d> Refa(2, Vector3d::Zero());
	  vector<Vector3d> RefaPartials(3, Vector3d::Zero());
	  Vector3d RefR;
	  
	  for(uint J = 0; J < 2; J++) { //J =0,1 (represents two surface coordinates
	    R = Vector3d::Zero();
	    RefR = Vector3d::Zero();
	    for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	      for(uint n = 0; n < numNodes; n++) {
		R(i) += _field[NodesID[n]*dim + i] * element->getN(q,n);
		a[J](i) += _field[NodesID[n]*dim + i] * element->getDN(q, n, J);
		aPartials[J](i) += _field[NodesID[n]*dim + i] * element->getDDN(q,n,J,J); //a_tt, a_pp

		RefR(i) += _myMesh->getX(NodesID[n],i) * element->getN(q,n);
		Refa[J](i) += _myMesh->getX(NodesID[n],i) * element->getDN(q, n, J);
		RefaPartials[J](i) += _myMesh->getX(NodesID[n],i) * element->getDDN(q,n,J,J); //a_tt, a_pp
	      }
	    }
	    
	  }
	  // a_tp component
	  for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	    for(uint n = 0; n < numNodes; n++) {
	      aPartials[2](i) += _field[NodesID[n]*dim + i]*element->getDDN(q,n,0,1); //a_tp
	      RefaPartials[2](i) += _myMesh->getX(NodesID[n],i) *element->getDDN(q,n,0,1); //a_tp
	    }
	  }
	  ShellGeometry geometry = ShellGeometry(a,aPartials);
	  ShellGeometry Refgeometry = ShellGeometry(Refa,RefaPartials);
	  _materials[e]->setRefGeometry(Refgeometry);
	  Real metric = geometry.metric();
	  _materials[e]->compute(ShRes, geometry);
	  Matrix2d gab = geometry.metricTensorInverse();
	  // Volume associated with QP q
	  Real Vol = element->getQPweights(q) * metric;
	  Real RefVol = element->getQPweights(q) * Refgeometry.metric();
	  double gamma = (area-RefArea)*epsilonA;
	  //TOCLEAN
	  //gaugeConstX=0; gaugeConstY=0;gaugeConstZ=0;

	  double gammaGaugeX = gaugeConstX*epsilonX;
	  double gammaGaugeY = gaugeConstY*epsilonY;
	  double gammaGaugeZ = gaugeConstZ*epsilonZ;
	  press = - (volume-Vbar)*epsilonV; // Pressure as a penulty

	  // Compute energy
	  if (Res.getRequest() & ENERGY) {
	    Res.addEnergy((ShRes.W)*Vol);
	   
	  }
	  // Compute Residual
	  if ( (Res.getRequest() & FORCE) ) {
	    for(uint n = 0; n < numNodes; n++) {
	      for(uint i = 0; i < dim; i++) {
		Real tempResidual = 0.0;
		for ( int alpha = 0; alpha < 2; alpha++){
		  // Stress Resultant part
		  tempResidual += (ShRes.n[alpha](i) + 2*gaugeEnergy*ShRes.ng[alpha](i)
				   + (gamma + (RefR(0)-trX)*gammaGaugeX + (RefR(1)-trY)*gammaGaugeY
				      + (RefR(2)-trZ)*gammaGaugeZ)*(geometry.aDual())[alpha](i)) *  element->getDN(q,n,alpha);
		  // Moment Resultant part
		  for(int beta=0; beta<2; beta++) {
		    uint index = -1;
		    if (alpha==0 && beta == 0) index = 0;
		    else if (alpha==1 && beta ==1) index = 1;
		    else if (alpha==1 && beta ==0) index =2;
		    else if (alpha==0 && beta ==1) index =3;
		    else exit(-1);
		    tempResidual +=  -  ShRes.m[alpha].dot(geometry.aDual()[beta]) *
		      ( element->getDDN(q,n,alpha,beta)*geometry.d()(i) + element->getDN(q,n,beta)*geometry.dPartials()[alpha](i) )
		      -  ShRes.m[alpha].dot(geometry.aDualPartials()[index])*element->getDN(q,n,beta)*geometry.d()(i);
		  }//beta loop
		  // ------------------ volume variation the alpha contribution
		  tempResidual += -press/3.0* (R.dot(geometry.d())*geometry.aDual()[alpha] - R.dot(geometry.aDual()[alpha])*geometry.d())(i)*element->getDN(q,n,alpha);
		}// alpha loop
		// Volume variation
		// ------------------ volume variation 
		tempResidual += -press/3.0*(geometry.d())(i) * element->getN(q,n);
		Res.addResidual(NodesID[n]*dim+i, tempResidual*Vol);
	      } //i loop
	    } // n loop
	  } // Internal force loop
	    // cout << setprecision(9) ;
	    // cout << gamma << " " << press << endl;
	}
      }// Element loop
    // Add global energy integral at the end
    if (Res.getRequest() & ENERGY) {
      Res.addEnergy(pow((area-RefArea),2)/(2)*epsilonA);
      Res.addEnergy(pow((volume-Vbar),2)/(2)*epsilonV); 
      Res.addEnergy(pow((gaugeConstX),2)/(2)*epsilonX);
      Res.addEnergy(pow((gaugeConstY),2)/(2)*epsilonY);
      Res.addEnergy(pow((gaugeConstZ),2)/(2)*epsilonZ);
      Res.addEnergy(pow(gaugeEnergy,2));
    }
  } // End: Compute Model  

  // Writing output
  void MembParticleModel::writeOutputVTK(const string OutputFile, int step) 
  {
    // Create outputFile name (surface)
    stringstream FileNameStream;
    FileNameStream << OutputFile << "_" << step << ".vtk";
    ofstream out;
    out.open( (FileNameStream.str()).c_str() );
    int NumNodes = _myMesh->getNumberOfNodes(); //_myMesh->getNumberOfNodes()-8*2; //56*2
    // Header
    char spchar = '#';
    out << spchar << " vtk DataFile Version 3.1" << endl;
    out << "vtk output" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << NumNodes << " FLOAT" << endl;

    int dim = 3;
    for (int i = 0; i < NumNodes; i++ ) {
      VectorXd x = VectorXd::Zero(dim);
      for (int j = 0; j < dim; j++) {
        x(j) = _field[i*dim + j];
      }
      out << (x).transpose()  << endl;
    }    

    int noOfElements =  _myMesh->getNumberOfElements(); 
    out << "POLYGONS " << noOfElements << " " 
        << noOfElements * 4 << endl;
    vector<GeomElement*> elements = _myMesh->getElements();
    for(int e=0; e< noOfElements; e++){
      vector<int> nodesID = elements[e]->getNodesID();
      out << "3 "<<nodesID[0]<< " " << nodesID[1] <<" "<< nodesID[2] << endl;
    }

    // Close file
    out.close();


    // Create outputFile (Particle)
    stringstream FileNameStreamPt;
    FileNameStreamPt << OutputFile << "_" << step << "_Ps.vtk";
    ofstream outPt;
    outPt.open( (FileNameStreamPt.str()).c_str() );
    //int NumNodes = _myMesh->getNumberOfNodes(); //_myMesh->getNumberOfNodes()-8*2; //56*2
    // Header
    char spcharPt = '#';
    outPt << spcharPt << " vtk DataFile Version 3.1" << endl;
    outPt << "vtk output" << endl;
    outPt << "ASCII" << endl;
    outPt << "DATASET POLYDATA" << endl;
    outPt << "POINTS " << _NP << " FLOAT" << endl;

    for (int k = 0; k < _NP; k++ ) {
      outPt << _particle3DPos[k](0) << " " << _particle3DPos[k](1) << " " << _particle3DPos[k](2) << endl;
    }    

    // Close file
    outPt.close();
    
  } // writeOutput

  double MembParticleModel::particleEnergy(Vector3d& fi, Vector3d& fj){
    
    
    double r = (fi-fj).norm();//sqrt(pow(fi(0)-fj(0),2)+pow(fi(1)-fj(1),2)+pow(fi(2)-fj(2),2));
    double sigma_r6 = pow(_sigma/r,6);
      
    // Lennard-Jones Potential
    return ( 4*_epsilon*( pow(sigma_r6,2) - sigma_r6 ) );
    // Harmonic
    //return ( _epsilon*pow(r-_sigma, 2) );

    // Electrostatic

    //return ( 1/r );
  }

  Vector3d MembParticleModel::gradParticleEnergy(Vector3d& fi, Vector3d& fj){

    double r = (fi-fj).norm();//sqrt(pow(fi(0)-fj(0),2)+pow(fi(1)-fj(1),2)+pow(fi(2)-fj(2),2));
    double sigma_r6 = pow(_sigma/r,6);
 
    // Lennard-Jones Potential
    return (-4*_epsilon*( 12*pow(_sigma/r,11) - 6*pow(_sigma/r,5) )*(_sigma/(r*r))/r )*(fi-fj);
    
    //return ((24 *_epsilon/r * sigma_r6 * (1-2*sigma_r6)/r)*(fi-fj)); //
    //Harmonic
    //return (2*_epsilon*(r-_sigma)*(fi-fj)/r);

    //Electrostatic
    //return (-1/pow(r,3)*(fi-fj));
    
  }

  int MembParticleModel::computeBarycentricCoords(Real theta, Real phi, Vector3d& a, Vector3d& b, Vector3d& c, Vector2d& vw, Vector2d& vw_theta, Vector2d& vw_phi){
    //int MembParticleModel::computeBarycentricCoords(Real theta, Real phi, Vector3d& a, Vector3d& b, Vector3d& c, Vector2d& vw, Vector2d& vw_theta, Vector2d& vw_phi, Real& p_ret, Real& p_theta_ret){

    Real cp = cos(phi); Real sp = sin(phi);
    Real ct = cos(theta); Real st = sin(theta);
    
    Vector3d er; er << cp*st, sp*st, ct;

    Vector3d er_theta;
    er_theta << cp*ct, sp*ct, -st;

    Vector3d er_phi;
    er_phi << -sp*st, cp*st, 0;
    
    Vector3d v0 = b - a, v1 = c - a;
    Vector3d n = v0.cross(v1);
    double erdotn = (er.dot(n));
    Real adotn = a.dot(n);
    Real t = adotn/erdotn;

    Real t_theta  = -adotn/(erdotn*erdotn)*(er_theta.dot(n));
    Real t_phi  = -adotn/(erdotn*erdotn)*(er_phi.dot(n));
    
    Vector3d p = t*er;
    Vector3d p_theta = t_theta*er + t*er_theta;
    Vector3d p_phi = t_phi*er + t*er_phi;
    
    
    Vector3d v2 = p - a;
    Real d00 = v0.dot(v0); //Dot(v0, v0);
    Real d01 = v0.dot(v1); //Dot(v0, v1);
    Real d11 = v1.dot(v1); //Dot(v1, v1);
    Real d20 = v2.dot(v0); //Dot(v2, v0);
    Real d21 = v2.dot(v1); //Dot(v2, v1);
    Real denom = d00 * d11 - d01 * d01;
    vw(0) = (d11 * d20 - d01 * d21) / denom;
    vw(1) = (d00 * d21 - d01 * d20) / denom;

    Real d20_theta = p_theta.dot(v0); Real d20_phi = p_phi.dot(v0);
    Real d21_theta = p_theta.dot(v1); Real d21_phi = p_phi.dot(v1);

    //cout << "denom = " << d11*d20_theta - d01*d21_theta << " " << denom << endl;

    Real v_theta = (d11*d20_theta - d01*d21_theta)/denom;
    Real v_phi = (d11*d20_phi - d01*d21_phi)/denom;

    Real w_theta = (d00*d21_theta - d01*d20_theta)/denom;
    Real w_phi = (d00*d21_phi - d01*d20_phi)/denom;

    vw_theta(0) = v_theta; vw_theta(1) =  w_theta;
    vw_phi(0) = v_phi; vw_phi(1) = w_phi;

    //p_ret = d20; p_theta_ret = d20_phi;
    
    return 0;
    //u = 1.0f - v - w;
    
  }

  bool MembParticleModel::RayIntersectsTriangle(Vector3d& rayVector,
						Vector3d& v0, Vector3d& v1, Vector3d& v2,
						Vector3d& outIntersectionPoint)
  {
    const Real EPSILON = 0.0000001;
    Vector3d edge1, edge2, h, s, q;
    Real a,f,u,v;
    edge1 = (v1 - v0);
    edge2 = (v2 - v0);
    h = rayVector.cross(edge2);
    a = edge1.dot(h);
    if (a > -EPSILON && a < EPSILON)
      return false;    // This ray is parallel to this triangle.
    f = 1.0/a;

    s =  -v0;
    u = f * s.dot(h);
    if (u < 0.0 || u > 1.0)
      return false;
    q = s.cross(edge1);
    v = f * rayVector.dot(q);
    if (v < 0.0 || u + v > 1.0)
      return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    Real t = f * edge2.dot(q);
    if (t > EPSILON) // ray intersection
      {
        outIntersectionPoint = t*rayVector;
        return true;
      }
    else // This means that there is a line intersection but not a ray intersection.
      return false;
  }

} // namespace voom
