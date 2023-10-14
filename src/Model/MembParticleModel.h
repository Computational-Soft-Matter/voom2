//-*-C++-*-
#ifndef __MembParticleModel_h__
#define __MembParticleModel_h__

#include "EllipticModel.h"
#include "GaugeLipid.h"
#include "ShellGeometry.h"
#include "EigenEllipticResult.h"
#include "LoopShellElement.h"

namespace voom{

  // Model Results
  class MembParticleModel: public EllipticModel {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs. 
     */
    MembParticleModel(Mesh* aMesh, vector<GaugeLipid * > _materials, 
		      const uint NodeDoF, int NParticles, Real sigmaLJ, Real epsilonLJ=1);
  
    //! Destructor
    ~MembParticleModel() {
      set<GaugeLipid *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      for (set<GaugeLipid *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++) 
	delete (*it);
    };



    //! Initialize field
    // From constant value
    void initializeField(const Real value = 1.0) {
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();
      if (_nodeDoF == dim ){
	//Case: Only displacement DOF
	for (uint i = 0; i < numNodes; i++) 
	  for (uint j = 0; j < dim; j++){
	    double value1 = 1+value*(double(rand())/double(RAND_MAX));
	    _field[i*dim+j] = _myMesh->getX(i,j)*value1; // value = isotropic expansion/shrinking
	  }
	for (int p = 0; p<_NP; p++) {
	  _field[_NMeshDoF+2*p] =  (double(rand())/double(RAND_MAX))*3.14/4 ;//(p+0.01)*3.00/(_NP-1) + value*(double(rand())/double(RAND_MAX))*3.14;
	  _field[_NMeshDoF+2*p+1] = (0.5-double(rand())/double(RAND_MAX))*2*3.14/4;// (6.0-0.02)/(_NP-1)*p+(-3.0+0.01) + value*(0.5-double(rand())/double(RAND_MAX))*2*3.14;

	}
      }
      else if (_nodeDoF == 1){
	// Case: Only scalar field DOF
	for (uint i = 0; i < numNodes; i++) 
	  _field[i] = 0; 
      }
      else if (_nodeDoF == 4){
	// Case: Displacement and scalar DOFs
	for (uint i = 0; i < numNodes; i++){ 
	  for (uint j = 0; j < dim; j++){
	    _field[i*(dim+1)+j] = _myMesh->getX(i,j)*value; // value = isotropic expansion/shrinking      
	  }
	  _field[i*(dim+1)+dim] = 0;
	}
      }
      else{
	cout << "Unknown initalization for _field"<<endl;
	exit(-1);
      }
    };

    //! From array
    void initializeField(const Real* value) {
      _field.assign(value, value+_field.size());
    };



    //! Linearized update
    void linearizedUpdate(const Real* localValues, Real fact) {
      const int nLocalDof = (_myMesh->getNumberOfNodes())*_nodeDoF;
      for(uint i = 0; i < nLocalDof; i++)
	_field[i] += fact*localValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      const uint dim = _myMesh->getDimension();
      // assert( id < _field.size() && dof < dim );
      _field[id*dim + dof] += value;
    }

    void linearizedUpdate(const int id, const int dof, int dim, const Real value) {
      // assert( id < _field.size() && dof < dim );
      _field[id*dim + dof] += value;
    }

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int dof, const Real value) {
      _field[dof] += value;
    }
    
    void setField(uint dof, Real value) {
      _field[dof] = value;
    }
    void setField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    void getField(vector<double> & x) {
      assert(x.size() == _field.size());
      x = _field;
    }

    void setPrevField() {
      _prevField = _field;
    };

    void printField() {
      // Assume dim = 3 - need to be changed
      int i = 0;
      while (i < _field.size()) {
	for (uint j = 0; j < 3; j++) {
	  cout << _field[i] << " ";
	  i++;
	}
	cout << endl;
      }
    }

    void writeField(string FileName) {
      ofstream out;
      out.open( FileName.c_str() );
      out << _field.size() << endl;
      for (uint i = 0; i < _field.size(); i++) {
	out << setprecision(15) << _field[i] << endl;
      }
      out.close();
    }

    uint getNumMat() {
      set<GaugeLipid *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      return UNIQUEmaterials.size();
    }

    uint getTotNumMatProp() {
      set<GaugeLipid *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      return ( UNIQUEmaterials.size() * (_materials[0]->getMaterialParameters()).size() );
    }

    vector<GaugeLipid * > getMaterials() {
      return _materials;
    }
    
    //! Write output
    void writeOutputVTK(const string OutputFile, int step); 

    //! Solve the system
    void compute(EllipticResult & R);

    double particleEnergy(Vector3d& f1, Vector3d& f2);

    Vector3d gradParticleEnergy(Vector3d& f1, Vector3d& f2);
    int computeBarycentricCoords(Real theta, Real phi, Vector3d& a, Vector3d& b, Vector3d& c, Vector2d& vw, Vector2d& vw_theta, Vector2d& vw_phi);
    bool RayIntersectsTriangle(Vector3d& rayVector,
			     Vector3d& v0, Vector3d& v1, Vector3d& v2,
			       Vector3d& outIntersectionPoint);

    /**
     * Given any random vector q, find if it's inside the mesh defined by @var:points and @var:connectivity.
     * param q is the ray vector
     */
    bool find_q_in_tri(Vector3d q, const vector<GeomElement*>& elements, Vector3d &intersect, int& eleID) {

      bool find_point = false;
      eleID = -1;

      // Search element in _hashTable
      int cubex = floor(_ax*q(0)+_bx);
      int cubey = floor(_ay*q(1)+_by);
      int cubez = floor(_az*q(2)+_bz);

      vector <int> cube{cubex, cubey, cubez};
      map<vector<int>, vector<int>>::iterator in;
      in = _hashTable.find(cube);
      vector<int>  subset_elements = in->second;
      // cout << "element:" << endl;
      // cout << cubex <<" " << cubey <<" " << cubez << endl;
      // cout << _ax*q(0)+_bx << "," << _ay*q(1)+_by << ", " << _az*q(2)+_bz << endl;
      // cout << q(0) << "," << q(1) <<"," << q(2) << endl;
      for (int j = 0; j < subset_elements.size(); j++) {

      	//cout<< subset_elements[j] << endl;
      	int i = subset_elements[j];
      	//cout << i << " ";
      	const vector<int  >& con = elements[i]->getNodesID();
      	Vector3d a; a << _myMesh->getX(con[0],0), _myMesh->getX(con[0],1), _myMesh->getX(con[0],2);
      	Vector3d b; b << _myMesh->getX(con[1],0), _myMesh->getX(con[1],1), _myMesh->getX(con[1],2);
      	Vector3d c; c << _myMesh->getX(con[2],0), _myMesh->getX(con[2],1), _myMesh->getX(con[2],2);

        //vector<int> con = elements.at(i);
        if (RayIntersectsTriangle(q, a, b, c, intersect)) {
      	  //cout << "here" << endl;
      	  find_point = true;
      	  eleID = i;
      	  return find_point;
        }
      }
      cout << endl;
     
      return find_point;


      // // Brute force search
      // for (int i = 0; i < elements.size(); i++) {
	
      // 	const vector<int  >& con = elements[i]->getNodesID();
      // 	Vector3d a; a << _myMesh->getX(con[0],0), _myMesh->getX(con[0],1), _myMesh->getX(con[0],2);
      // 	Vector3d b; b << _myMesh->getX(con[1],0), _myMesh->getX(con[1],1), _myMesh->getX(con[1],2);
      // 	Vector3d c; c << _myMesh->getX(con[2],0), _myMesh->getX(con[2],1), _myMesh->getX(con[2],2);
      
      //   //vector<int> con = elements.at(i);
      //   if (RayIntersectsTriangle(q, a, b, c, intersect)) {
      // 	  find_point = true;
      // 	  eleID = i;
      //   }
      // }
      // return find_point;



    }

  public:
    //! List of Material data at each element in the model
    vector<GaugeLipid * > _materials;
    map<vector<int>,vector<int>> _hashTable;
    Real _ax, _ay, _az;
    Real _bx, _by, _bz;

    vector<Real > _field;
    vector<Real > _prevField;

    vector<Vector3d> _particle3DPos; //vector containing particles' 3d positions
    vector< vector<Vector3d> > _aParticle; //vector containing particles' 3d positions
    vector<int> _particlesElementID;
    vector< vector<Real> > _NPart; //stores shape functions values at particle locations. One entry per particle, and entry consists of a vector of shape functions

    
    
    // LJ parameters
    Real _epsilon;
    Real _sigma;
    
    int _NMeshDoF;
    int _NP; // Number of Particles

    int _count;
    
  };

} // namespace voom
#endif
