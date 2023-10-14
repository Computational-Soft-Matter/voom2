
#ifndef __BiHarmParticle_h__
#define __BiHarmParticle_h__

#include "EllipticModel.h"
#include "BiHarmonic.h"
#include "ShellGeometry.h"
#include "EigenEllipticResult.h"
#include <map>

namespace voom{
  class BiHarmParticle: public EllipticModel {
  public:
    BiHarmParticle(Mesh* aMesh, vector<BiHarmonic * > _materials, 
                   const uint NodeDoF, string ghostBCfilename1,
            string ghostBCfilename2, string boundaryFilename);

~BiHarmParticle() {
  set<BiHarmonic *> UNIQUEmaterials;
  for (uint i = 0; i < _materials.size(); i++) 
    UNIQUEmaterials.insert(_materials[i]);
    
  for (set<BiHarmonic *>::iterator it = UNIQUEmaterials.begin();
       it != UNIQUEmaterials.end(); it++) 
            delete (*it);
};

    void initializeField(const Real value1, const Real value2) {
  const uint numNodes = _myMesh->getNumberOfNodes();
  const uint dim = _myMesh->getDimension();
  for (uint i = 0; i < numNodes; i++) {
    _field[2*i] =  value1;
    _field[2*i+1] =  value2;
  }
};

// Random
void initializeField() {
  const uint numNodes = _myMesh->getNumberOfNodes();
  const uint dim = _myMesh->getDimension();
  for (uint i = 0; i < numNodes; i++) {
    _field[2*i] = (0.5-double(rand())/double(RAND_MAX));
    _field[2*i+1] = (0.5-double(rand())/double(RAND_MAX));   
    //_prevField[i] = 0.0; //*********************************** Implicit Time stepping *****************
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
      _field[id] += value;
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
      //assert(x.size() == _field.size());
      x = _field;
    }
    void setPrevField() {
      _prevField = _field;
    };

    void printField() {
      int i = 0;
      while (i < _field.size()) {
        cout << _field[i] << endl;
        i ++;
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
      set<BiHarmonic *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
        UNIQUEmaterials.insert(_materials[i]);
        
      return UNIQUEmaterials.size();
    }

    uint getTotNumMatProp() {
      // set<BiHarmonic *> UNIQUEmaterials;
      // for (uint i = 0; i < _materials.size(); i++) 
      //        UNIQUEmaterials.insert(_materials[i]);
        
      // return ( UNIQUEmaterials.size() * (_materials[0]->getMaterialParameters()).size() );
    }

    vector<BiHarmonic * > getMaterials() {
      return _materials;
    }
    
    //! Write output
    void writeOutputVTK(const string OutputFile, int step); 
    
    int getNumberOfIrrelevantNodes(){
       return _NoOfIrrelevantNodes; 
    }
    //! Solve the system
    void compute(EllipticResult & R);
    
    void setImplicitDynamicsFlag(bool flag){
      
      _implicitDynamics == flag;
      if(_implicitDynamics == true)
      _prevField.resize( _myMesh->getNumberOfNodes()  ); //********************** Implicity time stepping **************
    
    }
#endif

protected:
    //! Compute Gradients
    void computeGradients(vector<ShellGeometry> defGeomVec, GeomElement* geomEl);

    //! List of Material data at each element in the model
    // (need to be modified for history dependent materials, e.g. plasticity)
    vector<BiHarmonic * > _materials;

    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]

    //vector<Vector3d> _displacements;

    vector<Real > _field;
    vector<Real > _prevField;

    bool _implicitDynamics;
    // Ghost BC Dictionaries
    map<int,int> ghostBC1;
    map<int,int> ghostBC2;
    map<int,int> Boundary;
    int _NoOfGhostFaces;
    int _NoOfGhostNodes;
    int _NoOfIrrelevantNodes;
  };

} // namespace voom
