//-*-C++-*-
#ifndef __EllipticResult_h__
#define __EllipticResult_h__
#include "Model.h"

namespace voom {
  //! Elliptic Result is a class which contain standard elliptic result and define basic interface to access/update them

  struct EllipticResult 
  {
    // Default constructor
    EllipticResult(): _request(0), _energy(0.0) {};

    // Return NumMat and PbDoF
    virtual int getNumMatProp() = 0;
    virtual int getPbDoF() = 0;

    // Virtual functions to set Elliptic results members (e.g. to be used in the EllipticModel compute function)
    virtual void setRequest(int request) { _request = request; };
    virtual void setEnergy(Real energy) {_energy = energy; };
    virtual void setResidual(int ind, Real value) = 0;

    // Reset function
    virtual void resetResidualToZero() = 0;
    virtual void resetStiffnessToZero() = 0;
    virtual void resetGradgToZero() = 0;
    virtual void resetHgToZero() = 0;
    
    // Virtual functions to add to Elliptic result members
    virtual void addEnergy(Real DeltaEnergy) {_energy += DeltaEnergy; };
    virtual void addResidual(int ind, Real value) = 0;
    virtual void addStiffness(int indRow, int indCol, Real value) = 0;
    virtual void FinalizeGlobalStiffnessAssembly() = 0;
    virtual void setStiffnessFromTriplets(vector<Triplet<Real > > &) = 0;

    virtual void addGradg(int ind, Real value) = 0;
    virtual void addHg(int indRow, int indCol, Real value) = 0;
    virtual void setHgFromTriplets(vector<Triplet<Real > > &) = 0;

    // Virtual functions to get Elliptic result members
    virtual int  getRequest() {return _request; };
    virtual Real getEnergy()  {return _energy; };
    virtual Real getResidual(int ind) = 0;
    virtual Real getStiffness(int indRow, int indCol) = 0;
    virtual Real getGradg(int ind) = 0;
    virtual Real getHg(int indRow, int indCol) = 0;

  private:
    int  _request;
    Real _energy;

  }; // End of declaration of EllipticResult base class

} // namespace voom

#endif
