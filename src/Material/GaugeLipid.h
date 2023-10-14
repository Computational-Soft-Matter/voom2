#if !defined(__GaugeLipid_h__)
#define __GaugeLipid_h__

#include "ShellMaterial.h"
#include "VoomMath.h"


namespace voom
{

  /*!  Class implementing the simple Helfrich spontaneous curvature
    model for an elastic lipid bilayer.
  */
  
  class GaugeLipid : public ShellMaterial
  {
  protected:

    double _kC;
    double _kG;
		
    double _C0;
    double _H; // mean curvature
    double _K; // gaussian curvature
    double _p;
    double _v;
    double _gamma;
    double _lambda;
    double _epsilon;

  public:

    GaugeLipid(const double kC, const double kG=0.0, const double C0=0.0, const double lambda=1)
      : _kC(kC), _kG(kG), _C0(C0), _lambda(lambda), _epsilon(1e3), _v(1.0) {};
    
    void compute(Shellresults & R, const ShellGeometry & defGeom);

    double meanCurvature() const {return _H;}
    double gaussianCurvature() const {return _K;}

    double spontaneousCurvature() const { return _C0; }
    double bendingModulus() const { return _kC; }
    double gaussianModulus() const { return _kG; }
    
    void setGamma(double gamma) {_gamma = gamma;}
    void setPressure(double p) {_p = p;}
    void setLambda(double lambda) {_lambda=lambda;}
    void setKappa(double k) {_kC = k;}
    double getLambda(){return _lambda;}
    double getGamma(){return _gamma;}
    double getPressure(){return _p;}
    void setEpsilon(double epsilon){_epsilon=epsilon;};
    double getEpsilon(){return _epsilon;};

    void setRedVol(double v){_v=v;};
    double getRedVol(){return _v;};

    void setMaterialParameters(const vector<Real > &) {;}
    void setInternalParameters(const vector<Real > &) {;}

    vector<Real > getMaterialParameters() {;}
    vector<Real > getInternalParameters() {;}

  };
  
} //namespace voom

#endif //  !defined(__GaugeLipid_h__)
