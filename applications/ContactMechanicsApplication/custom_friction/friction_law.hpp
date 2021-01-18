//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_FRICTION_LAW_H_INCLUDED)
#define KRATOS_FRICTION_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/properties.h"
#include "containers/flags.h"
#include "custom_utilities/contact_properties_extensions.hpp"

namespace Kratos
{
///@addtogroup ContactMechanicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**
 * Base class of friction laws.
 */

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) FrictionLaw
{
public:

  ///@name Type Definitions
  ///@{

  struct FrictionLawVariables {

    double FrictionCoefficient;
    double Alpha;
    double Area;

    double TangentPenalty;

    double PlasticSlipOld;
    double PlasticSlip;
    double Adhesion;

    bool Implex;

    FrictionLawVariables() {Implex = false; Area = 1.0;}

    void Initialize(const double & rTangentPenalty, double PS, double & rArea, bool rImplex = false )
    {
      PlasticSlipOld = PS;
      PlasticSlip = PS;
      Area = rArea;
      TangentPenalty = rTangentPenalty / Area;
      Implex = rImplex;
    };

  };

  /// Pointer definition of FrictionLaw
  KRATOS_CLASS_POINTER_DEFINITION(FrictionLaw);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  FrictionLaw()
  {
    mPlasticSlip = 0.0;
    mPlasticSlipNew = 0.0;
    mDeltaPlasticSlip = 0.0;
  };

  /// Copy constructor.
  FrictionLaw(FrictionLaw const& rOther)
      :mPlasticSlip(rOther.mPlasticSlip)
      ,mPlasticSlipNew(rOther.mPlasticSlipNew)
      ,mDeltaPlasticSlip(rOther.mDeltaPlasticSlip)
  {};

  /// Destructor.
  virtual ~FrictionLaw(){};

  /**
   * Clone function (has to be implemented by any derived class)
   * @return a pointer to a new instance of this friction law
   */
  virtual FrictionLaw::Pointer Clone() const
  {
    return Kratos::make_shared<FrictionLaw>(*this);
  };


  ///@}
  ///@name Operators
  ///@{


  ///@}
  ///@name Operations
  ///@{


  void InitializeSolutionStep();

  void FinalizeSolutionStep();

  // perform similar to a return mapping
  bool EvaluateFrictionLaw( double& rTangentForce, const double& rNormalForce, FrictionLawVariables& rTangentVariables);

  void EvaluateConstitutiveComponents( double& rNormalModulus, double & rTangentModulus, const double& rTangentForce, const double& rEffectiveNormalForce, FrictionLawVariables& rTangentVariables);

  double GetPlasticSlip() { return mPlasticSlip; };

  ///@}
  ///@name Access
  ///@{


  ///@}
  ///@name Inquiry
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  virtual std::string Info() const
  {
    std::stringstream buffer;
    buffer << "FrictionLaw";
    return buffer.str();
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const
  {
    rOStream << "FrictionLaw";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const
  {
    rOStream << "FrictionLaw Data";
  }

  ///@}
  ///@name Friends
  ///@{


  ///@}

protected:
  ///@name Protected static Member Variables
  ///@{


  ///@}
  ///@name Protected member Variables
  ///@{


  ///@}
  ///@name Protected Operators
  ///@{


  ///@}
  ///@name Protected Operations
  ///@{

  // Evaluate Hardening
  virtual double EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) { return 0; };

  // Evaluate Contact Yield Surface
  virtual double EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables ) { return 0; };

  // Evaluate Contact Yield Surface Stress Derivative
  virtual void EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, FrictionLawVariables&  rTangentVariables ) {};

  ///@}
  ///@name Protected  Access
  ///@{


  ///@}
  ///@name Protected Inquiry
  ///@{


  ///@}
  ///@name Protected LifeCycle
  ///@{


  ///@}


private:
  ///@name Static Member Variables
  ///@{


  ///@}
  ///@name Member Variables
  ///@{

  double mPlasticSlip;
  double mPlasticSlipNew;
  double mDeltaPlasticSlip;

  ///@}
  ///@name Private Operators
  ///@{


  ///@}
  ///@name Private Operations
  ///@{

  ///@}
  ///@name Serialization
  ///@{

  friend class Serializer;

  virtual void save( Serializer& rSerializer ) const
  {
    rSerializer.save("mPlasticSlip",mPlasticSlip);
    rSerializer.save("mPlasticSlipNew",mPlasticSlipNew);
    rSerializer.save("mDeltaPlasticSlip",mDeltaPlasticSlip);
  }

  virtual void load( Serializer& rSerializer )
  {
    rSerializer.load("mPlasticSlip",mPlasticSlip);
    rSerializer.load("mPlasticSlipNew",mPlasticSlipNew);
    rSerializer.load("mDeltaPlasticSlip",mDeltaPlasticSlip);
  }

  ///@}
  ///@name Private  Access
  ///@{


  ///@}
  ///@name Private Inquiry
  ///@{


  ///@}
  ///@name Un accessible methods
  ///@{

  /// Assignment operator.
  //FrictionLaw& operator=(FrictionLaw const& rOther);

  ///@}

public:

  DECLARE_HAS_THIS_TYPE_PROPERTIES
  DECLARE_ADD_THIS_TYPE_TO_PROPERTIES
  DECLARE_GET_THIS_TYPE_FROM_PROPERTIES

}; // Class FrictionLaw

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FrictionLaw& rThis)
{
  return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FrictionLaw& rThis)
{
  return rOStream << rThis.Info();
}
///@}

///@} addtogroup block

} // namespace Kratos

#endif // KRATOS_FRICTION_LAW_H_INCLUDED defined
