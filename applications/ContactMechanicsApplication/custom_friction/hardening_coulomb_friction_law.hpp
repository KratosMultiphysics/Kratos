//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_HARDENING_COULOMB_FRICTION_LAW_H_INCLUDED)
#define KRATOS_HARDENING_COULOMB_FRICTION_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_friction/coulomb_adhesion_friction_law.hpp"


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

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) HardeningCoulombFrictionLaw
    : public CoulombAdhesionFrictionLaw
{
public:

  ///@name Type Definitions
  ///@{

  /// Pointer definition of HardeningCoulombFrictionLaw
  KRATOS_CLASS_POINTER_DEFINITION( HardeningCoulombFrictionLaw );

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  HardeningCoulombFrictionLaw() : CoulombAdhesionFrictionLaw() {};

  /// Copy constructor.
  HardeningCoulombFrictionLaw(HardeningCoulombFrictionLaw const& rOther) : CoulombAdhesionFrictionLaw(rOther) {};

  /// Destructor.
  ~HardeningCoulombFrictionLaw() override {};


  /**
   * Clone function (has to be implemented by any derived class)
   * @return a pointer to a new instance of this friction law
   */
  FrictionLaw::Pointer Clone() const override
  {
    return Kratos::make_shared<HardeningCoulombFrictionLaw>(*this);
  }

  ///@}
  ///@name Operators
  ///@{


  ///@}
  ///@name Operations
  ///@{


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
  std::string Info() const override
  {
    std::stringstream buffer;
    buffer << "HardeningCoulombFrictionLaw";
    return buffer.str();
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "HardeningCoulombFrictionLaw";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override
  {
    rOStream << "HardeningCoulombFrictionLaw Data";
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

  double EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) override;

  double EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables) override;

  void EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, FrictionLawVariables& rTangentVariables) override;

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


  ///@}
  ///@name Private Operators
  ///@{


  ///@}
  ///@name Private Operations
  ///@{


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
  //HardeningCoulombFrictionLaw& operator=(HardeningCoulombFrictionLaw const& rOther);


  ///@}
  ///@name Serialization
  ///@{

  friend class Serializer;

  void save( Serializer& rSerializer ) const override
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, CoulombAdhesionFrictionLaw )
  }

  void load( Serializer& rSerializer ) override
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, CoulombAdhesionFrictionLaw )
  }

}; // Class HardeningCoulombFrictionLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  HardeningCoulombFrictionLaw& rThis)
{
  return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const HardeningCoulombFrictionLaw& rThis)
{
  return rOStream << rThis.Info();
}
///@}

///@} addtogroup block

} // namespace Kratos

#endif // define KRATOS_HARDENING_COULOMB_FRICTION_LAW_H_INCLUDED
