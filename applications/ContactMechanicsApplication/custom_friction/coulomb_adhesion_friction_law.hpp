//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_COULOMB_ADHESION_FRICTION_LAW_H_INCLUDED)
#define      KRATOS_COULOMB_ADHESION_FRICTION_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_friction/friction_law.hpp"


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

  class KRATOS_API(CONTACT_MECHANICS_APPLICATION) CoulombAdhesionFrictionLaw
    : public FrictionLaw
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of CoulombAdhesionFrictionLaw
    KRATOS_CLASS_POINTER_DEFINITION( CoulombAdhesionFrictionLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CoulombAdhesionFrictionLaw();

    /// Destructor.
    virtual ~CoulombAdhesionFrictionLaw();


    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this friction law
     */
    FrictionLaw::Pointer Clone() const override;

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
    //std::string Info() const override;

    /// Print information about this object.
    //void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    //void PrintData(std::ostream& rOStream) const override;


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
    //CoulombAdhesionFrictionLaw& operator=(CoulombAdhesionFrictionLaw const& rOther);

    /// Copy constructor.
    //CoulombAdhesionFrictionLaw(CoulombAdhesionFrictionLaw const& rOther);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FrictionLaw )
    }

    void load( Serializer& rSerializer ) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FrictionLaw )
    }

  }; // Class CoulombAdhesionFrictionLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  // inline std::istream& operator >> (std::istream& rIStream,
  // 				    CoulombAdhesionFrictionLaw& rThis);

  /// output stream function
  // inline std::ostream& operator << (std::ostream& rOStream,
  // 				    const CoulombAdhesionFrictionLaw& rThis)
  //   {
  //     rThis.PrintInfo(rOStream);
  //     rOStream << std::endl;
  //     rThis.PrintData(rOStream);

  //     return rOStream;
  //   }
  ///@}

  ///@} addtogroup block

} // namespace Kratos

#endif // define KRATOS_COULOMB_ADHESION_FRICTION_LAW_H_INCLUDED
