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
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLaw());
     *      return p_clone;
     */
    virtual FrictionLaw::Pointer Clone() const;

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
    //virtual std::string Info() const;

    /// Print information about this object.
    //virtual void PrintInfo(std::ostream& rOStream) const;
    
    /// Print object's data.
    //virtual void PrintData(std::ostream& rOStream) const;
    

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
    
    virtual double EvaluateHardening( const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables);

    virtual double EvaluateContactYield( const double& rTangentStress, const double& rNormalStress, const double& rPlasticSlip, FrictionLawVariables& rTangentVariables);

    virtual void EvaluateYieldDerivativeRespectStress( double& rdF_dt, double & rdF_dp, const double& rTangentStress, const double& rNormalStress, const double& Gamma, FrictionLawVariables& rTangentVariables);


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

    virtual void save( Serializer& rSerializer ) const
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FrictionLaw )
    }

    virtual void load( Serializer& rSerializer )
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
