//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CAM_CLAY_HARDENING_RULE_H_INCLUDED )
#define  KRATOS_CAM_CLAY_HARDENING_RULE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/hardening_rule.hpp"

namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
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
  /** Detail class definition.
   */
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CamClayHardeningRule : public HardeningRule
  {
  protected:

    constexpr static std::size_t VarSize = 3;

  public:

    typedef InternalVariables<VarSize>   InternalVariablesType;
    typedef PlasticModelData<VarSize>          PlasticDataType;

    /// Pointer definition of CamClayHardeningRule
    KRATOS_CLASS_POINTER_DEFINITION( CamClayHardeningRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CamClayHardeningRule();

    /// Copy constructor.
    CamClayHardeningRule(CamClayHardeningRule const& rOther);

    /// Assignment operator.
    CamClayHardeningRule& operator=(CamClayHardeningRule const& rOther);

    /// Clone.
    HardeningRule::Pointer Clone() const override;

    /// Destructor.
    ~CamClayHardeningRule() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /**
     * Calculate Hardening functions
     */

    virtual double& CalculateHardening(const PlasticDataType& rVariables, double& rHardening); //do not override -> it must hide the method

    /**
     * Calculate Hardening function derivatives
     */

    virtual double& CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening); //do not override -> it must hide the method


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
      buffer << "CamClayHardeningRule" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "CamClayHardeningRule";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "CamClayHardeningRule Data";
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

    using HardeningRule::CalculateHardening;
    using HardeningRule::CalculateDeltaHardening;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;


    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningRule )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningRule )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class CamClayHardeningRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CAM_CLAY_HARDENING_RULE_H_INCLUDED  defined


