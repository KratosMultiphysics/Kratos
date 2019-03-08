//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                February 2019 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CASM_CEMENTED_HARDENING_RULE_H_INCLUDED )
#define      KRATOS_CASM_CEMENTED_HARDENING_RULE_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CasmCementedHardeningRule : public HardeningRule
  {
  protected:

    constexpr static std::size_t VarSize = 11;
    
  public:
    
    typedef InternalVariables<VarSize>   InternalVariablesType;
    typedef PlasticModelData<VarSize>          PlasticDataType;

    /// Pointer definition of CasmCementedHardeningRule
    KRATOS_CLASS_POINTER_DEFINITION( CasmCementedHardeningRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CasmCementedHardeningRule();

    /// Copy constructor.
    CasmCementedHardeningRule(CasmCementedHardeningRule const& rOther);

    /// Assignment operator.
    CasmCementedHardeningRule& operator=(CasmCementedHardeningRule const& rOther);

    /// Clone.
    virtual HardeningRule::Pointer Clone() const override;

    /// Destructor.
    ~CasmCementedHardeningRule();

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

    virtual double& CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening, const MatrixType & rPlasticPotentialDerivative); //do not override -> it must hide the method

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
    virtual std::string Info() const override
    {
      std::stringstream buffer;
      buffer << "CasmCementedHardeningRule" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "CasmCementedHardeningRule";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "CasmCementedHardeningRule Data";
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


    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningRule )
    }

    virtual void load(Serializer& rSerializer) override
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

  }; // Class CasmCementedHardeningRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CASM_CEMENTED_HARDENING_RULE_H_INCLUDED  defined 


