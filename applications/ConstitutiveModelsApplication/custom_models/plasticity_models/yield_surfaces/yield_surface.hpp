//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_YIELD_SURFACE_H_INCLUDED )
#define  KRATOS_YIELD_SURFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/hardening_rule.hpp"
#include "custom_utilities/constitutive_model_utilities.hpp"

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
  template<class THardeningRule>
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) YieldSurface
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

    typedef THardeningRule                                       HardeningRuleType;
    typedef typename THardeningRule::PlasticDataType               PlasticDataType;
    typedef typename THardeningRule::InternalVariablesType   InternalVariablesType;

    /// Pointer definition of YieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( YieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    YieldSurface() {}

    /// Copy constructor.
    YieldSurface(YieldSurface const& rOther) : mHardeningRule(rOther.mHardeningRule) {}

    /// Assignment operator.
    YieldSurface& operator=(YieldSurface const& rOther)
    {
      mHardeningRule = rOther.mHardeningRule;
      return *this;
    }

    /// Clone.
    virtual typename YieldSurface::Pointer Clone() const
    {
      return Kratos::make_shared<YieldSurface>(*this);
    }

    /// Destructor.
    virtual ~YieldSurface() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Yield Condition
     */

    virtual double& CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rYieldCondition;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Yield Condition derivative
     */

    virtual double& CalculateDeltaYieldCondition(const PlasticDataType& rVariables, double & rDeltaYieldCondition)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rDeltaYieldCondition;

      KRATOS_CATCH(" ")
    }


    /**
     * Calculate Yield Condition Stresses derivative
     */

    virtual VectorType& CalculateDeltaStressYieldCondition(const PlasticDataType& rVariables, VectorType& rDeltaStressYieldCondition)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rDeltaStressYieldCondition;

      KRATOS_CATCH(" ")
    }


    /**
     * Calculate State Function
     */

    virtual double& CalculateStateFunction(const PlasticDataType& rVariables, double & rStateFunction)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rStateFunction;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate State Function derivative
     */

    virtual double& CalculateDeltaStateFunction(const PlasticDataType& rVariables, double & rDeltaStateFunction)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rDeltaStateFunction;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Plastic Dissipation
     */

    virtual double& CalculatePlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rPlasticDissipation;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Plastic Dissipation derivative
     */

    virtual double& CalculateDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rDeltaPlasticDissipation;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Implex Plastic Dissipation
     */

    virtual double& CalculateImplexPlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rPlasticDissipation;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Implex Plastic Dissipation derivative
     */

    virtual double& CalculateImplexDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the YieldSurface base class ... illegal operation" << std::endl;

      return rDeltaPlasticDissipation;

      KRATOS_CATCH(" ")
    }


    ///@}
    ///@name Access
    ///@{

    HardeningRuleType& GetHardeningRule() { return mHardeningRule; };

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
      buffer << "YieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "YieldSurface";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "YieldSurface Data";
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

    HardeningRuleType mHardeningRule;

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


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;


    virtual void save(Serializer& rSerializer) const
    {
      rSerializer.save("mHardeningRule", mHardeningRule);
    }

    virtual void load(Serializer& rSerializer)
    {
      rSerializer.load("mHardeningRule", mHardeningRule);
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class YieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_YIELD_SURFACE_H_INCLUDED  defined
