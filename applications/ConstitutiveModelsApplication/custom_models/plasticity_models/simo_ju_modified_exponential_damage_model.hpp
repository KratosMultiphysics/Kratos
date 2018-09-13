//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_SIMO_JU_MODIFIED_EXPONENTIAL_DAMAGE_MODEL_H_INCLUDED )
#define  KRATOS_SIMO_JU_MODIFIED_EXPONENTIAL_DAMAGE_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/damage_model.hpp"
#include "custom_models/plasticity_models/yield_surfaces/simo_ju_yield_surface.hpp"
#include "custom_models/plasticity_models/hardening_rules/modified_exponential_damage_hardening_rule.hpp"
#include "custom_models/elasticity_models/linear_elastic_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SimoJuModifiedExponentialDamageModel : public DamageModel<LinearElasticModel, SimoJuYieldSurface<ModifiedExponentialDamageHardeningRule> >
  {
  public:

    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef LinearElasticModel                             ElasticityModelType;
    typedef ElasticityModelType::Pointer                ElasticityModelPointer;

    //yield surface
    typedef ModifiedExponentialDamageHardeningRule           HardeningRuleType;
    typedef SimoJuYieldSurface<HardeningRuleType>             YieldSurfaceType;
    typedef YieldSurfaceType::Pointer                      YieldSurfacePointer;

    //base type
    typedef DamageModel<ElasticityModelType,YieldSurfaceType>         BaseType;

    //common types
    typedef BaseType::Pointer                         BaseTypePointer;
    typedef BaseType::SizeType                               SizeType;
    typedef BaseType::VoigtIndexType                   VoigtIndexType;
    typedef BaseType::MatrixType                           MatrixType;
    typedef BaseType::ModelDataType                     ModelDataType;
    typedef BaseType::MaterialDataType               MaterialDataType;
    typedef BaseType::PlasticDataType                 PlasticDataType;
    typedef BaseType::InternalVariablesType     InternalVariablesType;


    /// Pointer definition of SimoJuModifiedExponentialDamageModel
    KRATOS_CLASS_POINTER_DEFINITION( SimoJuModifiedExponentialDamageModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimoJuModifiedExponentialDamageModel() : BaseType() {}

    /// Copy constructor.
    SimoJuModifiedExponentialDamageModel(SimoJuModifiedExponentialDamageModel const& rOther)
      :BaseType(rOther) {}

    /// Assignment operator.
    SimoJuModifiedExponentialDamageModel& operator=(SimoJuModifiedExponentialDamageModel const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<SimoJuModifiedExponentialDamageModel>(*this);
    }

    /// Destructor.
    ~SimoJuModifiedExponentialDamageModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    /**
     * Has Values
     */
    bool Has(const Variable<double>& rThisVariable) override
    {
      return false;
    }


    /**
     * Get Values
     */
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
    {

      rValue=0;

      return rValue;
    }

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
      buffer << "SimoJuModifiedExponentialDamageModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "SimoJuModifiedExponentialDamageModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SimoJuModifiedExponentialDamageModel Data";
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


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class SimoJuModifiedExponentialDamageModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_SIMO_JU_MODIFIED_EXPONENTIAL_DAMAGE_MODEL_H_INCLUDED  defined
