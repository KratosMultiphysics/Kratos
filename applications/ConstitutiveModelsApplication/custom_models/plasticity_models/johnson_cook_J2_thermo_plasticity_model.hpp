//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_JOHNSON_COOK_J2_THERMO_PLASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_JOHNSON_COOK_J2_THERMO_PLASTICITY_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/non_linear_rate_dependent_plasticity_model.hpp"
#include "custom_models/plasticity_models/yield_surfaces/mises_huber_thermal_yield_surface.hpp"
#include "custom_models/plasticity_models/hardening_rules/johnson_cook_thermal_hardening_rule.hpp"
#include "custom_models/elasticity_models/incompressible_neo_hookean_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) JohnsonCookJ2ThermoPlasticityModel : public NonLinearRateDependentPlasticityModel<IncompressibleNeoHookeanModel, MisesHuberThermalYieldSurface<JohnsonCookThermalHardeningRule> >
  {
  public:

    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef IncompressibleNeoHookeanModel                  ElasticityModelType;
    typedef ElasticityModelType::Pointer                ElasticityModelPointer;

    //yield surface
    typedef JohnsonCookThermalHardeningRule                    HardeningRuleType;
    typedef MisesHuberThermalYieldSurface<HardeningRuleType>    YieldSurfaceType;
    typedef YieldSurfaceType::Pointer                        YieldSurfacePointer;

    //base type
    typedef NonLinearRateDependentPlasticityModel<ElasticityModelType,YieldSurfaceType>  BaseType;

    //common types
    typedef BaseType::Pointer                         BaseTypePointer;
    typedef BaseType::SizeType                               SizeType;
    typedef BaseType::VoigtIndexType                   VoigtIndexType;
    typedef BaseType::MatrixType                           MatrixType;
    typedef BaseType::ModelDataType                     ModelDataType;
    typedef BaseType::MaterialDataType               MaterialDataType;
    typedef BaseType::PlasticDataType                 PlasticDataType;
    typedef BaseType::InternalVariablesType     InternalVariablesType;


    /// Pointer definition of JohnsonCookJ2ThermoPlasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( JohnsonCookJ2ThermoPlasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    JohnsonCookJ2ThermoPlasticityModel() : BaseType() {}

    /// Copy constructor.
    JohnsonCookJ2ThermoPlasticityModel(JohnsonCookJ2ThermoPlasticityModel const& rOther) : BaseType(rOther) {}

    /// Assignment operator.
    JohnsonCookJ2ThermoPlasticityModel& operator=(JohnsonCookJ2ThermoPlasticityModel const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<JohnsonCookJ2ThermoPlasticityModel>(*this);
    }

    /// Destructor.
    ~JohnsonCookJ2ThermoPlasticityModel() override {}


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
     * method to ask the constitutive model the list of variables (dofs) needed from the domain
     * @param rScalarVariables : list of scalar dofs
     * @param rComponentVariables :  list of vector dofs
     */
    void GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
					std::vector<Variable<array_1d<double,3> > >& rComponentVariables) override
    {
      KRATOS_TRY

      PlasticityModel::GetDomainVariablesList(rScalarVariables, rComponentVariables);

      rScalarVariables.push_back(TEMPERATURE);

      KRATOS_CATCH(" ")
    }


    /**
     * Has Values
     */
    bool Has(const Variable<double>& rThisVariable) override
    {
      if(rThisVariable == PLASTIC_STRAIN || rThisVariable == DELTA_PLASTIC_STRAIN )
	return true;

      if(rThisVariable == PLASTIC_DISSIPATION || rThisVariable == DELTA_PLASTIC_DISSIPATION )
	return true;

      return false;
    }


    /**
     * Get Values
     */
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
    {

      rValue=0;

      if (rThisVariable==PLASTIC_STRAIN)
	{
	  rValue = this->mInternal.Variables[0];
	}


      if (rThisVariable==DELTA_PLASTIC_STRAIN)
	{
	  rValue = this->mInternal.Variables[0]-mPreviousInternal.Variables[0];
	}

      if (rThisVariable==PLASTIC_DISSIPATION)
	{
	  rValue = this->mThermalVariables.PlasticDissipation;
	}


      if (rThisVariable==DELTA_PLASTIC_DISSIPATION)
	{
	  rValue = this->mThermalVariables.DeltaPlasticDissipation;
	}

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
      buffer << "JohnsonCookJ2ThermoPlasticityModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "JohnsonCookJ2ThermoPlasticityModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "JohnsonCookJ2ThermoPlasticityModel Data";
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

  }; // Class JohnsonCookJ2ThermoPlasticityModel

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

#endif // KRATOS_JOHNSON_COOK_J2_THERMO_PLASTICITY_MODEL_H_INCLUDED  defined
