//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:               LMonforte $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                     -1.0 $
//
//

#if !defined(KRATOS_VON_MISES_LINEAR_ELASTIC_PLASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_VON_MISES_LINEAR_ELASTIC_PLASTICITY_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/non_linear_associative_plasticity_model.hpp"
#include "custom_models/plasticity_models/yield_surfaces/mises_huber_yield_surface.hpp"
#include "custom_models/plasticity_models/hardening_rules/simo_exponential_hardening_rule.hpp"
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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) VonMisesLinearElasticPlasticityModel : public NonLinearAssociativePlasticityModel<LinearElasticModel, MisesHuberYieldSurface<SimoExponentialHardeningRule> >
  {
  public:
    
    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef LinearElasticModel                             ElasticityModelType;
    typedef typename ElasticityModelType::Pointer       ElasticityModelPointer;

    //yield surface
    typedef SimoExponentialHardeningRule                     HardeningRuleType;
    typedef MisesHuberYieldSurface<HardeningRuleType>         YieldSurfaceType;
    typedef typename YieldSurfaceType::Pointer             YieldSurfacePointer;

    //base type
    typedef NonLinearAssociativePlasticityModel<ElasticityModelType,YieldSurfaceType>  BaseType;

    //common types
    typedef typename BaseType::Pointer                         BaseTypePointer;
    typedef typename BaseType::SizeType                               SizeType;
    typedef typename BaseType::VoigtIndexType                   VoigtIndexType;
    typedef typename BaseType::MatrixType                           MatrixType;
    typedef typename BaseType::ModelDataType                     ModelDataType;
    typedef typename BaseType::MaterialDataType               MaterialDataType;
    typedef typename BaseType::PlasticDataType                 PlasticDataType;
    typedef typename BaseType::InternalVariablesType     InternalVariablesType;


    /// Pointer definition of VonMisesLinearElasticPlasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( VonMisesLinearElasticPlasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VonMisesLinearElasticPlasticityModel() : BaseType() {}
    
    /// Copy constructor.
    VonMisesLinearElasticPlasticityModel(VonMisesLinearElasticPlasticityModel const& rOther) : BaseType(rOther) {}

    /// Assignment operator.
    VonMisesLinearElasticPlasticityModel& operator=(VonMisesLinearElasticPlasticityModel const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return ( VonMisesLinearElasticPlasticityModel::Pointer(new VonMisesLinearElasticPlasticityModel(*this)) );
    }
    
    /// Destructor.
    virtual ~VonMisesLinearElasticPlasticityModel() {}


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
    virtual bool Has(const Variable<double>& rThisVariable) override
    {
      if(rThisVariable == PLASTIC_STRAIN || rThisVariable == DELTA_PLASTIC_STRAIN )
	return true;

      return false;
    }
    

    /**
     * Get Values
     */
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
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

      
      return rValue;
    }
    
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
      buffer << "VonMisesLinearElasticPlasticityModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "VonMisesLinearElasticPlasticityModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "VonMisesLinearElasticPlasticityModel Data";
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

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }
    
    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class VonMisesLinearElasticPlasticityModel

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

#endif // KRATOS_VON_MISES_LINEAR_ELASTIC_PLASTICITY_MODEL_H_INCLUDED  defined 


