//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NEO_HOOKEAN_MODEL_H_INCLUDED)
#define  KRATOS_NEO_HOOKEAN_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/hyperelastic_models/hyperelastic_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NeoHookeanModel : public HyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of NeoHookeanModel
    KRATOS_CLASS_POINTER_DEFINITION(NeoHookeanModel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NeoHookeanModel() : HyperElasticModel() {}
    
    /// Copy constructor.
    NeoHookeanModel(NeoHookeanModel const& rOther) : HyperElasticModel(rOther) {}
    
    /// Assignment operator.
    NeoHookeanModel& operator=(NeoHookeanModel const& rOther)
    {
      HyperElasticModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual ElasticityModel::Pointer Clone() const override
    {
      return ( NeoHookeanModel::Pointer(new NeoHookeanModel(*this)) );      
    }
    
    /// Destructor.
    virtual ~NeoHookeanModel() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
       
    virtual void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
    {
      KRATOS_TRY
	
      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues, Variables);

      this->CalculateAndAddVolumetricStrainEnergy( Variables, rDensityFunction );
      
      rDensityFunction += rValues.MaterialParameters.GetModelParameters()[0] * ( Variables.Strain.Invariants.I1 - 3.0);
	
      KRATOS_CATCH(" ")
    }

       
    
    virtual int Check(const Properties& rMaterialProperties,
		      const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      HyperElasticModel::Check(rMaterialProperties,rCurrentProcessInfo);
	
      if( rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS].size() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"HYPERELASTIC_MODEL_PARAMETERS has an invalid size ", rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS].size() )

      if( rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS][0] <= 0.00 )
        KRATOS_THROW_ERROR( std::invalid_argument,"HYPERELASTIC_MODEL_PARAMETERS has an invalid value ", "" )

      return 0;
	  
      KRATOS_CATCH(" ")	  
    }
    
    
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
        buffer << "NeoHookeanModel";
        return buffer.str();
    }
    
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NeoHookeanModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NeoHookeanModel Data";
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

    //************// W
    
    virtual void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction)
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
	
      rVolumetricDensityFunction  = rMaterial.GetLameLambda() * 0.25 * ( rVariables.Strain.Invariants.I3 - 1.0);
      rVolumetricDensityFunction -= rMaterial.GetLameLambda() * 0.5 * std::log( rVariables.Strain.Invariants.J );
      rVolumetricDensityFunction -= rMaterial.GetLameMu() * std::log( rVariables.Strain.Invariants.J );
	
      KRATOS_CATCH(" ")
    }
    
    //************// dW
    
    virtual double& GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI1
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
      rDerivative = rMaterial.GetModelParameters()[0];

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI2
    {
      KRATOS_TRY
	
      rDerivative =  0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI3
    {
      KRATOS_TRY

      const MaterialDataType&  rMaterial = rVariables.GetMaterialParameters();
      
      rDerivative  = rMaterial.GetLameLambda();
      rDerivative  = 2.0 * rMaterial.GetLameMu();
      rDerivative /= -rVariables.Strain.Invariants.I3;
      rDerivative += rMaterial.GetLameLambda();
      rDerivative *= 0.25;
      
      return rDerivative;

      KRATOS_CATCH(" ")
    }


    virtual double& GetFunction2ndI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI1dI1
    {
      KRATOS_TRY
	
      rDerivative = 0.0;
      
      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction2ndI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI2dI2
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI3dI3
    {
      KRATOS_TRY

      const MaterialDataType&   rMaterial = rVariables.GetMaterialParameters();
	
      rDerivative  = 0.25 * rMaterial.GetLameLambda();
      rDerivative  = 0.5  * rMaterial.GetLameMu();
      rDerivative /= (rVariables.Strain.Invariants.I3 * rVariables.Strain.Invariants.I3);
      
      return rDerivative;

      KRATOS_CATCH(" ")
    }
    
    
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

    
    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticModel )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticModel )      
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NeoHookeanModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEO_HOOKEAN_MODEL_H_INCLUDED  defined 


