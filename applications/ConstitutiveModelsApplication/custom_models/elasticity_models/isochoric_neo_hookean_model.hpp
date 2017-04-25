//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_ISOCHORIC_NEO_HOOKEAN_MODEL_H_INCLUDED )
#define  KRATOS_ISOCHORIC_NEO_HOOKEAN_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/isochoric_hyperelastic_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IsochoricNeoHookeanModel : public IsochoricHyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsochoricNeoHookeanModel
    KRATOS_CLASS_POINTER_DEFINITION( IsochoricNeoHookeanModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsochoricNeoHookeanModel() : IsochoricHyperElasticModel() {}
    
    /// Copy constructor.
    IsochoricNeoHookeanModel(IsochoricNeoHookeanModel const& rOther) : IsochoricHyperElasticModel(rOther) {}

    /// Assignment operator.
    IsochoricNeoHookeanModel& operator=(IsochoricNeoHookeanModel const& rOther)
    {
      IsochoricHyperElasticModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual ConstitutiveModel::Pointer Clone() const override
    {
      return ( IsochoricNeoHookeanModel::Pointer(new IsochoricNeoHookeanModel(*this)) );      
    }
    
    /// Destructor.
    virtual ~IsochoricNeoHookeanModel() {}


    ///@}
    ///@name Operators
    ///@{

    
    ///@}
    ///@name Operations
    ///@{
  

    // Simplyfied methods must be implemented for performance purposes

    /**
     * Calculate Stresses
     */


    
    /**
     * Calculate Constitutive Components
     */    

    
    
    /**
     * Check
     */    

    virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      std::cout<<" Isochoric NeoHookean Check "<<std::endl;
	
      HyperElasticModel::Check(rMaterialProperties,rCurrentProcessInfo);
	
      if( rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS].size() != 1 )
        KRATOS_ERROR << "HYPERELASTIC_MODEL_PARAMETERS has an invalid size" << std::endl;

      if( rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS][0] <= 0.00 )
        KRATOS_ERROR << "HYPERELASTIC_MODEL_PARAMETERS has an invalid value" << std::endl;

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
        buffer << "IsochoricNeoHookeanModel";
        return buffer.str();
    }
    
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IsochoricNeoHookeanModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IsochoricNeoHookeanModel Data";
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
    
    virtual void CalculateAndAddIsochoricStrainEnergy(HyperElasticDataType& rVariables, double& rIsochoricDensityFunction)
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
	
      rIsochoricDensityFunction += rMaterial.GetModelParameters()[0] * ( rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.I1 - 3.0);
	
      KRATOS_CATCH(" ")
    }
    
    
    virtual void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction)
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
      
      rVolumetricDensityFunction += rMaterial.GetBulkModulus() * 0.25 * ( rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J - 1.0);
      rVolumetricDensityFunction -= rMaterial.GetBulkModulus() * 0.5 * std::log( rVariables.Strain.Invariants.J );
	
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
	
      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    virtual double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI3
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }

    
    virtual double& GetVolumetricFunctionJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //dU/dJ
    {
      KRATOS_TRY
	
      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
      
      rDerivative = 0.5 * rMaterial.GetBulkModulus() * ( rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J - 1.0 );

      rDerivative /= rVariables.Strain.Invariants.J;
      
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
	
      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    }
    

    virtual double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //ddU/dJdJ
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
      
      rDerivative = 0.5 * rMaterial.GetBulkModulus() * (rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J - 1.0 );

      rDerivative /= rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J ;

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


    virtual void save(Serializer& rSerializer) const  override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, IsochoricHyperElasticModel )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, IsochoricHyperElasticModel )      
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IsochoricNeoHookeanModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ISOCHORIC_NEO_HOOKEAN_MODEL_H_INCLUDED  defined 


