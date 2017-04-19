//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_MODEL_H_INCLUDED )
#define  KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_models/elasticity_models/hyperelastic_models/isochoric_neo_hookean_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IncompressibleNeoHookeanModel : public IsochoricNeoHookeanModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IncompressibleNeoHookeanModel
    KRATOS_CLASS_POINTER_DEFINITION( IncompressibleNeoHookeanModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IncompressibleNeoHookeanModel() {}
    
    /// Copy constructor.
    IncompressibleNeoHookeanModel(IncompressibleNeoHookeanModel const& rOther) {}

    /// Assignment operator.
    IncompressibleNeoHookeanModel& operator=(IncompressibleNeoHookeanModel const& rOther) {return *this;}

    /// Clone.
    virtual ElasticityModel::Pointer Clone() const
    {
      return ( IncompressibleNeoHookeanModel::Pointer(new IncompressibleNeoHookeanModel(*this)) );      
    }
    
    /// Destructor.
    virtual ~IncompressibleNeoHookeanModel() {}


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
    virtual int Check(const Properties& rMaterialProperties,
		      const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      HyperElasticModel::Check(rMaterialProperties,rCurrentProcessInfo);
	
      if( rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS].size() != 1 )
        KRATOS_THROW_ERROR( std::invalid_argument,"HYPERELASTIC_MODEL_PARAMETERS has an invalid size ", "" )

      if( rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS][0] <= 0.00 )
        KRATOS_THROW_ERROR( std::invalid_argument,"HYPERELASTIC_MODEL_PARAMETERS has an invalid value ", "" )

      return 0;
	  
      KRATOS_CATCH(" ")	  
    };
    
    
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
        buffer << "IncompressibleHyperElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IncompressibleHyperElasticModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IncompressibleHyperElasticModel Data";
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


    
    //************// dW
    
    virtual double& GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI1
    {
      KRATOS_TRY
	
      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      rDerivative = rMaterial.GetModelParameters()[0];

      return rDerivative;

      KRATOS_CATCH(" ")
    };

    virtual double& GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI2
    {
      KRATOS_TRY
	
      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    };

    virtual double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI3
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    };

    
    virtual double& GetVolumetricFunctionJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //dU/dJ
    {
      KRATOS_TRY

      const ModelDataType&  rValues = rVariables.GetModelData();
	
      rDerivative = rValues.GetDeterminantF();
      
      return rDerivative;

      KRATOS_CATCH(" ")
    };


    virtual double& GetFunction2ndI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI1dI1
    {
      KRATOS_TRY
	
      rDerivative = 0.0;
      
      return rDerivative;

      KRATOS_CATCH(" ")
    };

    virtual double& GetFunction2ndI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI2dI2
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    };

    virtual double& GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI3dI3
    {
      KRATOS_TRY
	
      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    };
    

    virtual double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //ddU/dJdJ
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    };

    
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, IsochoricNeoHookeanModel )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, IsochoricNeoHookeanModel )      
    }

 
    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IncompressibleNeoHookeanModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_MODEL_H_INCLUDED  defined 


