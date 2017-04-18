//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_YEOH_MODEL_H_INCLUDED)
#define  KRATOS_YEOH_MODEL_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) YeohModel : public HyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of YeohModel
    KRATOS_CLASS_POINTER_DEFINITION( YeohModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    YeohModel()
      {
      };
    
    /// Copy constructor.
    YeohModel(YeohModel const& rOther)
      {
      };

    /// Assignment operator.
    YeohModel& operator=(YeohModel const& rOther)
      {
	return *this;
      };


    /// Destructor.
    virtual ~YeohModel() {};


    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this yield criterion
     */
    virtual YeohModel::Pointer Clone() const
    {
      YeohModel::Pointer p_clone(new YeohModel(*this));
      return p_clone;
    }


    ///@}
    ///@name Operations
    ///@{
       
    virtual void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
    {
      KRATOS_TRY
	
      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues, Variables);

      rDensityFunction += rVariables.GetMaterialParameters()[0] * ( Variables.Strain.Invariants.I1 - 3.0) + rVariables.GetMaterialParameters()[1] * ( Variables.Strain.Invariants.I2 - 3.0) + rVariables.GetMaterialParameters()[2] * ( Variables.Strain.Invariants.I3 - 3.0);
	
      KRATOS_CATCH(" ")
    };

        
    
    virtual int Check(const Properties& rMaterialProperties,
		      const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      HyperElasticModel::Check(rMaterialProperties,rCurrentProcessInfo);
	
      if( rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS].size() != 3 )
        KRATOS_THROW_ERROR( std::invalid_argument,"HYPERELASTIC_MODEL_PARAMETERS has an invalid size ", "" )

      if(   rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS][0] <= 0.00
	 || rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS][1] >= 0.00
	 || rMaterialProperties[HYPERELASTIC_MODEL_PARAMETERS][2] <= 0.00 )
        KRATOS_THROW_ERROR( std::invalid_argument,"HYPERELASTIC_MODEL_PARAMETERS has an invalid value ", "" )
	  
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
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

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

    //************// dW
    
    virtual double& GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI1
    {
      KRATOS_TRY

      rDerivative  = rVariables.GetMaterialParameters()[0];
      rDerivative += 2.0 * rVariables.GetMaterialParameters()[1] * (rVariables.Strain.Invariants.I1-3);
      rDerivative += 3.0 * rVariables.GetMaterialParameters()[2] * (rVariables.Strain.Invariants.I1-3) * (rVariables.Strain.Invariants.I1-3);

      return rDerivative;

      KRATOS_CATCH(" ")
    };

    virtual double& GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI2
    {
      KRATOS_TRY
	
      rDerivative = 0;

      return rDerivative;

      KRATOS_CATCH(" ")
    };

    virtual double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI3
    {
      KRATOS_TRY
      
      rDerivative  = 0;
      
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
	
      rDerivative  = 0;
      
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

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const
    {
    };

    virtual void load(Serializer& rSerializer)
    {
    };

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class YeohModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_YEOH_MODEL_H_INCLUDED  defined 


