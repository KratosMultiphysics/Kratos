//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_COMPRESSIBLE_NEOHOOKEAN_MODEL_H_INCLUDED)
#define  KRATOS_COMPRESSIBLE_NEOHOOKEAN_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/hyperelastic_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CompressibleNeoHookeanModel : public HyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of CompressibleNeoHookeanModel
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNeoHookeanModel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CompressibleNeoHookeanModel() : HyperElasticModel() {}
    
    /// Copy constructor.
    CompressibleNeoHookeanModel(CompressibleNeoHookeanModel const& rOther) : HyperElasticModel(rOther) {}
    
    /// Assignment operator.
    CompressibleNeoHookeanModel& operator=(CompressibleNeoHookeanModel const& rOther)
    {
      HyperElasticModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual ConstitutiveModel::Pointer Clone() const override
    {
      return ( CompressibleNeoHookeanModel::Pointer(new CompressibleNeoHookeanModel(*this)) );      
    }
    
    /// Destructor.
    virtual ~CompressibleNeoHookeanModel() {}


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
		      const ProcessInfo& rCurrentProcessInfo) override
    { 
      KRATOS_TRY

      HyperElasticModel::Check(rMaterialProperties,rCurrentProcessInfo);
	
      if( C10.Key() == 0 )
	KRATOS_ERROR << "C10 has an invalid key or size" << std::endl;

      if( C10 <= 0.00 )
        KRATOS_ERROR << "C10 has an invalid value" << std::endl;

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
        buffer << "CompressibleNeoHookeanModel";
        return buffer.str();
    }
    
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CompressibleNeoHookeanModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "CompressibleNeoHookeanModel Data";
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


    // specialized methods:
    
    // virtual void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override
    // {
    //   KRATOS_TRY

    //   const ModelDataType&  rModelData        = rVariables.GetModelData();
    //   const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
    
    //   MatrixType StressMatrix;
      
    //   const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
      
    //   if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mStrainMatrix = RightCauchyGreen (C)

    // 	StressMatrix  = rMaterial.GetLameLambda() * std::log( rVariables.Strain.Invariants.J ) * rVariables.Strain.InverseMatrix;
    //     StressMatrix += rMaterial.GetLameMu() * ( msIdentityMatrix - rVariables.Strain.InverseMatrix);

    // 	rStressMatrix += StressMatrix;
      
    //   }
    //   else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mStrainMatrix = LeftCauchyGreen (b)

    // 	StressMatrix  = rMaterial.GetLameLambda() * std::log( rVariables.Strain.Invariants.J ) * msIdentityMatrix;
    //     StressMatrix += rMaterial.GetLameMu() * ( rVariables.Strain.Matrix - msIdentityMatrix );
	
    // 	rStressMatrix += StressMatrix;      
    //   }

    
    //   rVariables.State().Set(ConstitutiveModelData::COMPUTED_STRESS);
    
    //   KRATOS_CATCH(" ")
    // }

    
    // virtual double& AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
    // 					     const unsigned int& a, const unsigned int& b,
    // 					     const unsigned int& c, const unsigned int& d) override
    // {
    //   KRATOS_TRY
	
    
    //   double Cabcd = 0;
 
    //   const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

    //   const ModelDataType&  rModelData        = rVariables.GetModelData();
    //   const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
         
    //   if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mStrainMatrix = RightCauchyGreen (C)
    // 	Cabcd  = rMaterial.GetLameLambda() * (rVariables.Strain.InverseMatrix(a,b)*rVariables.Strain.InverseMatrix(c,d));
    // 	Cabcd += (rMaterial.GetLameMu() - rMaterial.GetLameLambda() * std::log(rVariables.Strain.Invariants.J)) * (rVariables.Strain.InverseMatrix(a,c)*rVariables.Strain.InverseMatrix(b,d)+rVariables.Strain.InverseMatrix(a,d)*rVariables.Strain.InverseMatrix(b,c));
    //   }
    //   else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mStrainMatrix = LeftCauchyGreen (b)
    // 	Cabcd  = rMaterial.GetLameLambda() * (msIdentityMatrix(a,b)*msIdentityMatrix(c,d));
    // 	Cabcd += (rMaterial.GetLameMu() - rMaterial.GetLameLambda() * std::log(rVariables.Strain.Invariants.J)) * (msIdentityMatrix(a,c)*msIdentityMatrix(b,d)+msIdentityMatrix(a,d)*msIdentityMatrix(b,c));  
    //   }
      
    //   rCabcd += Cabcd;
    
    //   rVariables.State().Set(ConstitutiveModelData::COMPUTED_CONSTITUTIVE_MATRIX);
    
    //   return rCabcd;
    
    //   KRATOS_CATCH(" ")
    // }
    

    //************// W
    
    virtual void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction)
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //g(J) = (lambda/2)*ln(J)² - (mu)*lnJ 
      rVolumetricDensityFunction = rMaterial.GetLameLambda() * 0.5 * std::log( rVariables.Strain.Invariants.J );
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

      //derivative of "g(J) = (lambda/2)*ln(J)² - (mu)*lnJ"
      //dg(J)/dI3 = (lambda/2)*(lnJ/J²) - mu*(1/J²)/2
      rDerivative  = 0.5 * rMaterial.GetLameLambda() * std::log( rVariables.Strain.Invariants.J );
      rDerivative -= 0.5 * rMaterial.GetLameMu();
      rDerivative /= rVariables.Strain.Invariants.I3;
      
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
      //ddg(J)/dI3dI3 = (lambda/4)*(1-lnJ)/J⁴ + (mu/2)*(1/J⁴)
      rDerivative  = 0.25 * rMaterial.GetLameLambda() * (1.0 - 2.0 *  std::log( rVariables.Strain.Invariants.J ) );
      rDerivative += 0.5 * rMaterial.GetLameMu();
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

  }; // Class CompressibleNeoHookeanModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_NEOHOOKEAN_MODEL_H_INCLUDED  defined 


