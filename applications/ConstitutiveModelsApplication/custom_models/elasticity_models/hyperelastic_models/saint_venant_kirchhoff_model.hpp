//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_SAINT_VENANT_KIRCHHOFF_MODEL_H_INCLUDED )
#define  KRATOS_SAINT_VENANT_KIRCHHOFF_MODEL_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SaintVenantKirchhoffModel : public HyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of SaintVenantKirchhoffModel
    KRATOS_CLASS_POINTER_DEFINITION( SaintVenantKirchhoffModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SaintVenantKirchhoffModel() : HyperElasticModel() {}
    
    /// Copy constructor.
    SaintVenantKirchhoffModel(SaintVenantKirchhoffModel const& rOther) : HyperElasticModel(rOther) {}

    /// Assignment operator.
    SaintVenantKirchhoffModel& operator=(SaintVenantKirchhoffModel const& rOther)
    {
	HyperElasticModel::operator=(rOther);
	return *this;
    }

    /// Clone.
    virtual ElasticityModel::Pointer Clone() const override
    {
      return ( SaintVenantKirchhoffModel::Pointer(new SaintVenantKirchhoffModel(*this)) );      
    }
 
    /// Destructor.
    virtual ~SaintVenantKirchhoffModel() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
  

    virtual void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction) override
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues, Variables);

      rDensityFunction = 0;
      this->CalculateAndAddStrainEnergy( Variables, rDensityFunction );

	
      KRATOS_CATCH(" ")
    }


    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      this->CalculateAndAddConstiutiveTensor(Variables);

      VectorType StrainVector;
      StrainVector = ConstitutiveLawUtilities::StrainTensorToVector(Variables.CauchyGreenMatrix, StrainVector);

      VectorType StressVector;
      this->CalculateAndAddStressTensor(Variables, StrainVector, StressVector);

      rStressMatrix = ConstitutiveLawUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);
      
      
      KRATOS_CATCH(" ")
    }
    
    
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY
     
      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      VectorType StrainVector;
      StrainVector = ConstitutiveLawUtilities::StrainTensorToVector(Variables.CauchyGreenMatrix, StrainVector);
    
      VectorType StressVector;
      this->CalculateAndAddStressTensor(Variables,StrainVector,StressVector);

      rStressMatrix = ConstitutiveLawUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

    
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
        buffer << "SaintVenantKirchhoffModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SaintVenantKirchhoffModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SaintVenantKirchhoffModel Data";
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

    void HyperElasticModel::CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables)
    {
      KRATOS_TRY

      //set model data pointer
      rVariables.SetModelData(rValues);
      rVariables.SetState(rValues.State);
    
      //cauchy green tensor
      const MatrixType& rStrainMatrix         = rValues.GetStrainMatrix();
      const MatrixType& rDeformationGradientF = rValues.GetDeformationGradientF();
      
      const StrainMeasureType& rStrainMeasure = rValues.GetStrainMeasure();
      const StressMeasureType& rStressMeasure = rValues.GetStressMeasure();
    
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //mCauchyGreenMatrix = GreenLagrangeTensor
	
	if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_Right ){
	  ConstitutiveLawUtilities::RightCauchyToGreenLagrangeStrain( rStrainMatrix, rVariables.Strain.CauchyGreenMatrix);  
	  rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
	}
	else if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_None ){
	  ConstitutiveLawUtilities::CalculateGreenLagrangeStrain( rDeformationGradient, rVariables.Strain.CauchyGreenMatrix);
	  rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
	}
	else{
	  KRATOS_ERROR << "calling initialize HyperElasticModel .. StrainMeasure provided is inconsistent" << std::endl;
	}

	rStrainMatrix = rVariables.Strain.CauchyGreenMatrix;
	
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mCauchyGreenMatrix = GreenLagrangeTensor

	if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_Left ){
	  ConstitutiveLawUtilities::LeftCauchyGreenToAlmansiStrain( rStrainMatrix , rVariables.Strain.CauchyGreenMatrix);

	  //rVariables.Strain.InverseCauchyGreenMatrix used as an auxiliar matrix
	  noalias( rVariables.Strain.InverseCauchyGreenMatrix ) = prod( trans(rDeformationGradientF), rVariables.Strain.CauchyGreenMatrix );
	  noalias( rVariables.Strain.CauchyGreenMatrix)  = prod( rVariables.Strain.InverseCauchyGreenMatrix, rDeformationGradientF );

	  rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
	}
	else if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_None ){
	  ConstitutiveLawUtilities::CalculateGreenLagrangeStrain( rDeformationGradient, rVariables.Strain.CauchyGreenMatrix);
	  rValues.StrainMatrix = rVariables.Strain.CauchyGreenMatrix;
	  rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
	}
	else{
	  KRATOS_ERROR << "calling initialize HyperElasticModel .. StrainMeasure provided is inconsistent" << std::endl;
	}
      
      }
      else{
	KRATOS_ERROR << "calling initialize HyperElasticModel .. StressMeasure required is inconsistent"  << std::endl;
      }

      
      KRATOS_CATCH(" ")
    }

    
    void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, VectorType& rStressVector)
    {
      KRATOS_TRY

      noalias(rStressVector) = prod(rVariables.ConstitutiveMatrix,rStrainVector);
      
      rVariables.State().Set(ConstitutiveModelData::COMPUTED_STRESS);
    
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

  }; // Class SaintVenantKirchhoffModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SAINT_VENANT_KIRCHHOFF_MODEL_H_INCLUDED  defined 


