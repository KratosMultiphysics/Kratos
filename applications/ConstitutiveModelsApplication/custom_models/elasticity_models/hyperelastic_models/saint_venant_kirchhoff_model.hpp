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

      //Calculate HyperElastic Saint Venant Kirchhoff density function
      const MaterialDataType& rMaterial = rValues.GetMaterialParameters();

      // Lame constants
      const double& rLameLambda = rMaterial.GetLameLambda();
      const double& rLameMu     = rMaterial.GetLameMu();

      double trace = (Variables.Strain.CauchyGreenMatrix(0,0)+Variables.Strain.CauchyGreenMatrix(1,1)+Variables.Strain.CauchyGreenMatrix(2,2));

      rDensityFunction = 0.5*rLameLambda*trace;
      
      trace = Variables.Strain.CauchyGreenMatrix(0,0)*Variables.Strain.CauchyGreenMatrix(0,0)
	    + Variables.Strain.CauchyGreenMatrix(0,1)*Variables.Strain.CauchyGreenMatrix(1,0)
	    + Variables.Strain.CauchyGreenMatrix(0,2)*Variables.Strain.CauchyGreenMatrix(2,0)
	    + Variables.Strain.CauchyGreenMatrix(1,0)*Variables.Strain.CauchyGreenMatrix(0,1)
	    + Variables.Strain.CauchyGreenMatrix(1,1)*Variables.Strain.CauchyGreenMatrix(1,1)
	    + Variables.Strain.CauchyGreenMatrix(1,2)*Variables.Strain.CauchyGreenMatrix(2,1)
	    + Variables.Strain.CauchyGreenMatrix(2,0)*Variables.Strain.CauchyGreenMatrix(0,2)
	    + Variables.Strain.CauchyGreenMatrix(2,1)*Variables.Strain.CauchyGreenMatrix(1,2)
	    + Variables.Strain.CauchyGreenMatrix(2,2)*Variables.Strain.CauchyGreenMatrix(2,2);
      
      trace *= trace;
      
      rDensityFunction += rLameMu*trace;
      
      KRATOS_CATCH(" ")
    }


    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      // bounded_matrix<double,6,6> ConstitutiveTensor;
      // this->CalculateAndAddConstitutiveTensor(Variables,ConstitutiveTensor);

      // VectorType StrainVector;
      // StrainVector = ConstitutiveLawUtilities::StrainTensorToVector(Variables.Strain.CauchyGreenMatrix,StrainVector);

      // VectorType StressVector;
      // this->CalculateAndAddStressTensor(Variables,ConstitutiveTensor,StrainVector,StressVector);

      // rStressMatrix = ConstitutiveLawUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);
      
      this->CalculateAndAddStressTensor(Variables,rStressMatrix);
      
      const StressMeasureType& rStressMeasure = rValues.GetStressMeasure();
   
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){
	
	const MatrixType& rDeformationGradientF = rValues.GetDeformationGradientF();

	//Variables.Strain.InverseCauchyGreenMatrix used as an auxiliar matrix (contravariant push forward)
	noalias( Variables.Strain.InverseCauchyGreenMatrix ) = prod( trans(rDeformationGradientF), rStressMatrix );
	noalias( rStressMatrix )  = prod( Variables.Strain.InverseCauchyGreenMatrix, rDeformationGradientF );
	
      }
      
      
      KRATOS_CATCH(" ")
    }
    
    void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix)
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      bounded_matrix<double,6,6> ConstitutiveTensor;
      this->CalculateAndAddConstitutiveTensor(Variables,ConstitutiveTensor);
      
      rConstitutiveMatrix = ConstitutiveLawUtilities::ConstitutiveTensorToMatrix(ConstitutiveTensor,rConstitutiveMatrix);

      // if StressMeasure_Kirchhoff, a push forward of the ConstitutiveMatrix must be done, but it is avoided
      // it is computationally expensive but not relevant for the convegence of the method
	            
      KRATOS_CATCH(" ")
    }

    
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY
     
      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      bounded_matrix<double,6,6> ConstitutiveTensor;
      this->CalculateAndAddConstitutiveTensor(Variables,ConstitutiveTensor);

      VectorType StrainVector;
      StrainVector = ConstitutiveLawUtilities::StrainTensorToVector(Variables.Strain.CauchyGreenMatrix,StrainVector);
    
      VectorType StressVector;
      this->CalculateAndAddStressTensor(Variables,ConstitutiveTensor,StrainVector,StressVector);

      rStressMatrix = ConstitutiveLawUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);
      rConstitutiveMatrix = ConstitutiveLawUtilities::ConstitutiveTensorToMatrix(ConstitutiveTensor,rConstitutiveMatrix);

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

    void CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables)
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
	  ConstitutiveLawUtilities::CalculateGreenLagrangeStrain( rDeformationGradientF, rVariables.Strain.CauchyGreenMatrix);
	  rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
	}
	else{
	  KRATOS_ERROR << "calling initialize HyperElasticModel .. StrainMeasure provided is inconsistent" << std::endl;
	}

	rValues.StrainMatrix = rVariables.Strain.CauchyGreenMatrix;
	
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //mCauchyGreenMatrix = GreenLagrangeTensor

	if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_Left ){
	  ConstitutiveLawUtilities::LeftCauchyToAlmansiStrain( rStrainMatrix , rVariables.Strain.CauchyGreenMatrix);

	  //rVariables.Strain.InverseCauchyGreenMatrix used as an auxiliar matrix (covariant pull back)
	  noalias( rVariables.Strain.InverseCauchyGreenMatrix ) = prod( trans(rDeformationGradientF), rVariables.Strain.CauchyGreenMatrix );
	  noalias( rVariables.Strain.CauchyGreenMatrix)  = prod( rVariables.Strain.InverseCauchyGreenMatrix, rDeformationGradientF );

	  rValues.State.Set(ConstitutiveModelData::COMPUTED_STRAIN);
	}
	else if( rStrainMeasure == ConstitutiveModelData::CauchyGreen_None ){
	  ConstitutiveLawUtilities::CalculateGreenLagrangeStrain( rDeformationGradientF, rVariables.Strain.CauchyGreenMatrix);
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


    void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, bounded_matrix<double,6,6>& rConstitutiveTensor, VectorType& rStrainVector, VectorType& rStressVector)
    {
      KRATOS_TRY

      noalias(rStressVector) = prod(rConstitutiveTensor,rStrainVector);
      
      rVariables.State().Set(ConstitutiveModelData::COMPUTED_STRESS);
    
      KRATOS_CATCH(" ")
    }
    
    
    void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      const ModelDataType&  rModelData  = rVariables.GetModelData();
      const MaterialDataType& rMaterial = rModelData.GetMaterialParameters();

      // Lame constants
      const double& rLameLambda = rMaterial.GetLameLambda();
      const double& rLameMu     = rMaterial.GetLameMu();
	
      rStressMatrix  = rVariables.Strain.CauchyGreenMatrix;
      rStressMatrix *= 2.0 * rLameMu;
      
      double trace = (rVariables.Strain.CauchyGreenMatrix(0,0)+rVariables.Strain.CauchyGreenMatrix(1,1)+rVariables.Strain.CauchyGreenMatrix(2,2));
      trace *= rLameLambda;

      rStressMatrix(0,0) += trace;
      rStressMatrix(1,1) += trace;
      rStressMatrix(2,2) += trace;
	    
      rVariables.State().Set(ConstitutiveModelData::COMPUTED_STRESS);
    
      KRATOS_CATCH(" ")
    }

    
    void CalculateAndAddConstitutiveTensor(HyperElasticDataType& rVariables, bounded_matrix<double,6,6>& rConstitutiveTensor)
    {
      KRATOS_TRY
              
      //Calculate HyperElastic ConstitutiveMatrix
      const ModelDataType&  rModelData  = rVariables.GetModelData();
      const MaterialDataType& rMaterial = rModelData.GetMaterialParameters();

      // Lame constants
      const double& rYoungModulus       = rMaterial.GetYoungModulus();
      const double& rPoissonCoefficient = rMaterial.GetPoissonCoefficient();

      rConstitutiveTensor.clear();
      
      // 3D linear elastic constitutive matrix
      rConstitutiveTensor ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
      rConstitutiveTensor ( 1 , 1 ) = rConstitutiveTensor ( 0 , 0 );
      rConstitutiveTensor ( 2 , 2 ) = rConstitutiveTensor ( 0 , 0 );

      rConstitutiveTensor ( 3 , 3 ) = rConstitutiveTensor ( 0 , 0 )*(1.0-2.0*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));
      rConstitutiveTensor ( 4 , 4 ) = rConstitutiveTensor ( 3 , 3 );
      rConstitutiveTensor ( 5 , 5 ) = rConstitutiveTensor ( 3 , 3 );

      rConstitutiveTensor ( 0 , 1 ) = rConstitutiveTensor ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
      rConstitutiveTensor ( 1 , 0 ) = rConstitutiveTensor ( 0 , 1 );

      rConstitutiveTensor ( 0 , 2 ) = rConstitutiveTensor ( 0 , 1 );
      rConstitutiveTensor ( 2 , 0 ) = rConstitutiveTensor ( 0 , 1 );

      rConstitutiveTensor ( 1 , 2 ) = rConstitutiveTensor ( 0 , 1 );
      rConstitutiveTensor ( 2 , 1 ) = rConstitutiveTensor ( 0 , 1 );

    
      rVariables.State().Set(ConstitutiveModelData::COMPUTED_CONSTITUTIVE_MATRIX);

    
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


