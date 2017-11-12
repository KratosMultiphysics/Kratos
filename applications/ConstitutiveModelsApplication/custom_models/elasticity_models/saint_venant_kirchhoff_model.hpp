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
    virtual ConstitutiveModel::Pointer Clone() const override
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

      double trace = (Variables.Strain.Matrix(0,0)+Variables.Strain.Matrix(1,1)+Variables.Strain.Matrix(2,2));

      rDensityFunction = 0.5*rLameLambda*trace;
      
      trace = Variables.Strain.Matrix(0,0)*Variables.Strain.Matrix(0,0)
	    + Variables.Strain.Matrix(0,1)*Variables.Strain.Matrix(1,0)
	    + Variables.Strain.Matrix(0,2)*Variables.Strain.Matrix(2,0)
	    + Variables.Strain.Matrix(1,0)*Variables.Strain.Matrix(0,1)
	    + Variables.Strain.Matrix(1,1)*Variables.Strain.Matrix(1,1)
	    + Variables.Strain.Matrix(1,2)*Variables.Strain.Matrix(2,1)
	    + Variables.Strain.Matrix(2,0)*Variables.Strain.Matrix(0,2)
	    + Variables.Strain.Matrix(2,1)*Variables.Strain.Matrix(1,2)
	    + Variables.Strain.Matrix(2,2)*Variables.Strain.Matrix(2,2);
      
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
      // this->CalculateAndAddConstitutiveMatrix(Variables,ConstitutiveTensor);

      // VectorType StrainVector;
      // StrainVector = ConstitutiveModelUtilities::StrainTensorToVector(Variables.Strain.Matrix,StrainVector);

      // VectorType StressVector;
      // this->CalculateAndAddStressTensor(Variables,ConstitutiveTensor,StrainVector,StressVector);

      // rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);
      
      this->CalculateAndAddStressTensor(Variables,rStressMatrix);
      
      const StressMeasureType& rStressMeasure = rValues.GetStressMeasure();
   
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){
	
	const MatrixType& rTotalDeformationMatrix = rValues.GetTotalDeformationMatrix();

	//Variables.Strain.InverseMatrix used as an auxiliar matrix (contravariant push forward)
	noalias( Variables.Strain.InverseMatrix ) = prod( trans(rTotalDeformationMatrix), rStressMatrix );
	noalias( rStressMatrix )  = prod( Variables.Strain.InverseMatrix, rTotalDeformationMatrix );
	
      }
      
      
      KRATOS_CATCH(" ")
    }
    
    void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      bounded_matrix<double,6,6> ConstitutiveTensor;
      this->CalculateAndAddConstitutiveMatrix(Variables,ConstitutiveTensor);
      
      rConstitutiveMatrix = ConstitutiveModelUtilities::ConstitutiveTensorToMatrix(ConstitutiveTensor,rConstitutiveMatrix);

      // if StressMeasure_Kirchhoff, a push forward of the ConstitutiveMatrix must be done, but it is avoided
      // it is computationally expensive but not relevant for the convergence of the method
	            
      KRATOS_CATCH(" ")
    }

    
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY
     
      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      bounded_matrix<double,6,6> ConstitutiveTensor;
      this->CalculateAndAddConstitutiveMatrix(Variables,ConstitutiveTensor);

      // VectorType StrainVector;
      // StrainVector = ConstitutiveModelUtilities::StrainTensorToVector(Variables.Strain.Matrix,StrainVector);
      
      // VectorType StressVector;
      // this->CalculateAndAddStressTensor(Variables,ConstitutiveTensor,StrainVector,StressVector);
      // rStressMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(StressVector,rStressMatrix);

      this->CalculateAndAddStressTensor(Variables,rStressMatrix);
      
      const StressMeasureType& rStressMeasure = rValues.GetStressMeasure();
      
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){
	
	const MatrixType& rTotalDeformationMatrix = rValues.GetTotalDeformationMatrix();
	
	//Variables.Strain.InverseMatrix used as an auxiliar matrix (contravariant push forward)
	noalias( Variables.Strain.InverseMatrix ) = prod( trans(rTotalDeformationMatrix), rStressMatrix );
	noalias( rStressMatrix )  = prod( Variables.Strain.InverseMatrix, rTotalDeformationMatrix );
	
      }     

      rConstitutiveMatrix = ConstitutiveModelUtilities::ConstitutiveTensorToMatrix(ConstitutiveTensor,rConstitutiveMatrix);

      KRATOS_CATCH(" ")
    }

    
    /**
     * Check
     */    
    virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
	KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value" << std::endl;

      if(POISSON_RATIO.Key() == 0){
	KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value" << std::endl;
      }
      else{
	const double& nu = rMaterialProperties[POISSON_RATIO];
	if( (nu > 0.499 && nu < 0.501) || (nu < -0.999 && nu > -1.01) )
	  KRATOS_ERROR << "POISSON_RATIO has an invalid value" << std::endl;
      }
      	
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

    void CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables) override
    {
      KRATOS_TRY

      //set model data pointer
      rVariables.SetModelData(rValues);
      rVariables.SetState(rValues.State);
    
      //deformation gradient
      const MatrixType& rDeltaDeformationMatrix  = rValues.GetDeltaDeformationMatrix();
      const MatrixType& rTotalDeformationMatrix  = rValues.GetTotalDeformationMatrix();

      const StressMeasureType& rStressMeasure  = rValues.GetStressMeasure();
      
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Strain.Matrix = GreenLagrangeTensor

	//set working strain measure
	rValues.SetStrainMeasure(ConstitutiveModelData::CauchyGreen_Right);
	
	//historical strain matrix
	rValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(this->mHistoryVector,rValues.StrainMatrix);

	//current strain matrix b
	noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rDeltaDeformationMatrix));
	noalias(rValues.StrainMatrix) = prod(rDeltaDeformationMatrix, rVariables.Strain.Matrix);
	
	//inverted total deformation gradient
	double detF = 0;
	ConstitutiveModelUtilities::InvertMatrix3( rTotalDeformationMatrix, rVariables.Strain.InverseMatrix, detF ); //InverseMatrix used as wildcard here (InverseMatrix = InverseTotalDeformationGradient)

	//strain measure C
	noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rVariables.Strain.InverseMatrix));
	rVariables.Strain.InverseMatrix = prod(trans(rTotalDeformationMatrix), rVariables.Strain.Matrix); //InverseMatrix used as a wildcard here (InverseMatrix = RightCauchyGreenMatrix )
      

	ConstitutiveModelUtilities::RightCauchyToGreenLagrangeStrain( rVariables.Strain.InverseMatrix, rVariables.Strain.Matrix);  
	
	rValues.State.Set(ConstitutiveModelData::STRAIN_COMPUTED);       
	
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Strain.Matrix = GreenLagrangeTensor
	
	//set working strain measure
	rValues.SetStrainMeasure(ConstitutiveModelData::CauchyGreen_Left);
	
	//historical strain matrix
	rValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(this->mHistoryVector,rValues.StrainMatrix);

	//current strain matrix b
	noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rDeltaDeformationMatrix));
	noalias(rValues.StrainMatrix) = prod(rDeltaDeformationMatrix, rVariables.Strain.Matrix);
	
	ConstitutiveModelUtilities::LeftCauchyToAlmansiStrain( rValues.StrainMatrix, rVariables.Strain.Matrix);
	
	//rVariables.Strain.InverseMatrix used as an auxiliar matrix (covariant pull back)  to GreenLagrangeStrain
	noalias( rVariables.Strain.InverseMatrix ) = prod( trans(rTotalDeformationMatrix), rVariables.Strain.Matrix );
	noalias( rVariables.Strain.Matrix )  = prod( rVariables.Strain.InverseMatrix, rTotalDeformationMatrix );

	//set as the current strain
	rValues.State.Set(ConstitutiveModelData::STRAIN_COMPUTED);


      }
      else{
	
	//set working strain measure
	rValues.SetStrainMeasure(ConstitutiveModelData::CauchyGreen_None);
	KRATOS_ERROR << "calling initialize SaintVenantKirchhoffModel .. StressMeasure is inconsistent"  << std::endl;
	
      }

      
      KRATOS_CATCH(" ")
    }


    void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, bounded_matrix<double,6,6>& rConstitutiveTensor, VectorType& rStrainVector, VectorType& rStressVector)
    {
      KRATOS_TRY

      noalias(rStressVector) = prod(rConstitutiveTensor,rStrainVector);
      
      rVariables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);
    
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
	
      rStressMatrix  = rVariables.Strain.Matrix;
      rStressMatrix *= 2.0 * rLameMu;
      
      double trace = (rVariables.Strain.Matrix(0,0)+rVariables.Strain.Matrix(1,1)+rVariables.Strain.Matrix(2,2));
      trace *= rLameLambda;

      rStressMatrix(0,0) += trace;
      rStressMatrix(1,1) += trace;
      rStressMatrix(2,2) += trace;
	    
      rVariables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);
    
      KRATOS_CATCH(" ")
    }

    
    void CalculateAndAddConstitutiveMatrix(HyperElasticDataType& rVariables, bounded_matrix<double,6,6>& rConstitutiveTensor)
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

    
      rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED);

    
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


