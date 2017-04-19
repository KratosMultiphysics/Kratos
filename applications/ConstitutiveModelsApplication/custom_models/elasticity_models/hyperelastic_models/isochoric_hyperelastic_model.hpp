//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_ISOCHORIC_HYPERELASTIC_MODEL_H_INCLUDED )
#define  KRATOS_ISOCHORIC_HYPERELASTIC_MODEL_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IsochoricHyperElasticModel : public HyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsochoricHyperElasticModel
    KRATOS_CLASS_POINTER_DEFINITION( IsochoricHyperElasticModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsochoricHyperElasticModel()
      {
      }
    
    /// Copy constructor.
    IsochoricHyperElasticModel(IsochoricHyperElasticModel const& rOther)
      {
      }

    /// Assignment operator.
    IsochoricHyperElasticModel& operator=(IsochoricHyperElasticModel const& rOther)
      {
	return *this;
      }

    /// Clone.
    virtual ElasticityModel::Pointer Clone() const override
    {
      return ( IsochoricHyperElasticModel::Pointer(new IsochoricHyperElasticModel(*this)) );      
    }
 
    /// Destructor.
    virtual ~IsochoricHyperElasticModel() {}


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

      rDensityFunction = 0;
      this->CalculateAndAddIsochoricStrainEnergy( Variables, rDensityFunction );
      this->CalculateAndAddVolumetricStrainEnergy( Variables, rDensityFunction );

	
      KRATOS_CATCH(" ")
    }

        
    virtual void CalculateStress(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_TRY
	
      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      this->CalculateAndAddIsochoricStress(Variables, rStressMatrix);
      
      rValues.StressMatrix = rStressMatrix; //store isochoric stress matrix as StressMatrix

      this->CalculateAndAddVolumetricStress(Variables, rStressMatrix);
            
      KRATOS_CATCH(" ")
    }


    virtual double& ConstitutiveComponent(ModelDataType& rValues, double &rCabcd,
					  const unsigned int& a, const unsigned int& b,
					  const unsigned int& c, const unsigned int& d)
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);
      this->CalculateScalingFactors(Variables);

      rValues.StressMatrix.clear();
      this->CalculateAndAddIsochoricStress(Variables, rValues.StressMatrix); //store isochoric stress in StressMatrix
      
      rCabcd = this->AddIsochoricConstitutiveComponent(Variables,rCabcd,a,b,c,d);
      rCabcd = this->AddVolumetricConstitutiveComponent(Variables,rCabcd,a,b,c,d);

      return rCabcd;
	
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
        buffer << "IsochoricHyperElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IsochoricHyperElasticModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IsochoricHyperElasticModel Data";
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

	
    virtual void CalculateAndAddIsochoricStress(HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      const ModelDataType&  rModelData        = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
	
      MatrixType StressPartMatrix;
      MatrixType StressMatrix;
      
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Variables.Strain.CauchyGreenMatrix = RightCauchyGreen (C)

	StressPartMatrix = GetI1RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(rStressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

	StressPartMatrix = GetI2RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(rStressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

	StressPartMatrix = GetI3RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(rStressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

	StressPartMatrix = GetIsochoricRightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	StressMatrix = prod(StressMatrix, StressPartMatrix);

	StressMatrix *= 2.0;
	
	rStressMatrix += StressMatrix;
	
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Variables.Strain.CauchyGreenMatrix = LeftCauchyGreen (b)

	StressPartMatrix = GetI1LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

	StressPartMatrix = GetI2LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

	StressPartMatrix = GetI3LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

	StressPartMatrix = GetIsochoricLeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	StressMatrix = prod(StressMatrix, StressPartMatrix);

	StressMatrix *= 2.0 * rVariables.Strain.Invariants.J;

	rStressMatrix += StressMatrix;
      }


      
      KRATOS_CATCH(" ")
    }


    virtual void CalculateAndAddVolumetricStress(HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      const ModelDataType&  rModelData        = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
	
      MatrixType StressMatrix;
            
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Variables.Strain.CauchyGreenMatrix = RightCauchyGreen (C)

	StressMatrix  = GetJRightCauchyGreenDerivative(rVariables.Strain,StressMatrix);
	StressMatrix *= rVariables.Factors.Alpha4;

	StressMatrix *= 2.0;

	rStressMatrix += StressMatrix;
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Variables.Strain.CauchyGreenMatrix = LeftCauchyGreen (b)

	StressMatrix  = GetJLeftCauchyGreenDerivative(rVariables.Strain,StressMatrix);
	StressMatrix *= rVariables.Factors.Alpha4;
	
	StressMatrix *= 2.0 * rVariables.Strain.Invariants.J;

	rStressMatrix += StressMatrix;
      }
      	
	
      KRATOS_CATCH(" ")
    }


    virtual double& AddIsochoricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						      const unsigned int& a, const unsigned int& b,
						      const unsigned int& c, const unsigned int& d)
    {
      KRATOS_TRY

      const ModelDataType& rModelData         = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
	
      double Cabcd = 0;
      double Cabef = 0;
      double Ccdef = 0;
      double Cefmn = 0;
      
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Variables.Strain.CauchyGreenMatrix = RightCauchyGreen (C)

	
	for(unsigned int e=0; e<3; e++)
	  {
	    for(unsigned int f=0; f<3; f++)
	      {

		Cabef = GetIsochoricRightCauchyGreenDerivative(rVariables.Strain,Cabef,a,b,e,f);
		Ccdef = GetIsochoricRightCauchyGreenDerivative(rVariables.Strain,Ccdef,c,d,e,f);

		rCabcd -= ( rVariables.Strain.InverseCauchyGreenMatrix(c,d) * Cabef + rVariables.Strain.InverseCauchyGreenMatrix(a,b) * Ccdef + rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13 * ( (rVariables.Strain.InverseCauchyGreenMatrix(a,b)*rVariables.Strain.InverseCauchyGreenMatrix(d,c)/3) - rVariables.Strain.InverseCauchyGreenMatrix(a,c)*rVariables.Strain.InverseCauchyGreenMatrix(b,d) ) * rVariables.Strain.CauchyGreenMatrix(e,f) ) * (2.0/3.0) * rModelData.GetStressMatrix()(e,f);

		for(unsigned int m=0; m<3; m++)
		  {
		    for(unsigned int n=0; n<3; n++)
		      {
			//2nd derivatives
			Cabcd = GetI1RightCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Alpha1 * Cabcd;
			
			Cabcd = GetI2RightCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Alpha2 * Cabcd;
			
			Cabcd = GetI3RightCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Alpha3 * Cabcd;
			
			//1st derivatives
			Cabcd = GetI1RightCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Beta1 * Cabcd;
			
			Cabcd = GetI2RightCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Beta2 * Cabcd;
			
			Cabcd = GetI3RightCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Beta3 * Cabcd;
			
			Cefmn *= 4.0;
		      }
		  }

		rCabcd += Cabef * Cefmn * Ccdef; 
	      }
	  }
		
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Variables.Strain.CauchyGreenMatrix = LeftCauchyGreen (b)
	
	for(unsigned int e=0; e<3; e++)
	  {
	    for(unsigned int f=0; f<3; f++)
	      {

		Cabef = GetIsochoricLeftCauchyGreenDerivative(rVariables.Strain,Cabef,a,b,e,f);
		Ccdef = GetIsochoricLeftCauchyGreenDerivative(rVariables.Strain,Ccdef,c,d,e,f);

		rCabcd -= ( msIdentityMatrix(c,d) * Cabef + msIdentityMatrix(a,b) * Ccdef + rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13 * ( (msIdentityMatrix(a,b)*msIdentityMatrix(d,c)/3) - msIdentityMatrix(a,c)*msIdentityMatrix(b,d) ) * rVariables.Strain.CauchyGreenMatrix(e,f) ) * (2.0/3.0) * rModelData.GetStressMatrix()(e,f);

		for(unsigned int m=0; m<3; m++)
		  {
		    for(unsigned int n=0; n<3; n++)
		      {
			//2nd derivatives
			Cabcd = GetI1LeftCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Alpha1 * Cabcd;
			
			Cabcd = GetI2LeftCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Alpha2 * Cabcd;
			
			Cabcd = GetI3LeftCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Alpha3 * Cabcd;
			
			//1st derivatives
			Cabcd = GetI1LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Beta1 * Cabcd;
			
			Cabcd = GetI2LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Beta2 * Cabcd;
			
			Cabcd = GetI3LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
			Cefmn += rVariables.Factors.Beta3 * Cabcd;
			
			Cefmn *= 4.0;
		      }
		  }

		rCabcd += Cabef * Cefmn * Ccdef; 
	      }
	  }
		
      }

      return rCabcd;
	
      KRATOS_CATCH(" ")
    }
    

    virtual double& AddVolumetricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						       const unsigned int& a, const unsigned int& b,
						       const unsigned int& c, const unsigned int& d)
    {
      KRATOS_TRY

      const ModelDataType& rModelData         = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();
	
      double Cabcd = 0;
      double nCabcd = 0;
      
      if( rStressMeasure == ConstitutiveModelData::StressMeasure_PK2 ){ //Variables.Strain.CauchyGreenMatrix = RightCauchyGreen (C)

	//2nd derivatives
	Cabcd = GetJRightCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,a,b,c,d);
	nCabcd += rVariables.Factors.Alpha4 * Cabcd;
	
	//1st derivatives
	Cabcd = GetJRightCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,a,b,c,d);
	nCabcd += rVariables.Factors.Beta4 * Cabcd;

	nCabcd *= 4.0;

	rCabcd += nCabcd;
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasure_Kirchhoff ){ //Variables.Strain.CauchyGreenMatrix = LeftCauchyGreen (b)
	//2nd derivatives
	Cabcd = GetJLeftCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,a,b,c,d);
	nCabcd += rVariables.Factors.Alpha4 * Cabcd;
	
	//1st derivatives
	Cabcd = GetJLeftCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,a,b,c,d);
	nCabcd += rVariables.Factors.Beta4 * Cabcd;

	nCabcd *= 4.0 * rVariables.Strain.Invariants.J;
	
	rCabcd += nCabcd;
      }	
	
      return rCabcd;
      
      KRATOS_CATCH(" ")
    }
    
    // set the default volumetric function for the incompressible case
    
    virtual void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction)
    {
      KRATOS_TRY

      const ModelDataType&  rValues = rVariables.GetModelData();
	
      rVolumetricDensityFunction += rValues.GetPressure() * rVariables.Strain.Invariants.J;
	
      KRATOS_CATCH(" ")
    }
    

    virtual void CalculateInvariants(HyperElasticDataType& rVariables)
    {
      KRATOS_TRY

      //isochoric-volumetric splitted law
      rVariables.Strain.CauchyGreenMatrix *=  rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13;
      
      //invariants
      this->CalculateStrainInvariants( rVariables.Strain.CauchyGreenMatrix, rVariables.Strain.Invariants.I1, rVariables.Strain.Invariants.I2, rVariables.Strain.Invariants.I3 );

      rVariables.Strain.Invariants.I3 = rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13; //for volumetric consistency in splitted law

      
      KRATOS_CATCH(" ")
    }

    
    virtual void CalculateScalingFactors(HyperElasticDataType& rVariables)
    {
      KRATOS_TRY

      HyperElasticModel::CalculateScalingFactors(rVariables);
	
      rVariables.Factors.Alpha4 = this->GetVolumetricFunctionJDerivative(rVariables,rVariables.Factors.Alpha4);
      rVariables.Factors.Beta4  = this->GetVolumetricFunction2ndJDerivative(rVariables,rVariables.Factors.Beta4);
		
      KRATOS_CATCH(" ")
    }

    // set the default volumetric function for the incompressible case
    
    virtual double& GetVolumetricFunctionJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //dU/dJ
    {
      KRATOS_TRY
	
      const ModelDataType&  rValues = rVariables.GetModelData();
	
      rDerivative = rValues.GetPressure();

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

  }; // Class IsochoricHyperElasticModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ISOCHORIC_HYPERELASTIC_MODEL_H_INCLUDED  defined 


