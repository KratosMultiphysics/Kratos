//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_ISOCHORIC_MOONEY_RIVLIN_MODEL_H_INCLUDED )
#define  KRATOS_ISOCHORIC_MOONEY_RIVLIN_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/mooney_rivlin_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IsochoricMooneyRivlinModel : public MooneyRivlinModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsochoricMooneyRivlinModel
    KRATOS_CLASS_POINTER_DEFINITION( IsochoricMooneyRivlinModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsochoricMooneyRivlinModel() : MooneyRivlinModel() {}

    /// Copy constructor.
    IsochoricMooneyRivlinModel(IsochoricMooneyRivlinModel const& rOther) : MooneyRivlinModel(rOther) {}

    /// Assignment operator.
    IsochoricMooneyRivlinModel& operator=(IsochoricMooneyRivlinModel const& rOther)
    {
	MooneyRivlinModel::operator=(rOther);
	return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<IsochoricMooneyRivlinModel>(*this);
    }

    /// Destructor.
    ~IsochoricMooneyRivlinModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction) override
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues, Variables);

      rDensityFunction = 0;
      this->CalculateAndAddIsochoricStrainEnergy( Variables, rDensityFunction );
      this->CalculateAndAddVolumetricStrainEnergy( Variables, rDensityFunction );


      KRATOS_CATCH(" ")
    }


    void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      this->CalculateAndAddIsochoricStressTensor(Variables, rStressMatrix);

      rValues.StressMatrix = rStressMatrix; //store isochoric stress matrix as StressMatrix

      this->CalculateAndAddVolumetricStressTensor(Variables, rStressMatrix);

      Variables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

      KRATOS_CATCH(" ")
    }


    void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY

      HyperElasticDataType Variables;
      this->CalculateStrainData(rValues,Variables);

      //Calculate Stress Matrix
      this->CalculateAndAddIsochoricStressTensor(Variables, rStressMatrix);

      rValues.StressMatrix = rStressMatrix; //store isochoric stress matrix as StressMatrix

      this->CalculateAndAddVolumetricStressTensor(Variables, rStressMatrix);

      //Calculate Constitutive Matrix
      this->CalculateAndAddConstitutiveTensor(Variables,rConstitutiveMatrix);

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IsochoricMooneyRivlinModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IsochoricMooneyRivlinModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IsochoricMooneyRivlinModel Data";
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


    void CalculateAndAddIsochoricStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      const ModelDataType&  rModelData        = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

      MatrixType StressPartMatrix;
      MatrixType StressMatrix;

      if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //Variables.Strain.Matrix = RightCauchyGreen (C)

	StressPartMatrix = GetI1RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

	StressPartMatrix = GetI2RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

	StressPartMatrix = GetI3RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

	//option a:
	StressPartMatrix = GetIsochoricRightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	StressMatrix = prod(StressMatrix, StressPartMatrix);

	//option b:
	// StressPartMatrix = StressMatrix;
	// StressMatrix.clear();
	// double Derivative = 0;
	// for (unsigned int i = 0; i < 3; i++) {
	//   for (unsigned int j = 0; j < 3; j++) {
	//     for (unsigned int k = 0; k < 3; k++) {
	//       for (unsigned int l = 0; l < 3; l++) {
	// 	Derivative = GetIsochoricRightCauchyGreenDerivative( rVariables.Strain, Derivative, i, j, k, l);
	// 	StressMatrix(i,j) += Derivative * StressPartMatrix(k, l);
	//       }
	//     }
	//   }
	// }

	StressMatrix *= 2.0;

	rStressMatrix += StressMatrix;

      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //Variables.Strain.Matrix = LeftCauchyGreen (b)

         StressPartMatrix = GetI1LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
         noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

         StressPartMatrix = GetI2LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
         noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

         StressPartMatrix = GetI3LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
         noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

	 //option a:
         StressPartMatrix = GetIsochoricLeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
         StressMatrix = prod(StressMatrix, StressPartMatrix);

	 //option b:
	 // StressPartMatrix = StressMatrix;
	 // StressMatrix.clear();
	 // double Derivative = 0;
	 // for (unsigned int i = 0; i < 3; i++) {
	 //   for (unsigned int j = 0; j < 3; j++) {
	 //     for (unsigned int k = 0; k < 3; k++) {
	 //       for (unsigned int l = 0; l < 3; l++) {
	 // 	Derivative = GetIsochoricLeftCauchyGreenDerivative( rVariables.Strain, Derivative, i, j, k, l);
	 // 	StressMatrix(i,j) += Derivative * StressPartMatrix(k, l);
	 //       }
	 //     }
	 //   }
	 // }

         StressMatrix *= 2.0;

         rStressMatrix += StressMatrix;
      }


      KRATOS_CATCH(" ")
    }


    void CalculateAndAddVolumetricStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      const ModelDataType&  rModelData        = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

      MatrixType StressMatrix;

      if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //Variables.Strain.Matrix = RightCauchyGreen (C)

	StressMatrix  = GetJRightCauchyGreenDerivative(rVariables.Strain,StressMatrix);
	StressMatrix *= rVariables.Factors.Alpha4;

	StressMatrix *= 2.0;

	noalias(rStressMatrix) += StressMatrix;
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //Variables.Strain.Matrix = LeftCauchyGreen (b)

	StressMatrix  = GetJLeftCauchyGreenDerivative(rVariables.Strain,StressMatrix);
	StressMatrix *= rVariables.Factors.Alpha4;

	StressMatrix *= 2.0;

	noalias(rStressMatrix) += StressMatrix;
      }


      KRATOS_CATCH(" ")
    }



    double& AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
					     const unsigned int& a, const unsigned int& b,
					     const unsigned int& c, const unsigned int& d) override
    {
      KRATOS_TRY

      rCabcd = this->AddIsochoricConstitutiveComponent(rVariables,rCabcd,a,b,c,d);

      rCabcd = this->AddVolumetricConstitutiveComponent(rVariables,rCabcd,a,b,c,d);

      //std::cout<<" Cabcd ["<<a<<","<<b<<","<<c<<","<<d<<"] :"<<rCabcd<<std::endl;

      return rCabcd;

      KRATOS_CATCH(" ")
    }




    double& AddIsochoricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						      const unsigned int& a, const unsigned int& b,
						      const unsigned int& c, const unsigned int& d) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData         = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

      double Dabcd = 0;
      double Cabef = 0;
      double Ccdef = 0;
      //double Cabcd = 0;
      //double Cefmn = 0;
      //double Ccdmn = 0;

      MatrixType StressMatrix;
      MatrixType StressPartMatrix;


      if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //Variables.Strain.CauchyGreenMatrix = RightCauchyGreen (C)

	//constitutive tensor for PK2 still not giving the correct numbers in the 4 and 5 row, something is wrong

	StressPartMatrix = GetI1RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

	StressPartMatrix = GetI2RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

	StressPartMatrix = GetI3RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

	StressMatrix *= 2.0;

	Dabcd = ( rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13 * ( (rVariables.Strain.InverseMatrix(a,b)*rVariables.Strain.InverseMatrix(d,c)/3.0) - 0.5 * (rVariables.Strain.InverseMatrix(a,c)*rVariables.Strain.InverseMatrix(b,d)+rVariables.Strain.InverseMatrix(a,d)*rVariables.Strain.InverseMatrix(b,c)) ) );


	for(unsigned int e=0; e<3; e++)
	  {
	    for(unsigned int f=0; f<3; f++)
	      {

		Cabef = GetIsochoricRightCauchyGreenDerivative(rVariables.Strain,Cabef,a,b,e,f);
		Ccdef = GetIsochoricRightCauchyGreenDerivative(rVariables.Strain,Ccdef,c,d,e,f);

		rCabcd += ( rVariables.Strain.InverseMatrix(c,d) * Cabef + rVariables.Strain.InverseMatrix(a,b) * Ccdef + Dabcd * rVariables.Strain.Matrix(e,f) ) * (-2.0/3.0) * StressMatrix(e,f);

		// splitted option:
		//rCabcd += ( rVariables.Strain.InverseMatrix(c,d) * Cabef + rVariables.Strain.InverseMatrix(a,b) * Ccdef ) * (-2.0/3.0) * StressMatrix(e,f);
		//rCabcd += ( Dabcd * rVariables.Strain.Matrix(e,f) ) * (-2.0/3.0) * StressMatrix(e,f);

		//some terms are not included in the ddW/dCdC: (generally it is not needed)
		// for(unsigned int m=0; m<3; m++)
		//   {
		//     for(unsigned int n=0; n<3; n++)
		//       {
                //         Ccdmn = GetIsochoricRightCauchyGreenDerivative( rVariables.Strain, Ccdmn,c,d,m,n);

		// 	//2nd derivatives
		// 	Cabcd = GetI1RightCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Alpha1 * Cabcd * Cabef * Ccdmn;

		// 	Cabcd = GetI2RightCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Alpha2 * Cabcd * Cabef * Ccdmn;

		// 	Cabcd = GetI3RightCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Alpha3 * Cabcd * Cabef * Ccdmn;

		// 	//1st derivatives
		// 	Cabcd = GetI1RightCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Beta1 * Cabcd * Cabef * Ccdmn;

		// 	Cabcd = GetI2RightCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Beta2 * Cabcd * Cabef * Ccdmn;

		// 	Cabcd = GetI3RightCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Beta3 * Cabcd * Cabef * Ccdmn;

		//       }
		//   }

	      }
	  }

	//rCabcd += Cefmn * 4;

      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //Variables.Strain.Matrix = LeftCauchyGreen (b)

	//constitutive tensor for kirchhoff still not giving the correct numbers, something is wrong

	StressPartMatrix = GetI1LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

	StressPartMatrix = GetI2LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

	StressPartMatrix = GetI3LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
	noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

	StressMatrix *= 2.0;

	Dabcd = ( rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13 * ( (msIdentityMatrix(a,b)*msIdentityMatrix(d,c)/3.0) - 0.5 * (msIdentityMatrix(a,c)*msIdentityMatrix(b,d)+msIdentityMatrix(a,d)*msIdentityMatrix(b,c)) ) );

	for(unsigned int e=0; e<3; e++)
	  {
	    for(unsigned int f=0; f<3; f++)
	      {

		Cabef = GetIsochoricLeftCauchyGreenDerivative(rVariables.Strain,Cabef,a,b,e,f);
		Ccdef = GetIsochoricLeftCauchyGreenDerivative(rVariables.Strain,Ccdef,c,d,e,f);

		rCabcd += ( msIdentityMatrix(c,d) * Cabef + msIdentityMatrix(a,b) * Ccdef + Dabcd * msIdentityMatrix(e,f) ) * (-2.0/3.0) * StressMatrix(e,f);

		// splitted option:
		//rCabcd += ( msIdentityMatrix(c,d) * Cabef + msIdentityMatrix(a,b) * Ccdef ) * (-2.0/3.0) * StressMatrix(e,f);
		//rCabcd += ( Dabcd * msIdentityMatrix(e,f) ) * (-2.0/3.0) * StressMatrix(e,f);


		// some terms are not included in the ddW/dbdb: (generally it is not needed)
		// for(unsigned int m=0; m<3; m++)
		//   {
		//     for(unsigned int n=0; n<3; n++)
		//       {
                //         Ccdmn = GetIsochoricLeftCauchyGreenDerivative( rVariables.Strain, Ccdmn,c,d,m,n);

		// 	//2nd derivatives
		//      //check why this term is not needed
		// 	//Cabcd = GetI1LeftCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	//Cefmn += rVariables.Factors.Alpha1 * Cabcd * Cabef * Ccdmn;

		// 	Cabcd = GetI2LeftCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Alpha2 * Cabcd * Cabef * Ccdmn;

		// 	Cabcd = GetI3LeftCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Alpha3 * Cabcd * Cabef * Ccdmn;

		// 	//1st derivatives
		// 	Cabcd = GetI1LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Beta1 * Cabcd * Cabef * Ccdmn;

		// 	Cabcd = GetI2LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Beta2 * Cabcd * Cabef * Ccdmn;

		// 	Cabcd = GetI3LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,e,f,m,n);
		// 	Cefmn += rVariables.Factors.Beta3 * Cabcd * Cabef * Ccdmn;

		//       }
		//   }

	      }
	  }

	//rCabcd += Cefmn * 4;

      }

      return rCabcd;

      KRATOS_CATCH(" ")
    }


    double& AddVolumetricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						       const unsigned int& a, const unsigned int& b,
						       const unsigned int& c, const unsigned int& d) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData         = rVariables.GetModelData();
      const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

      double Cabcd = 0;
      double nCabcd = 0;

      if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //Variables.Strain.Matrix = RightCauchyGreen (C)

	//2nd derivatives
	Cabcd = GetJRightCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,a,b,c,d);
	nCabcd += rVariables.Factors.Alpha4 * Cabcd;

	//1st derivatives
	Cabcd = GetJRightCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,a,b,c,d);
	nCabcd += rVariables.Factors.Beta4 * Cabcd;

	nCabcd *= 4.0;

	rCabcd += nCabcd;
      }
      else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //Variables.Strain.Matrix = LeftCauchyGreen (b)
	//2nd derivatives
	Cabcd = GetJLeftCauchyGreen2ndDerivative(rVariables.Strain,Cabcd,a,b,c,d);
	nCabcd += rVariables.Factors.Alpha4 * Cabcd;

	//1st derivatives
	Cabcd = GetJLeftCauchyGreenSquare1stDerivative(rVariables.Strain,Cabcd,a,b,c,d);
	nCabcd += rVariables.Factors.Beta4 * Cabcd;

	nCabcd *= 4.0;

	rCabcd += nCabcd;
      }

      return rCabcd;

      KRATOS_CATCH(" ")
    }

    // set the default volumetric function for the incompressible case

    void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction) override
    {
      KRATOS_TRY

      const ModelDataType&  rValues = rVariables.GetModelData();

      rVolumetricDensityFunction += rValues.GetPressure() * rVariables.Strain.Invariants.J;

      KRATOS_CATCH(" ")
    }


    void CalculateScalingFactors(HyperElasticDataType& rVariables) override
    {
      KRATOS_TRY

      MooneyRivlinModel::CalculateScalingFactors(rVariables);

      rVariables.Factors.Alpha4 = this->GetVolumetricFunction1stJDerivative(rVariables,rVariables.Factors.Alpha4);
      rVariables.Factors.Beta4  = this->GetVolumetricFunction2ndJDerivative(rVariables,rVariables.Factors.Beta4);

      KRATOS_CATCH(" ")
    }

    // set the default volumetric function for the incompressible case

    double& GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //dU/dJ
    {
      KRATOS_TRY

      const ModelDataType&  rValues = rVariables.GetModelData();

      rDerivative = rValues.GetPressure();

      return rDerivative;

      KRATOS_CATCH(" ")
    };


    double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddU/dJdJ
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


    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MooneyRivlinModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MooneyRivlinModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IsochoricMooneyRivlinModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ISOCHORIC_MOONEY_RIVLIN_MODEL_H_INCLUDED  defined
