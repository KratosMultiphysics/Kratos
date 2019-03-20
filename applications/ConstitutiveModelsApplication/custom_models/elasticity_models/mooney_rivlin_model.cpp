//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/mooney_rivlin_model.hpp"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  MooneyRivlinModel::MooneyRivlinModel()
    : HyperElasticModel()
  {
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  MooneyRivlinModel::MooneyRivlinModel(const MooneyRivlinModel& rOther)
    : HyperElasticModel(rOther)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveModel::Pointer MooneyRivlinModel::Clone() const
  {
    return Kratos::make_shared<MooneyRivlinModel>(*this);
  }

  //********************************ASSIGNMENT******************************************
  //************************************************************************************
  MooneyRivlinModel& MooneyRivlinModel::operator=(MooneyRivlinModel const& rOther)
  {
    HyperElasticModel::operator=(rOther);
    return *this;
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  MooneyRivlinModel::~MooneyRivlinModel()
  {
  }


  //************************************************************************************
  //************************************************************************************

  void MooneyRivlinModel::CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

    MatrixType StressPartMatrix;
    MatrixType StressMatrix;

    if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //Strain.Matrix = RightCauchyGreen (C)

      StressPartMatrix = GetI1RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

      StressPartMatrix = GetI2RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

      StressPartMatrix = GetI3RightCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

      StressMatrix *= 2.0;

      rStressMatrix += StressMatrix;

    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //Strain.Matrix = LeftCauchyGreen (b)

      StressPartMatrix = GetI1LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix)  = rVariables.Factors.Alpha1 * StressPartMatrix;

      StressPartMatrix = GetI2LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix) += rVariables.Factors.Alpha2 * StressPartMatrix;

      StressPartMatrix = GetI3LeftCauchyGreenDerivative(rVariables.Strain,StressPartMatrix);
      noalias(StressMatrix) += rVariables.Factors.Alpha3 * StressPartMatrix;

      StressMatrix *= 2.0;

      rStressMatrix += StressMatrix;
    }


    rVariables.State().Set(ConstitutiveModelData::STRESS_COMPUTED);

    KRATOS_CATCH(" ")
  }


  //***********************PROTECTED OPERATIONS FROM BASE CLASS*************************
  //************************************************************************************

  void MooneyRivlinModel::CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables)
  {
    KRATOS_TRY

    //set model data pointer
    rVariables.SetModelData(rValues);
    rVariables.SetState(rValues.State);

    //deformation gradient
    const MatrixType& rDeltaDeformationMatrix = rValues.GetDeltaDeformationMatrix();
    const MatrixType& rTotalDeformationMatrix = rValues.GetTotalDeformationMatrix();

    const StressMeasureType& rStressMeasure = rValues.GetStressMeasure();

    if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //Strain.Matrix = RightCauchyGreen (C=FT*F)  C^-1=(FT*F)^-1=F^-1*FT^-1

      //set working strain measure
      rValues.SetStrainMeasure(ConstitutiveModelData::StrainMeasureType::CauchyGreen_Right);

      //historical strain matrix
      rValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(this->mHistoryVector,rValues.StrainMatrix);

      //current strain matrix b
      noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rDeltaDeformationMatrix));
      noalias(rValues.StrainMatrix) = prod(rDeltaDeformationMatrix, rVariables.Strain.Matrix);

      //inverted total deformation gradient
      ConstitutiveModelUtilities::InvertMatrix3( rTotalDeformationMatrix, rVariables.Strain.InverseMatrix, rVariables.Strain.Invariants.I3 ); //InverseMatrix and I3 is used as wildcard here (InverseMatrix = InverseTotalDeformationGradient)

      //strain measure C
      noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rVariables.Strain.InverseMatrix));
      rVariables.Strain.Matrix = prod(trans(rTotalDeformationMatrix), rVariables.Strain.Matrix);


      //inverted strain measure
      ConstitutiveModelUtilities::InvertMatrix3( rVariables.Strain.Matrix, rVariables.Strain.InverseMatrix, rVariables.Strain.Invariants.I3 );

      rValues.State.Set(ConstitutiveModelData::STRAIN_COMPUTED);


    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //Strain.Matrix = LeftCauchyGreen (b=F*FT)

      //set working strain measure
      rValues.SetStrainMeasure(ConstitutiveModelData::StrainMeasureType::CauchyGreen_Left);

      //historical strain matrix
      rValues.StrainMatrix = ConstitutiveModelUtilities::VectorToSymmetricTensor(this->mHistoryVector,rValues.StrainMatrix);

      //current strain matrix b
      noalias(rVariables.Strain.Matrix) = prod(rValues.StrainMatrix,trans(rDeltaDeformationMatrix));
      noalias(rValues.StrainMatrix) = prod(rDeltaDeformationMatrix, rVariables.Strain.Matrix);


      noalias(rVariables.Strain.Matrix) = rValues.StrainMatrix;

      //inverted strain measure
      ConstitutiveModelUtilities::InvertMatrix3( rVariables.Strain.Matrix, rVariables.Strain.InverseMatrix, rVariables.Strain.Invariants.I3 );
      rValues.State.Set(ConstitutiveModelData::STRAIN_COMPUTED);

    }
    else{

      //set working strain measure
      rValues.SetStrainMeasure(ConstitutiveModelData::StrainMeasureType::CauchyGreen_None);
      KRATOS_ERROR << "calling initialize MooneyRivlinModel .. StressMeasure is inconsistent"  << std::endl;

    }

    //Calculate Invariants
    this->CalculateInvariants(rVariables);

    //Calculate LameMuBar
    rValues.MaterialParameters.LameMuBar = rValues.MaterialParameters.LameMu * ( rVariables.Strain.Matrix(0,0) + rVariables.Strain.Matrix(1,1) + rVariables.Strain.Matrix(2,2) ) * rVariables.Strain.Invariants.J_13 * rVariables.Strain.Invariants.J_13 * (1.0/3.0) ;

    //std::cout<<" LameMuBar "<<rValues.MaterialParameters.LameMuBar<<std::endl;

    //Algorithmic moduli factors
    this->CalculateScalingFactors(rVariables);


    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  double& MooneyRivlinModel::AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						      const unsigned int& a, const unsigned int& b,
						      const unsigned int& c, const unsigned int& d)
  {
    KRATOS_TRY

    //this is a simplified version, more terms of the derivatives are needed to be general:

    const ModelDataType&  rModelData        = rVariables.GetModelData();
    const StressMeasureType& rStressMeasure = rModelData.GetStressMeasure();

    double Cabcd = 0;
    double Dabcd = 0;

    if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_PK2 ){ //Strain.Matrix = LeftCauchyGreen (C)

      //2nd derivatives
      Dabcd = GetI1RightCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha1 * Dabcd;

      Dabcd = GetI2RightCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha2 * Dabcd;

      Dabcd = GetI3RightCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha3 * Dabcd;

      //1st derivatives
      Dabcd = GetI1RightCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta1 * Dabcd;

      Dabcd = GetI2RightCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta2 * Dabcd;

      Dabcd = GetI3RightCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta3 * Dabcd;

      Cabcd *= 4.0;


    }
    else if( rStressMeasure == ConstitutiveModelData::StressMeasureType::StressMeasure_Kirchhoff ){ //Strain.Matrix = LeftCauchyGreen (b)

      //2nd derivatives

      // check why this term is not needed
      // Dabcd = GetI1LeftCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      // Cabcd += rVariables.Factors.Alpha1 * Dabcd;

      Dabcd = GetI2LeftCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha2 * Dabcd;

      Dabcd = GetI3LeftCauchyGreen2ndDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Alpha3 * Dabcd;

      //1st derivatives
      Dabcd = GetI1LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta1 * Dabcd;

      Dabcd = GetI2LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta2 * Dabcd;

      Dabcd = GetI3LeftCauchyGreenSquare1stDerivative(rVariables.Strain,Dabcd,a,b,c,d);
      Cabcd += rVariables.Factors.Beta3 * Dabcd;


      Cabcd *= 4.0;

    }

    rCabcd += Cabcd;

    return rCabcd;


    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void MooneyRivlinModel::CalculateScalingFactors(HyperElasticDataType& rVariables)
  {
    KRATOS_TRY

    rVariables.Factors.Alpha1 = this->GetFunction1stI1Derivative(rVariables, rVariables.Factors.Alpha1);
    rVariables.Factors.Alpha2 = this->GetFunction1stI2Derivative(rVariables, rVariables.Factors.Alpha2);
    rVariables.Factors.Alpha3 = this->GetFunction1stI3Derivative(rVariables, rVariables.Factors.Alpha3);

    rVariables.Factors.Beta1  = this->GetFunction2ndI1Derivative(rVariables, rVariables.Factors.Beta1);
    rVariables.Factors.Beta2  = this->GetFunction2ndI2Derivative(rVariables, rVariables.Factors.Beta2);
    rVariables.Factors.Beta3  = this->GetFunction2ndI3Derivative(rVariables, rVariables.Factors.Beta3);

    // std::cout<<" Alpha["<<rVariables.Factors.Alpha1<<" "<<rVariables.Factors.Alpha2<<" "<<rVariables.Factors.Alpha3<<"]"<<std::endl;
    // std::cout<<" Beta ["<<rVariables.Factors.Beta1<<" "<<rVariables.Factors.Beta2<<" "<<rVariables.Factors.Beta3<<"]"<<std::endl;

    KRATOS_CATCH(" ")
  }


  //************// dW

  double& MooneyRivlinModel::GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI1
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in MooneyRivlinModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& MooneyRivlinModel::GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI2
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in MooneyRivlinModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //dW/dI3
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in MooneyRivlinModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************


  double& MooneyRivlinModel::GetFunction2ndI1Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI1dI1
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in MooneyRivlinModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetFunction2ndI2Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI2dI2
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in MooneyRivlinModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative) //ddW/dI3dI3
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class function in MooneyRivlinModel ... illegal operation" << std::endl;

    return rDerivative;

    KRATOS_CATCH(" ")
  }



  //************************************************************************************
  //************************************************************************************

  //isochoric volumetric split

  MooneyRivlinModel::MatrixType& MooneyRivlinModel::GetIsochoricRightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //dC'/dC
  {
    KRATOS_TRY

    noalias(rDerivative)  = rStrain.InverseMatrix;
    rDerivative *= -rStrain.Invariants.I1/3.0;
    noalias(rDerivative) += this->msIdentityMatrix;

    rDerivative *= rStrain.Invariants.J_13 * rStrain.Invariants.J_13;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  double& MooneyRivlinModel::GetIsochoricRightCauchyGreenDerivative(const StrainData& rStrain,
								    double& rDerivative,
								    const double& a,
								    const double& b,
								    const double& c,
								    const double& d) //dC'/dC
  {
    KRATOS_TRY

    rDerivative  = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative -= rStrain.InverseMatrix(a,b) * rStrain.Matrix(c,d)/3.0;

    rDerivative *= rStrain.Invariants.J_13 * rStrain.Invariants.J_13;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  MooneyRivlinModel::MatrixType& MooneyRivlinModel::GetIsochoricLeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //db'/db
  {
    KRATOS_TRY

    noalias(rDerivative)  = rStrain.InverseMatrix;
    rDerivative *= -rStrain.Invariants.I1/3.0;
    noalias(rDerivative) += this->msIdentityMatrix;

    rDerivative *= rStrain.Invariants.J_13 * rStrain.Invariants.J_13;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetIsochoricLeftCauchyGreenDerivative(const StrainData& rStrain,
								   double& rDerivative,
								   const double& a,
								   const double& b,
								   const double& c,
								   const double& d) //db'/db
  {
    KRATOS_TRY

    rDerivative  = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative -= this->msIdentityMatrix(a,b) * this->msIdentityMatrix(c,d)/3.0;

    rDerivative *= rStrain.Invariants.J_13 * rStrain.Invariants.J_13;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************


  //************// right cauchy green: C

  MooneyRivlinModel::MatrixType& MooneyRivlinModel::GetI1RightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //dI1/dC
  {
    KRATOS_TRY

    noalias(rDerivative) = this->msIdentityMatrix;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  MooneyRivlinModel::MatrixType& MooneyRivlinModel::GetI2RightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //dI2/dC
  {
    KRATOS_TRY

    noalias(rDerivative)  = this->msIdentityMatrix;
    rDerivative *= rStrain.Invariants.I1;
    noalias(rDerivative) += rStrain.Matrix;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  MooneyRivlinModel::MatrixType& MooneyRivlinModel::GetI3RightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //dI3/dC
  {
    KRATOS_TRY

    noalias(rDerivative)  = rStrain.InverseMatrix;
    rDerivative *= rStrain.Invariants.I3;

    return rDerivative;

    KRATOS_CATCH(" ")
  }



  double& MooneyRivlinModel::GetInverseRightCauchyGreenDerivative(const StrainData& rStrain,
								  double& rDerivative,
								  const double& a,
								  const double& b,
								  const double& c,
								  const double& d) //dC^-1/dC
  {
    KRATOS_TRY

    rDerivative = -0.5*(rStrain.InverseMatrix(a,c)*rStrain.InverseMatrix(b,d)+rStrain.InverseMatrix(a,d)*rStrain.InverseMatrix(b,c));

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //Invariants 1st derivatives by components
  double& MooneyRivlinModel::GetI1RightCauchyGreen1stDerivative(const StrainData& rStrain,
								double& rDerivative,
								const double& a,
								const double& b) //dI1/dC
  {
    KRATOS_TRY

    rDerivative = this->msIdentityMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI2RightCauchyGreen1stDerivative(const StrainData& rStrain,
								double& rDerivative,
								const double& a,
								const double& b) //dI2/dC
  {
    KRATOS_TRY

    rDerivative = rStrain.Invariants.I1 * this->msIdentityMatrix(a,b) + rStrain.Matrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI3RightCauchyGreen1stDerivative(const StrainData& rStrain,
								double& rDerivative,
								const double& a,
								const double& b) //dI3/dC
  {
    KRATOS_TRY

    rDerivative  = rStrain.Invariants.I3 * rStrain.InverseMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }



  //************************************************************************************
  //************************************************************************************

  //Invariants Square of the 1st derivatives by components
  double& MooneyRivlinModel::GetI1RightCauchyGreenSquare1stDerivative(const StrainData& rStrain,
								      double& rDerivative,
								      const double& a,
								      const double& b,
								      const double& c,
								      const double& d) //dI1/dC * dI1/dC
  {
    KRATOS_TRY

    rDerivative  = this->msIdentityMatrix(a,b);
    rDerivative *= this->msIdentityMatrix(c,d);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI2RightCauchyGreenSquare1stDerivative(const StrainData& rStrain,
								      double& rDerivative,
								      const double& a,
								      const double& b,
								      const double& c,
								      const double& d) //dI2/dC * dI2/dC
  {
    KRATOS_TRY

    rDerivative  = rStrain.Invariants.I1 * this->msIdentityMatrix(a,b) + rStrain.Matrix(a,b);
    rDerivative *= (rStrain.Invariants.I1 * this->msIdentityMatrix(c,d) + rStrain.Matrix(c,d));

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI3RightCauchyGreenSquare1stDerivative(const StrainData& rStrain,
								      double& rDerivative,
								      const double& a,
								      const double& b,
								      const double& c,
								      const double& d) //dI3/dC * dI3/dC
  {
    KRATOS_TRY

    rDerivative  = rStrain.Invariants.I3 * rStrain.InverseMatrix(a,b);
    rDerivative *= (rStrain.Invariants.I3 * rStrain.InverseMatrix(c,d));

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  //Invariants 2nd derivatives by components
  double& MooneyRivlinModel::GetI1RightCauchyGreen2ndDerivative(const StrainData& rStrain,
								double& rDerivative,
								const double& a,
								const double& b,
								const double& c,
								const double& d) //ddI1/dCdC
  {
    KRATOS_TRY

    rDerivative = 0.0;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI2RightCauchyGreen2ndDerivative(const StrainData& rStrain,
								double& rDerivative,
								const double& a,
								const double& b,
								const double& c,
								const double& d) //ddI2/dCdC
  {
    KRATOS_TRY

    rDerivative  = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative *= rStrain.Invariants.I1;
    rDerivative -= rStrain.Matrix(a,b)*this->msIdentityMatrix(c,d);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI3RightCauchyGreen2ndDerivative(const StrainData& rStrain,
								double& rDerivative,
								const double& a,
								const double& b,
								const double& c,
								const double& d) //ddI3/dCdC
  {
    KRATOS_TRY

    rDerivative  = GetInverseRightCauchyGreenDerivative(rStrain,rDerivative,a,b,c,d);
    rDerivative += rStrain.InverseMatrix(a,b)*rStrain.InverseMatrix(c,d);
    rDerivative *= rStrain.Invariants.I3;
    return rDerivative;

    KRATOS_CATCH(" ")
  }



  //************************************************************************************
  //************************************************************************************


  //************// left cauchy green : b

  MooneyRivlinModel::MatrixType& MooneyRivlinModel::GetI1LeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //dI1/db
  {
    KRATOS_TRY

    noalias(rDerivative) = rStrain.Matrix;

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  MooneyRivlinModel::MatrixType& MooneyRivlinModel::GetI2LeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative)  //dI2/db
  {
    KRATOS_TRY

    noalias(rDerivative)  = rStrain.Matrix;
    rDerivative *= rStrain.Invariants.I1;
    noalias(rDerivative) -= prod(rStrain.Matrix, rStrain.Matrix);


    return rDerivative;

    KRATOS_CATCH(" ")
  }

  MooneyRivlinModel::MatrixType& MooneyRivlinModel::GetI3LeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative) //dI3/db
  {
    KRATOS_TRY

    noalias(rDerivative) = this->msIdentityMatrix;
    rDerivative *= rStrain.Invariants.I3;

    return rDerivative;

    KRATOS_CATCH(" ")
  }


  //Invariants 1st derivatives by components
  double& MooneyRivlinModel::GetI1LeftCauchyGreen1stDerivative(const StrainData& rStrain,
							       double& rDerivative,
							       const double& a,
							       const double& b) //dI1/db
  {
    KRATOS_TRY

    rDerivative = rStrain.Matrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI2LeftCauchyGreen1stDerivative(const StrainData& rStrain,
							       double& rDerivative,
							       const double& a,
							       const double& b) //dI2/db
  {
    KRATOS_TRY

    rDerivative  = rStrain.Invariants.I1 * rStrain.Matrix(a,b) + rStrain.Matrix(a,1)*rStrain.Matrix(1,b) + rStrain.Matrix(a,2)*rStrain.Matrix(2,b) + rStrain.Matrix(a,3)*rStrain.Matrix(3,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI3LeftCauchyGreen1stDerivative(const StrainData& rStrain,
							       double& rDerivative,
							       const double& a,
							       const double& b) //dI3/db
  {
    KRATOS_TRY

    rDerivative  = rStrain.Invariants.I3 * this->msIdentityMatrix(a,b);

    return rDerivative;

    KRATOS_CATCH(" ")
  }



  //Invariants Square of the 1st derivatives by components
  double& MooneyRivlinModel::GetI1LeftCauchyGreenSquare1stDerivative(const StrainData& rStrain,
								     double& rDerivative,
								     const double& a,
								     const double& b,
								     const double& c,
								     const double& d) //dI1/db * dI1/db
  {
    KRATOS_TRY

    rDerivative  = rStrain.Matrix(a,b);
    rDerivative *= rStrain.Matrix(c,d);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI2LeftCauchyGreenSquare1stDerivative(const StrainData& rStrain,
								     double& rDerivative,
								     const double& a,
								     const double& b,
								     const double& c,
								     const double& d) //dI2/db * dI2/db
  {
    KRATOS_TRY

    rDerivative  = rStrain.Invariants.I1 * rStrain.Matrix(a,b) + rStrain.Matrix(a,0)*rStrain.Matrix(0,b) + rStrain.Matrix(a,1)*rStrain.Matrix(1,b) + rStrain.Matrix(a,2)*rStrain.Matrix(2,b);
    rDerivative *= (rStrain.Invariants.I1 * rStrain.Matrix(c,d) + rStrain.Matrix(c,0)*rStrain.Matrix(0,d) + rStrain.Matrix(c,1)*rStrain.Matrix(1,d) + rStrain.Matrix(c,2)*rStrain.Matrix(2,d));

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI3LeftCauchyGreenSquare1stDerivative(const StrainData& rStrain,
								     double& rDerivative,
								     const double& a,
								     const double& b,
								     const double& c,
								     const double& d) //dI3/db * dI3/db
  {
    KRATOS_TRY

    rDerivative  =  rStrain.Invariants.I3 * this->msIdentityMatrix(a,b);
    rDerivative *= (rStrain.Invariants.I3 * this->msIdentityMatrix(c,d));

    return rDerivative;

    KRATOS_CATCH(" ")
  }



  //************************************************************************************
  //************************************************************************************

  //Invariants 2nd derivatives by components
  double& MooneyRivlinModel::GetI1LeftCauchyGreen2ndDerivative(const StrainData& rStrain,
							       double& rDerivative,
							       const double& a,
							       const double& b,
							       const double& c,
							       const double& d) //ddI1/dbdb
  {
    KRATOS_TRY

    rDerivative = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);

    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI2LeftCauchyGreen2ndDerivative(const StrainData& rStrain,
							       double& rDerivative,
							       const double& a,
							       const double& b,
							       const double& c,
							       const double& d) //ddI2/dbdb
  {
    KRATOS_TRY

    rDerivative  = GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative *= rStrain.Invariants.I1;
    rDerivative += rStrain.Matrix(a,b)*rStrain.Matrix(c,d);
    rDerivative -= (this->msIdentityMatrix(a,c)*this->msIdentityMatrix(b,d)+this->msIdentityMatrix(a,d)*this->msIdentityMatrix(b,c)) * rStrain.Matrix(b,d);
    return rDerivative;

    KRATOS_CATCH(" ")
  }

  double& MooneyRivlinModel::GetI3LeftCauchyGreen2ndDerivative(const StrainData& rStrain,
							       double& rDerivative,
							       const double& a,
							       const double& b,
							       const double& c,
							       const double& d) //ddI3/dbdb
  {
    KRATOS_TRY

    rDerivative  = (-1.0) * GetFourthOrderUnitTensor(rDerivative,a,b,c,d);
    rDerivative += this->msIdentityMatrix(a,b)*this->msIdentityMatrix(c,d);
    rDerivative *= rStrain.Invariants.I3;

    return rDerivative;

    KRATOS_CATCH(" ")
  }




  //************************************************************************************
  //************************************************************************************

  double& MooneyRivlinModel::GetFourthOrderUnitTensor(double& rValue,
						      const double& a,
						      const double& b,
						      const double& c,
						      const double& d) //ddC/dCdC or ddb/dbdb
  {
    KRATOS_TRY

    rValue = ConstitutiveModelUtilities::CalculateFourthOrderUnitTensor(this->msIdentityMatrix,rValue,a,b,c,d);

    return rValue;

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  double& MooneyRivlinModel::GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //dU/dJ
  {
      KRATOS_TRY

      // const ModelDataType&  rValues = rVariables.GetModelData();

      // rDerivative = rValues.GetPressure();

      // return rDerivative;

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //derivative of "U(J) = (K/2)*ln(J)²"
      //dU(J)/dJ = (K)*(lnJ/J)
      rDerivative = rMaterial.GetBulkModulus() * std::log( rVariables.Strain.Invariants.J );

      rDerivative /= rVariables.Strain.Invariants.J;

      return rDerivative;

      KRATOS_CATCH(" ")
  };

  //************************************************************************************
  //************************************************************************************

  double& MooneyRivlinModel::GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) //ddU/dJdJ
  {
      KRATOS_TRY

      // rDerivative = 0.0;

      // return rDerivative;

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      //derivative of "dU(J)/dJ = (K)*(lnJ/J)"
      //ddU(J)/dJdJ = (K)*(1-lnJ)/J²
      rDerivative = rMaterial.GetBulkModulus() * (1.0 -std::log(rVariables.Strain.Invariants.J)) / (rVariables.Strain.Invariants.J * rVariables.Strain.Invariants.J);

      return rDerivative;

      KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  int MooneyRivlinModel::Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    return 0;

    KRATOS_CATCH(" ")
  }



} // Namespace Kratos
