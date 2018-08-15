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
#include "custom_models/constitutive_model.hpp"

namespace Kratos
{

  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModel, ADD_HISTORY_VECTOR,        0 );
  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModel, HISTORY_STRAIN_MEASURE,    1 );
  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModel, HISTORY_STRESS_MEASURE,    2 );

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  ConstitutiveModel::ConstitutiveModel()
      :mHistoryVector(ZeroVector(6))
  {
    //this->mHistoryVector.clear();
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  ConstitutiveModel::ConstitutiveModel(const ConstitutiveModel& rOther)
    :mOptions(rOther.mOptions)
    ,mHistoryVector(rOther.mHistoryVector)
  {
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  ConstitutiveModel& ConstitutiveModel::operator=(const ConstitutiveModel& rOther)
  {
    mOptions = rOther.mOptions;
    mHistoryVector = rOther.mHistoryVector;
    return *this;
  }


  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveModel::Pointer ConstitutiveModel::Clone() const
  {
    return Kratos::make_shared<ConstitutiveModel>(*this);
  }



  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  ConstitutiveModel::~ConstitutiveModel()
  {
  }


  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool ConstitutiveModel::Has(const Variable<double>& rThisVariable)
  {
    KRATOS_TRY

    return false;

    KRATOS_CATCH(" ")
  }


  //***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************
  void ConstitutiveModel::SetValue(const Variable<double>& rVariable, const double& rValue,
				   const ProcessInfo& rCurrentProcessInfo)

  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  void ConstitutiveModel::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue,
				   const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  void ConstitutiveModel::SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue,
				   const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }


  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& ConstitutiveModel::GetValue(const Variable<double>& rThisVariable, double& rValue)
  {
    KRATOS_TRY

    rValue=0;

    return rValue;

    KRATOS_CATCH(" ")
  }


  void ConstitutiveModel::GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
						 std::vector<Variable<array_1d<double,3> > >& rComponentVariables)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the Constitutive Model base class Variables List... illegal operation" << std::endl;

    KRATOS_CATCH(" ")
  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************

  void ConstitutiveModel::InitializeMaterial(const Properties& rMaterialProperties)
  {
    KRATOS_TRY

    //KRATOS_ERROR << "calling ConstitutiveModel Initialize base class " << std::endl;

    KRATOS_CATCH(" ")
  }



  void ConstitutiveModel::InitializeModel(ModelDataType& rValues)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Initialize base class " << std::endl;

    KRATOS_CATCH(" ")
  }


  void ConstitutiveModel::FinalizeModel(ModelDataType& rValues)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

  /**
   * Calculate Strain Energy Density Functions
   */
  void ConstitutiveModel::CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }


  /**
   * Calculate Stresses
   */
  void ConstitutiveModel::CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }

  void ConstitutiveModel::CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }

  void ConstitutiveModel::CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }


  /**
   * Calculate Constitutive Tensor
   */
  void ConstitutiveModel::CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }

  void ConstitutiveModel::CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }

  void ConstitutiveModel::CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }


  /**
   * Calculate Stress and Constitutive Tensor
   */
  void ConstitutiveModel::CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }

  void ConstitutiveModel::CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }

  void ConstitutiveModel::CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel Finalize base class " << std::endl;

    KRATOS_CATCH(" ")
  }


  int ConstitutiveModel::Check(const Properties& rMaterialProperties,
			       const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;

    return 0;

    KRATOS_CATCH(" ")
  }

} // Namespace Kratos
