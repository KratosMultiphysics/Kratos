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
#include "custom_laws/constitutive_3D_law.hpp"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  Constitutive3DLaw::Constitutive3DLaw() : ConstitutiveLaw()
  {
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  Constitutive3DLaw::Constitutive3DLaw(const Constitutive3DLaw& rOther) : ConstitutiveLaw(rOther)
  {
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  Constitutive3DLaw& Constitutive3DLaw::operator=(const Constitutive3DLaw& rOther)
  {
    return *this;
  }


  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer Constitutive3DLaw::Clone() const
  {
    return Kratos::make_shared<Constitutive3DLaw>(*this);
  }



  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  Constitutive3DLaw::~Constitutive3DLaw()
  {
  }

  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool Constitutive3DLaw::Has( const Variable<double>& rThisVariable )
  {
    KRATOS_TRY

    return false;

    KRATOS_CATCH(" ")
  }

  bool Constitutive3DLaw::Has( const Variable<Vector>& rThisVariable )
  {
    KRATOS_TRY

    return false;

    KRATOS_CATCH(" ")
  }

  bool Constitutive3DLaw::Has( const Variable<Matrix>& rThisVariable )
  {
    KRATOS_TRY

    return false;

    KRATOS_CATCH(" ")
  }

  bool Constitutive3DLaw::Has( const Variable<array_1d<double,3> >& rThisVariable )
  {
    KRATOS_TRY

    return false;

    KRATOS_CATCH(" ")
  }

  bool Constitutive3DLaw::Has( const Variable<array_1d<double,6> >& rThisVariable )
  {
    KRATOS_TRY

    return false;

    KRATOS_CATCH(" ")
  }


  //***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************


  void Constitutive3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  void Constitutive3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  void Constitutive3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  void Constitutive3DLaw::SetValue( const Variable<array_1d<double,3> >& rThisVariable, const array_1d<double,3>& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }


  void Constitutive3DLaw::SetValue( const Variable<array_1d<double,6> >& rThisVariable, const array_1d<double,6>& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& Constitutive3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
    KRATOS_TRY

    return rValue;

    KRATOS_CATCH(" ")
  }

  Vector& Constitutive3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
  {
    KRATOS_TRY

    return rValue;

    KRATOS_CATCH(" ")
  }

  Matrix& Constitutive3DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
  {
    KRATOS_TRY

    return rValue;

    KRATOS_CATCH(" ")
  }

  array_1d<double,3>& Constitutive3DLaw::GetValue( const Variable<array_1d<double,3> >& rThisVariable, array_1d<double,3>& rValue )
  {
    KRATOS_TRY

    return rValue;

    KRATOS_CATCH(" ")
  }

  array_1d<double,6>& Constitutive3DLaw::GetValue( const Variable<array_1d<double,6> >& rThisVariable, array_1d<double,6>& rValue )
  {
    KRATOS_TRY

    return rValue;

    KRATOS_CATCH(" ")
  }

  //***********************CALCULATE VALUE: DOUBLE - VECTOR - MATRIX********************
  //************************************************************************************

  int& Constitutive3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<int>& rThisVariable, int& rValue)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    ModelValues.SetIntVariableData(rThisVariable,rValue);

    this->CalculateValue(rParameterValues,ModelValues);

    return rValue;

    KRATOS_CATCH(" ")
  }

  double& Constitutive3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    ModelValues.SetDoubleVariableData(rThisVariable,rValue);

    this->CalculateValue(rParameterValues,ModelValues);

    return rValue;

    KRATOS_CATCH(" ")
  }

  Vector& Constitutive3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    ModelValues.SetVectorVariableData(rThisVariable,rValue);

    this->CalculateValue(rParameterValues,ModelValues);

    return rValue;

    KRATOS_CATCH(" ")
  }

  Matrix& Constitutive3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    ModelValues.SetMatrixVariableData(rThisVariable,rValue);

    this->CalculateValue(rParameterValues,ModelValues);

    return rValue;

    KRATOS_CATCH(" ")
  }

  array_1d<double, 3 > & Constitutive3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double,3> >& rThisVariable, array_1d<double,3> & rValue)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    ModelValues.SetArray3VariableData(rThisVariable,rValue);

    this->CalculateValue(rParameterValues,ModelValues);

    return rValue;

    KRATOS_CATCH(" ")
  }

  array_1d<double, 6 > & Constitutive3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double,6> >& rThisVariable, array_1d<double,6> & rValue)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    ModelValues.SetArray6VariableData(rThisVariable,rValue);

    this->CalculateValue(rParameterValues,ModelValues);

    return rValue;

    KRATOS_CATCH(" ")
  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,
					      const GeometryType& rElementGeometry,
					      const Vector& rShapeFunctionsValues )
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Constitutive3DLaw::InitializeSolutionStep( const Properties& rMaterialProperties,
						  const GeometryType& rElementGeometry, //this is just to give the array of nodes
						  const Vector& rShapeFunctionsValues,
						  const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Constitutive3DLaw::FinalizeSolutionStep( const Properties& rMaterialProperties,
						const GeometryType& rElementGeometry, //this is just to give the array of nodes
						const Vector& rShapeFunctionsValues,
						const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::InitializeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::FinalizeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }


  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void Constitutive3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponseKirchhoff(rValues);

    //1.- Obtain parameters
    Flags & Options                      = rValues.GetOptions();

    Vector& StressVector                 = rValues.GetStressVector();

    const Matrix& DeltaDeformationMatrix = rValues.GetDeformationGradientF();
    const double& DeltaDeformationDet    = rValues.GetDeterminantF();

    Matrix& ConstitutiveMatrix           = rValues.GetConstitutiveMatrix();


    //2.-Calculate Total PK2 stress
    if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
      {
	TransformStresses(StressVector, DeltaDeformationMatrix, DeltaDeformationDet, StressMeasure_Kirchhoff, StressMeasure_PK2);
      }

    //3.-Calculate PK2 constitutive tensor
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
      {
	PullBackConstitutiveMatrix(ConstitutiveMatrix, DeltaDeformationMatrix);
      }

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************


  void Constitutive3DLaw::CalculateMaterialResponsePK1(Parameters& rValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponsePK2(rValues);

    Vector& rStressVector                 = rValues.GetStressVector();
    const Matrix& rDeltaDeformationMatrix = rValues.GetDeformationGradientF();
    const double& rDeltaDeformationDet    = rValues.GetDeterminantF();

    TransformStresses(rStressVector,rDeltaDeformationMatrix,rDeltaDeformationDet,StressMeasure_PK2,StressMeasure_PK1);

    KRATOS_CATCH(" ")

  }

  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponsePK2(rValues);

    //1.- Obtain parameters
    Flags & Options                      = rValues.GetOptions();

    Vector& StressVector                 = rValues.GetStressVector();

    const Matrix& DeltaDeformationMatrix = rValues.GetDeformationGradientF();
    const double& DeltaDeformationDet    = rValues.GetDeterminantF();

    Matrix& ConstitutiveMatrix           = rValues.GetConstitutiveMatrix();


    //2.-Calculate Total Kirchhoff stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
      {
        TransformStresses(StressVector, DeltaDeformationMatrix, DeltaDeformationDet, StressMeasure_PK2, StressMeasure_Kirchhoff);
      }

    //3.-Calculate Kirchhoff constitutive tensor
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      {
        PushForwardConstitutiveMatrix(ConstitutiveMatrix, DeltaDeformationMatrix);
      }

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponseKirchhoff(rValues);

    const double& rDeltaDeformationDet   = rValues.GetDeterminantF();
    Vector& rStressVector                = rValues.GetStressVector();
    Matrix& rConstitutiveMatrix          = rValues.GetConstitutiveMatrix();

    //Set to cauchy Stress:
    rStressVector       /= rDeltaDeformationDet;
    rConstitutiveMatrix /= rDeltaDeformationDet;

    KRATOS_CATCH(" ")
  }


  //***************************MATERIAL RESPONSES WITH MODELS***************************
  //************************************************************************************

  void Constitutive3DLaw::CalculateValue(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponseKirchhoff(rValues,rModelValues);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::CalculateMaterialResponsePK2(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponseKirchhoff(rValues, rModelValues);

    //1.- Obtain parameters
    Flags & Options                      = rValues.GetOptions();

    Vector& StressVector                 = rValues.GetStressVector();

    const Matrix& DeltaDeformationMatrix = rValues.GetDeformationGradientF();
    const double& DeltaDeformationDet    = rValues.GetDeterminantF();

    Matrix& ConstitutiveMatrix           = rValues.GetConstitutiveMatrix();


    //2.-Calculate Total PK2 stress
    if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
      {
	TransformStresses(StressVector, DeltaDeformationMatrix, DeltaDeformationDet, StressMeasure_Kirchhoff, StressMeasure_PK2);
      }

    //3.-Calculate PK2 constitutive tensor
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
      {
	PullBackConstitutiveMatrix(ConstitutiveMatrix, DeltaDeformationMatrix);
      }

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************


  void Constitutive3DLaw::CalculateMaterialResponsePK1(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponsePK2(rValues, rModelValues);

    Vector& rStressVector                 = rValues.GetStressVector();
    const Matrix& rDeltaDeformationMatrix = rValues.GetDeformationGradientF();
    const double& rDeltaDeformationDet    = rValues.GetDeterminantF();

    TransformStresses(rStressVector,rDeltaDeformationMatrix,rDeltaDeformationDet,StressMeasure_PK2,StressMeasure_PK1);

    KRATOS_CATCH(" ")

  }

  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponsePK2(rValues, rModelValues);

    //1.- Obtain parameters
    Flags & Options                      = rValues.GetOptions();

    Vector& StressVector                 = rValues.GetStressVector();

    const Matrix& DeltaDeformationMatrix = rValues.GetDeformationGradientF();
    const double& DeltaDeformationDet    = rValues.GetDeterminantF();

    Matrix& ConstitutiveMatrix           = rValues.GetConstitutiveMatrix();


    //2.-Calculate Total Kirchhoff stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
      {
        TransformStresses(StressVector, DeltaDeformationMatrix, DeltaDeformationDet, StressMeasure_PK2, StressMeasure_Kirchhoff);
      }

    //3.-Calculate Kirchhoff constitutive tensor
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      {
        PushForwardConstitutiveMatrix(ConstitutiveMatrix, DeltaDeformationMatrix);
      }

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    this->CalculateMaterialResponseKirchhoff(rValues, rModelValues);

    const double& rDeltaDeformationDet   = rValues.GetDeterminantF();
    Vector& rStressVector                = rValues.GetStressVector();
    Matrix& rConstitutiveMatrix          = rValues.GetConstitutiveMatrix();

    //Set to cauchy Stress:
    rStressVector       /= rDeltaDeformationDet;
    rConstitutiveMatrix /= rDeltaDeformationDet;

    KRATOS_CATCH(" ")
  }


  //***********************************INITIALIZE***************************************
  //************************************************************************************

  void Constitutive3DLaw::InitializeMaterialResponsePK2(Parameters& rValues)
  {
    KRATOS_TRY

    // nothing to be done
    // rValues.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);
    // rValues.Reset(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Constitutive3DLaw::InitializeMaterialResponsePK1(Parameters& rValues)
  {
    KRATOS_TRY

    // nothing to be done

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Constitutive3DLaw::InitializeMaterialResponseKirchhoff(Parameters& rValues)
  {
    KRATOS_TRY

    // nothing to be done

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::InitializeMaterialResponseCauchy(Parameters& rValues)
  {
    KRATOS_TRY

    // nothing to be done

    KRATOS_CATCH(" ")
  }



  //***********************************FINALIZE*****************************************
  //************************************************************************************

  void Constitutive3DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
  {
    KRATOS_TRY

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK2(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Constitutive3DLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
  {
    KRATOS_TRY

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK1(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Constitutive3DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
  {
    KRATOS_TRY

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseKirchhoff(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void Constitutive3DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
  {
    KRATOS_TRY

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseCauchy(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  // void Constitutive3DLaw::CalculateStrainEnergy(Parameters& rValues)
  // {
  //   KRATOS_TRY

  //   mpElasticModel->CalculateStrainEnergyFunction(mStrainEnergy);

  //   // double ln_J = std::log(rVariables.DeltaDeformationDet);
  //   // double trace_C = 0.0;

  //   // for(unsigned int i = 0; i<RightCauchyGreen.size1();i++)
  //   //   {
  //   // 	trace_C += RightCauchyGreen(i,i);
  //   //   }

  //   // mStrainEnergy =  0.5*ElasticVariables.LameLambda*ln_J*ln_J - ElasticVariables.LameMu*ln_J + 0.5*ElasticVariables.LameMu*(trace_C-3); //see Belytschko page 239

  //   KRATOS_CATCH(" ")
  // }


  //******************************* COMPUTE DOMAIN VARIABLES  **************************
  //************************************************************************************


  void Constitutive3DLaw::CalculateDomainVariables(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    LawDataType& rVariables = rModelValues.rConstitutiveLawData();

    rVariables.Temperature = this->CalculateDomainTemperature(rValues, rVariables.Temperature);
    rVariables.Pressure    = this->CalculateDomainPressure(rValues, rVariables.Pressure);

    const GeometryType& rDomainGeometry = rValues.GetElementGeometry();
    ConstitutiveModelUtilities::CalculateCharacteristicSize(rDomainGeometry, rVariables.CharacteristicSize);

    KRATOS_CATCH(" ")
  }

  //****************************** COMPUTE DOMAIN VARIABLE ****************************
  //************************************************************************************

  double& Constitutive3DLaw::CalculateDomainVariable(Parameters& rValues, const Variable<double>& rThisVariable, double& rVariable)
  {
    KRATOS_TRY

    const GeometryType& DomainGeometry = rValues.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rValues.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    rVariable = 0;

    if( DomainGeometry[0].SolutionStepsDataHas(rThisVariable) ){

      for( unsigned int j = 0; j < number_of_nodes; j++ )
	{
	  rVariable += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(rThisVariable);
	}
    }
    else{

      rVariable = 0;
      //KRATOS_ERROR << "Constitutive3DLaw : Asking a variable that is not in the SolutionStepsData .." << rThisVariable << std::endl;

    }

    return rVariable;

    KRATOS_CATCH(" ")
  }

  //******************************* COMPUTE DOMAIN TEMPERATURE  ************************
  //************************************************************************************


  double& Constitutive3DLaw::CalculateDomainTemperature(Parameters& rValues, double& rTemperature)
  {
    KRATOS_TRY

    rTemperature = this->CalculateDomainVariable(rValues,TEMPERATURE,rTemperature);

    return rTemperature;

    KRATOS_CATCH(" ")
  }

  //******************************* COMPUTE DOMAIN PRESSURE ****************************
  //************************************************************************************


  double& Constitutive3DLaw::CalculateDomainPressure(Parameters& rValues, double& rPressure)
  {
    KRATOS_TRY

    rPressure = this->CalculateDomainVariable(rValues,PRESSURE,rPressure);

    return rPressure;

    KRATOS_CATCH(" ")
  }


  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void Constitutive3DLaw::GetLawFeatures(Features& rFeatures)
  {
    KRATOS_TRY

    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( ISOTROPIC );


    //Set strain measure required by the consitutive law

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    KRATOS_CATCH(" ")
  }


  //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
  //************************************************************************************

  bool Constitutive3DLaw::CheckParameters(Parameters& rValues)
  {
    KRATOS_TRY

    return rValues.CheckAllParameters();

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  int Constitutive3DLaw::Check(const Properties& rMaterialProperties,
			       const GeometryType& rElementGeometry,
			       const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    return 0;

    KRATOS_CATCH(" ")
  }

} // Namespace Kratos
