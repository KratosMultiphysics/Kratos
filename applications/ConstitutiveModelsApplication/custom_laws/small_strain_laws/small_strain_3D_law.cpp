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
#include "custom_laws/small_strain_laws/small_strain_3D_law.hpp"
#include "custom_utilities/constitutive_model_utilities.hpp"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SmallStrain3DLaw::SmallStrain3DLaw()
    : Constitutive3DLaw()
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  //******************************CONSTRUCTOR WITH THE MODEL****************************
  //************************************************************************************

  SmallStrain3DLaw::SmallStrain3DLaw(ModelTypePointer pModel)
    : Constitutive3DLaw()
  {
    KRATOS_TRY

    //model
    mpModel = pModel->Clone();

    KRATOS_CATCH(" ")
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  SmallStrain3DLaw::SmallStrain3DLaw(const SmallStrain3DLaw& rOther)
    : Constitutive3DLaw(rOther)
  {
    mpModel = rOther.mpModel->Clone();
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SmallStrain3DLaw& SmallStrain3DLaw::operator=(const SmallStrain3DLaw& rOther)
  {
    Constitutive3DLaw::operator=(rOther);
    mpModel = rOther.mpModel->Clone();
    return *this;
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer SmallStrain3DLaw::Clone() const
  {
    return Kratos::make_shared<SmallStrain3DLaw>(*this);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  SmallStrain3DLaw::~SmallStrain3DLaw()
  {
  }


  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************
  //*************************** SET VALUE: VECTOR **************************************
  //************************************************************************************


  void SmallStrain3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    mpModel->SetValue(rThisVariable,rValue, rCurrentProcessInfo);

    KRATOS_CATCH(" ")
  }



  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

  void SmallStrain3DLaw::InitializeMaterial( const Properties& rMaterialProperties,
					     const GeometryType& rElementGeometry,
					     const Vector& rShapeFunctionsValues )
  {
    KRATOS_TRY

    ConstitutiveLaw::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);

    mpModel->InitializeMaterial(rMaterialProperties);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void SmallStrain3DLaw::InitializeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY

    rModelValues.SetOptions(rValues.GetOptions());
    rModelValues.SetMaterialProperties(rValues.GetMaterialProperties());
    rModelValues.SetProcessInfo(rValues.GetProcessInfo());
    rModelValues.SetVoigtSize(this->GetStrainSize());
    rModelValues.SetVoigtIndexTensor(this->GetVoigtIndexTensor());

    Vector& rStrainVector = rValues.GetStrainVector();
    MatrixType& rStrainMatrix = rModelValues.rStrainMatrix();
    rStrainMatrix = ConstitutiveModelUtilities::StrainVectorToTensor(rStrainVector,rStrainMatrix);

    if( rValues.GetOptions().Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE) )
      rModelValues.State.Set(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES);

    //initialize model
    mpModel->InitializeModel(rModelValues);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void SmallStrain3DLaw::FinalizeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY

    //Finalize Material response
    if(rValues.GetOptions().Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE)){

      //finalize model (update total strain measure)
      mpModel->FinalizeModel(rModelValues);

    }

    KRATOS_CATCH(" ")
  }

  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void SmallStrain3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {
    KRATOS_TRY

    ModelDataType ModelValues;

    this->CalculateMaterialResponseKirchhoff(rValues,ModelValues);

    KRATOS_CATCH(" ")
  }

  void SmallStrain3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    const Flags& rOptions = rValues.GetOptions();

    //1.- Initialize hyperelastic model parameters
    LawDataType& rVariables = rModelValues.rConstitutiveLawData();
    rVariables.StressMeasure = ConstitutiveModelData::StressMeasure_Kirchhoff; //set required stress measure

    this->InitializeModelData(rValues, rModelValues);

    //2.-Calculate domain variables (Temperature, Pressure, Size) and calculate material parameters
    this->CalculateDomainVariables(rValues, rModelValues);

    ConstitutiveModelData::CalculateMaterialParameters(rModelValues);

    //3.-Calculate Total kirchhoff stress and  Constitutive Matrix related to Total Kirchhoff stress

    if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS) && rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      Vector& rStressVector       = rValues.GetStressVector();
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      this->CalculateStressVectorAndConstitutiveMatrix(rModelValues, rStressVector, rConstitutiveMatrix);

    }
    else{

      //4.-Calculate Total Kirchhoff stress

      if(rOptions.Is(ConstitutiveLaw::COMPUTE_STRESS)){

	Vector& rStressVector       = rValues.GetStressVector();
	this->CalculateStressVector(rModelValues, rStressVector);

      }

      //5.-Calculate Constitutive Matrix related to Total Kirchhoff stress

      if(rOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){

      	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	this->CalculateConstitutiveMatrix(rModelValues, rConstitutiveMatrix);

      }

    }

    //6.- Finalize hyperelastic model parameters
    this->FinalizeModelData(rValues,rModelValues);

    // std::cout<<" ConstitutiveMatrix "<<rValues.GetConstitutiveMatrix()<<std::endl;
    // std::cout<<" StrainVector "<<rValues.GetStrainVector()<<std::endl;
    // std::cout<<" StressVector "<<rValues.GetStressVector()<<std::endl;

    //-----------------------------//
    // const Properties& rMaterialProperties  = rValues.GetMaterialProperties();

    // Vector& rStrainVector                  = rValues.GetStrainVector();
    // Vector& rStressVector                  = rValues.GetStressVector();


    // // Calculate total Kirchhoff stress

    // if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

    //   Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

    //   this->CalculateConstitutiveMatrix( rConstitutiveMatrix, rMaterialProperties);

    //   this->CalculateStress( rStrainVector, rConstitutiveMatrix, rStressVector );

    // }
    // else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

    //   Matrix& rConstitutiveMatrix  = rValues.GetConstitutiveMatrix();
    //   this->CalculateConstitutiveMatrix(rConstitutiveMatrix, rMaterialProperties);

    // }
    //-----------------------------//


    KRATOS_CATCH(" ")

  }

  //*******************************COMPUTE STRESS VECTOR********************************
  //************************************************************************************

  void SmallStrain3DLaw::CalculateStressVector(ModelDataType& rModelValues, Vector& rStressVector)
  {
    KRATOS_TRY

    MatrixType StressMatrix;
    StressMatrix.clear();

    if(rModelValues.GetOptions().Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)){
      mpModel->CalculateIsochoricStressTensor(rModelValues, StressMatrix);
    }
    else if(rModelValues.GetOptions().Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){
      mpModel->CalculateVolumetricStressTensor(rModelValues, StressMatrix);
    }
    else{
      mpModel->CalculateStressTensor(rModelValues, StressMatrix);
    }

    rStressVector = ConstitutiveModelUtilities::StressTensorToVector(StressMatrix, rStressVector);

    KRATOS_CATCH(" ")
  }

  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  void SmallStrain3DLaw::CalculateConstitutiveMatrix(ModelDataType& rModelValues, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    //Calculate ConstitutiveMatrix
    if(rModelValues.GetOptions().Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)){

      mpModel->CalculateIsochoricConstitutiveTensor(rModelValues, rConstitutiveMatrix);
    }
    else if(rModelValues.GetOptions().Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){

      mpModel->CalculateVolumetricConstitutiveTensor(rModelValues, rConstitutiveMatrix);
    }
    else{

      mpModel->CalculateConstitutiveTensor(rModelValues, rConstitutiveMatrix);
    }

    KRATOS_CATCH(" ")
  }

  //******************COMPUTE STRESS AND ALGORITHMIC CONSTITUTIVE MATRIX****************
  //************************************************************************************

  void SmallStrain3DLaw::CalculateStressVectorAndConstitutiveMatrix(ModelDataType& rModelValues, Vector& rStressVector, Matrix& rConstitutiveMatrix)
  {
    KRATOS_TRY

    MatrixType StressMatrix;
    StressMatrix.clear();

    //Calculate Stress and ConstitutiveMatrix
    if(rModelValues.GetOptions().Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)){

      mpModel->CalculateIsochoricStressAndConstitutiveTensors(rModelValues, StressMatrix, rConstitutiveMatrix);
    }
    else if(rModelValues.GetOptions().Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){

      mpModel->CalculateVolumetricStressAndConstitutiveTensors(rModelValues, StressMatrix, rConstitutiveMatrix);
    }
    else{

      mpModel->CalculateStressAndConstitutiveTensors(rModelValues, StressMatrix, rConstitutiveMatrix);
    }

    rStressVector = ConstitutiveModelUtilities::StressTensorToVector(StressMatrix, rStressVector);


    KRATOS_CATCH(" ")
  }

  //***********************COMPUTE TOTAL STRESS ****************************************
  //************************************************************************************


  void SmallStrain3DLaw::CalculateStress(const Vector & rStrainVector,
					 const Matrix & rConstitutiveMatrix,
					 Vector& rStressVector)
  {
    KRATOS_TRY

    // Cauchy StressVector
    noalias(rStressVector) = prod(rConstitutiveMatrix,rStrainVector);

    KRATOS_CATCH(" ")

  }



  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************


  void SmallStrain3DLaw::CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix,
						     const Properties& rMaterialProperties)
  {
    KRATOS_TRY

    // Lame constants
    const double& rYoungModulus          = rMaterialProperties[YOUNG_MODULUS];
    const double& rPoissonCoefficient    = rMaterialProperties[POISSON_RATIO];


    // 3D linear elastic constitutive matrix
    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );
    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 3 , 3 ) = rConstitutiveMatrix ( 0 , 0 )*(1.0-2.0*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));
    rConstitutiveMatrix ( 4 , 4 ) = rConstitutiveMatrix ( 3 , 3 );
    rConstitutiveMatrix ( 5 , 5 ) = rConstitutiveMatrix ( 3 , 3 );

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

    rConstitutiveMatrix ( 0 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
    rConstitutiveMatrix ( 2 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

    rConstitutiveMatrix ( 1 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
    rConstitutiveMatrix ( 2 , 1 ) = rConstitutiveMatrix ( 0 , 1 );


    //initialize to zero other values
    for(unsigned int i=0; i<3; i++)
    {
          for(unsigned int j=3; i<6; i++)
          {
            rConstitutiveMatrix ( i , j ) = 0;
            rConstitutiveMatrix ( j , i ) = 0;
          }

    }

    KRATOS_CATCH(" ")
  }



  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void SmallStrain3DLaw::GetLawFeatures(Features& rFeatures)
  {
    KRATOS_TRY

    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void SmallStrain3DLaw::GetModelFeatures(Features& rFeatures)
  {
    KRATOS_TRY

    //Get model variables and set law characteristics
    if( mpModel != NULL ){

      std::vector<Variable<double> > ScalarVariables;
      std::vector<Variable<array_1d<double,3> > > ComponentVariables;

      mpModel->GetDomainVariablesList(ScalarVariables, ComponentVariables);

      for(std::vector<Variable<array_1d<double,3> > >::iterator cv_it=ComponentVariables.begin(); cv_it != ComponentVariables.end(); cv_it++)
	{
	  if( *cv_it == DISPLACEMENT ){
	    for(std::vector<Variable<double> >::iterator sv_it=ScalarVariables.begin(); sv_it != ScalarVariables.end(); sv_it++)
	      {

		if( *sv_it == PRESSURE )
		  rFeatures.mOptions.Set( U_P_LAW );
	      }
	  }
	  // if( *cv_it == VELOCITY ){
	  //   for(std::vector<Variables<double> >::iterator sv_it=ScalarVariables.begin(); sv_it != ScalarVariables.end(); )
	  //     {
	  // 	if( *sv_it == PRESSURE )
	  // 	  rFeatures.mOptions.Set( V_P_LAW );
	  //     }
	  // }
	}

      //...
    }



    KRATOS_CATCH(" ")
  }

} // Namespace Kratos
