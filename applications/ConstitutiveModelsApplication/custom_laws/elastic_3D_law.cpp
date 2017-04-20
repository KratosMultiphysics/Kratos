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
#include "custom_laws/elastic_3D_law.hpp"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  Elastic3DLaw::Elastic3DLaw()
    : ConstitutiveLaw()
  {
  }
  
  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  Elastic3DLaw::Elastic3DLaw(const Elastic3DLaw& rOther)
    : ConstitutiveLaw(rOther)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer Elastic3DLaw::Clone() const
  {
    return ( Elastic3DLaw::Pointer(new Elastic3DLaw(*this)) );
  }



  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  Elastic3DLaw::~Elastic3DLaw()
  {
  }
  
  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool Elastic3DLaw::Has( const Variable<double>& rThisVariable )
  {
    KRATOS_TRY

    return false;
    
    KRATOS_CATCH(" ")
  }

  bool Elastic3DLaw::Has( const Variable<Vector>& rThisVariable )
  {
    KRATOS_TRY

    return false;
    
    KRATOS_CATCH(" ")
  }

  bool Elastic3DLaw::Has( const Variable<Matrix>& rThisVariable )
  {
    KRATOS_TRY

    return false;
    
    KRATOS_CATCH(" ")
  }


  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& Elastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
    KRATOS_TRY

    return( rValue );
    
    KRATOS_CATCH(" ")   
  }

  Vector& Elastic3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
  {
    KRATOS_TRY

    return( rValue );
    
    KRATOS_CATCH(" ")   
  }

  Matrix& Elastic3DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
  {
    KRATOS_TRY

    return( rValue );
    
    KRATOS_CATCH(" ")   
  }


  //***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************


  void Elastic3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY
      
    KRATOS_CATCH(" ")
  }

  void Elastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY
      
    KRATOS_CATCH(" ")
  }

  void Elastic3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY
      
    KRATOS_CATCH(" ")
  }


  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void Elastic3DLaw::InitializeMaterial( const Properties& rMaterialProperties,
					      const GeometryType& rElementGeometry,
					      const Vector& rShapeFunctionsValues )
  {
    KRATOS_TRY
      
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Elastic3DLaw::InitializeSolutionStep( const Properties& rMaterialProperties,
						  const GeometryType& rElementGeometry, //this is just to give the array of nodes
						  const Vector& rShapeFunctionsValues,
						  const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Elastic3DLaw::FinalizeSolutionStep( const Properties& rMaterialProperties,
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
  
  void  Elastic3DLaw::InitializeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY
    
    KRATOS_CATCH(" ")     
  }

  //************************************************************************************
  //************************************************************************************

  void  Elastic3DLaw::FinalizeModelData(Parameters& rValues,ModelDataType& rModelValues)
  {
    KRATOS_TRY
           
    KRATOS_CATCH(" ")      
  }
  
  
  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void  Elastic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
  {
    KRATOS_TRY
 
    this->CalculateMaterialResponseKirchhoff(rValues);

    //1.- Obtain parameters
    Flags & Options                    = rValues.GetOptions();    

    Vector& StressVector               = rValues.GetStressVector();
    Vector& StrainVector               = rValues.GetStrainVector();

    const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
    const double& DeterminantF         = rValues.GetDeterminantF();

    Matrix& ConstitutiveMatrix         = rValues.GetConstitutiveMatrix();

    //2.-Green-Lagrange Strain:
    if(Options.Is(ConstitutiveLaw::COMPUTE_STRAIN))
      {
	TransformStrains (StrainVector, DeformationGradientF, StrainMeasure_Almansi, StrainMeasure_GreenLagrange);
      }

    //3.-Calculate Total PK2 stress
    if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
      {
	TransformStresses(StressVector, DeformationGradientF, DeterminantF, StressMeasure_Kirchhoff, StressMeasure_PK2);
      }

    //4.-Calculate PK2 constitutive tensor
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
      {
	PullBackConstitutiveMatrix(ConstitutiveMatrix, DeformationGradientF);
      }
    
    KRATOS_CATCH(" ")      
  }


  //************************************************************************************
  //************************************************************************************


  void Elastic3DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
  {
    KRATOS_TRY
 
    this->CalculateMaterialResponsePK2 (rValues);

    Vector& rStressVector               = rValues.GetStressVector();
    const Matrix& rDeformationGradientF = rValues.GetDeformationGradientF();
    const double& rDeterminantF         = rValues.GetDeterminantF();

    TransformStresses(rStressVector,rDeformationGradientF,rDeterminantF,StressMeasure_PK2,StressMeasure_PK1);
    
    KRATOS_CATCH(" ")
	  
  }

  //************************************************************************************
  //************************************************************************************

  void Elastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
  {
    KRATOS_TRY
      
    this->CalculateMaterialResponsePK2 (rValues);

    //1.- Obtain parameters
    Flags & Options                    = rValues.GetOptions();    

    Vector& StressVector               = rValues.GetStressVector();
    Vector& StrainVector               = rValues.GetStrainVector();

    const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
    const double& DeterminantF         = rValues.GetDeterminantF();

    Matrix& ConstitutiveMatrix         = rValues.GetConstitutiveMatrix();

    //2.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
      {
        TransformStrains (StrainVector, DeformationGradientF, StrainMeasure_GreenLagrange, StrainMeasure_Almansi);
      }

    //3.-Calculate Total Kirchhoff stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
      {
        TransformStresses(StressVector, DeformationGradientF, DeterminantF, StressMeasure_PK2, StressMeasure_Kirchhoff);
      }

    //4.-Calculate Kirchhoff constitutive tensor
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      {
        PushForwardConstitutiveMatrix(ConstitutiveMatrix, DeformationGradientF);
      }

    KRATOS_CATCH(" ")      
  }


  //************************************************************************************
  //************************************************************************************

  void Elastic3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
  {
    KRATOS_TRY
      
    this->CalculateMaterialResponseKirchhoff (rValues);

    const double& rDeterminantF          = rValues.GetDeterminantF();
    Vector& rStressVector                = rValues.GetStressVector();
    Matrix& rConstitutiveMatrix          = rValues.GetConstitutiveMatrix();

    //Set to cauchy Stress:
    rStressVector       /= rDeterminantF;
    rConstitutiveMatrix /= rDeterminantF;

    KRATOS_CATCH(" ")
  }


  //***********************************UPDATE*******************************************
  //************************************************************************************

  void Elastic3DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
  {
    KRATOS_TRY
      
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK2(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Elastic3DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
  {
    KRATOS_TRY
      
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK1(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void Elastic3DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
  {
    KRATOS_TRY
      
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseKirchhoff(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void Elastic3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
  {
    KRATOS_TRY
      
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseCauchy(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  // void Elastic3DLaw::CalculateStrainEnergy(Parameters& rValues)
  // {
  //   KRATOS_TRY

  //   mpElasticModel->CalculateStrainEnergyFunction(mStrainEnergy);
      
  //   // double ln_J = std::log(ElasticVariables.DeterminantF);
  //   // double trace_C = 0.0;
    
  //   // for (unsigned int i = 0; i<RightCauchyGreen.size1();i++)
  //   //   {
  //   // 	trace_C += RightCauchyGreen(i,i);
  //   //   }
    
  //   // mStrainEnergy =  0.5*ElasticVariables.LameLambda*ln_J*ln_J - ElasticVariables.LameMu*ln_J + 0.5*ElasticVariables.LameMu*(trace_C-3); //see Belytschko page 239

  //   KRATOS_CATCH(" ")
  // }


  //******************************* COMPUTE DOMAIN VARIABLES  **************************
  //************************************************************************************


  void Elastic3DLaw::CalculateDomainVariables (Parameters& rValues, ModelDataType& rModelValues)
  {
    KRATOS_TRY

    ConstitutiveLawDataType& rVariables = rModelValues.rConstitutiveLawData();
    
    rVariables.Temperature = this->CalculateDomainTemperature(rValues, rVariables.Temperature);
    rVariables.Pressure    = this->CalculateDomainPressure(rValues, rVariables.Pressure);

    const GeometryType& rDomainGeometry = rValues.GetElementGeometry();
    ConstitutiveLawUtilities::CalculateCharacteristicSize(rDomainGeometry, rVariables.CharacteristicSize);

    KRATOS_CATCH(" ")
  }

  //****************************** COMPUTE DOMAIN VARIABLE ****************************
  //************************************************************************************
  
  double& Elastic3DLaw::CalculateDomainVariable(Parameters& rValues, const Variable<double>& rThisVariable, double& rVariable)
  {
    KRATOS_TRY
      
    const GeometryType& DomainGeometry = rValues.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rValues.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();
    
    rVariable = 0;    

    if( DomainGeometry[0].SolutionStepsDataHas(rThisVariable) )
      KRATOS_THROW_ERROR( std::logic_error, "Elastic3DLaw : Asking a variable that is not in the SolutionStepsData ..", rThisVariable )
    
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
	  rVariable += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(rThisVariable);
      }
    
    
    return rVariable;

    KRATOS_CATCH(" ")
  }
  
  //******************************* COMPUTE DOMAIN TEMPERATURE  ************************
  //************************************************************************************


  double& Elastic3DLaw::CalculateDomainTemperature (Parameters& rValues, double& rTemperature)
  {
    KRATOS_TRY
      
      //rTemperature = this->CalculateDomainVariable(rValues,TEMPERATURE,rTemperature);
				  
    return rTemperature;

    KRATOS_CATCH(" ")
  }

  //******************************* COMPUTE DOMAIN PRESSURE ****************************
  //************************************************************************************


  double& Elastic3DLaw::CalculateDomainPressure (Parameters& rValues, double& rPressure)
  {
    KRATOS_TRY
      
      //rPressure = this->CalculateDomainVariable(rValues,PRESSURE,rPressure);
				  
    return rPressure;

    KRATOS_CATCH(" ")
  }


  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void Elastic3DLaw::GetLawFeatures(Features& rFeatures)
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

  bool Elastic3DLaw::CheckParameters(Parameters& rValues)
  {
    KRATOS_TRY

    return rValues.CheckAllParameters();

    KRATOS_CATCH(" ")    
  }

  //************************************************************************************
  //************************************************************************************


  int Elastic3DLaw::Check(const Properties& rMaterialProperties,
			       const GeometryType& rElementGeometry,
			       const ProcessInfo& rCurrentProcessInfo)
  {    
    KRATOS_TRY
      
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
      KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ", "" )

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
      KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO has Key zero invalid value ", "" )
		
    return 0;
    
    KRATOS_CATCH(" ")
  }

} // Namespace Kratos
