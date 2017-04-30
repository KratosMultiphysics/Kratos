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

    //member variables initialization
    mInitialStrainVector.clear();

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
      
    //member variables initialization
    mInitialStrainVector.clear();
    
    KRATOS_CATCH(" ")    
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  SmallStrain3DLaw::SmallStrain3DLaw(const SmallStrain3DLaw& rOther)
    : Constitutive3DLaw(rOther)
    ,mInitialStrainVector(rOther.mInitialStrainVector)
  {
    mpModel = rOther.mpModel->Clone();
  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SmallStrain3DLaw& SmallStrain3DLaw::operator=(const SmallStrain3DLaw& rOther)
  {
    Constitutive3DLaw::operator=(rOther);
    mpModel = rOther.mpModel->Clone();
    mInitialStrainVector = rOther.mInitialStrainVector;
    return *this;
  } 
  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer SmallStrain3DLaw::Clone() const
  {
    return ( SmallStrain3DLaw::Pointer(new SmallStrain3DLaw(*this)) );
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

    // A method to compute the initial linear strain from the stress is needed
    //if(rThisVariable == INITIAL_STRESS_VECTOR)

    // A method to compute the initial linear strain from the stress is needed
    // if(rThisVariable == INITIAL_STRAIN_VECTOR){
    //   mInitialStrainVector = rValue;
    // }
         
    KRATOS_CATCH(" ")
  }

   

  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

 
  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void SmallStrain3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
  {
    KRATOS_TRY
     
    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();
    
    const Properties& rMaterialProperties  = rValues.GetMaterialProperties();    

    Vector& rStrainVector                  = rValues.GetStrainVector();
    Vector& rStressVector                  = rValues.GetStressVector();

    AddInitialStrainVector(rStrainVector);
    
    //-----------------------------//

    
    //2.-Calculate total Kirchhoff stress

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){
      
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
      
      this->CalculateConstitutiveMatrix(rConstitutiveMatrix, rMaterialProperties);
      
      this->CalculateStress(rStrainVector, rConstitutiveMatrix, rStressVector);
      
    }
    else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

      Matrix& rConstitutiveMatrix  = rValues.GetConstitutiveMatrix();
      this->CalculateConstitutiveMatrix(rConstitutiveMatrix, rMaterialProperties);
      
    }
    
    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************


  void SmallStrain3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {    
    KRATOS_TRY

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    
    const Properties& rMaterialProperties  = rValues.GetMaterialProperties();    

    Vector& rStrainVector                  = rValues.GetStrainVector();
    Vector& rStressVector                  = rValues.GetStressVector();

    AddInitialStrainVector(rStrainVector);
    
    //-----------------------------//

    // Calculate total Kirchhoff stress

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){
      
      Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
      
      this->CalculateConstitutiveMatrix( rConstitutiveMatrix, rMaterialProperties);
      
      this->CalculateStress( rStrainVector, rConstitutiveMatrix, rStressVector );
      
    }
    else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

      Matrix& rConstitutiveMatrix  = rValues.GetConstitutiveMatrix();
      this->CalculateConstitutiveMatrix(rConstitutiveMatrix, rMaterialProperties);
      
    }

    KRATOS_CATCH(" ")

  }


  //***********************COMPUTE TOTAL STRESS PK2*************************************
  //************************************************************************************


  void SmallStrain3DLaw::CalculateStress(const Vector & rStrainVector,
					   const Matrix & rConstitutiveMatrix,
					   Vector& rStressVector)
  {
    KRATOS_TRY
          
    //1.-2nd Piola Kirchhoff StressVector or Cauchy StressVector
    if( rStressVector.size() != rStrainVector.size() )
      rStressVector.resize(rStrainVector.size(),false);

    noalias(rStressVector) = prod(rConstitutiveMatrix,rStrainVector);

    KRATOS_CATCH(" ")

  }



  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************


  void SmallStrain3DLaw::CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix,
							const Properties& rMaterialProperties)
  {
    KRATOS_TRY
          
    rConstitutiveMatrix.clear();

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
    
    KRATOS_CATCH(" ")
  }


  //***********************COMPUTE TOTAL STRESS PK2*************************************
  //************************************************************************************


  void SmallStrain3DLaw::AddInitialStrainVector(Vector & rStrainVector)
  {
    KRATOS_TRY
          
    for(unsigned int i=0; i<rStrainVector.size(); i++)
      {
	rStrainVector[i] += mInitialStrainVector[i];
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

    //Get model variables and set law characteristics
    if( mpModel != NULL ){

      std::vector<Variable<double> > ScalarVariables;
      std::vector<Variable<array_1d<double,3> > > ComponentVariables;

      mpModel->GetDomainVariablesList(ScalarVariables, ComponentVariables);
      
      for(std::vector<Variable<array_1d<double,3> > >::iterator cv_it=ComponentVariables.begin(); cv_it != ComponentVariables.end(); )
	{
	  if( *cv_it == DISPLACEMENT ){
	    for(std::vector<Variable<double> >::iterator sv_it=ScalarVariables.begin(); sv_it != ScalarVariables.end(); )
	      {
		if( *sv_it == PRESSURE )
		  rFeatures.mOptions.Set( U_P_LAW );
	      }
	  }
	  // if( *cv_it == VELOCITY ){
	  //   for(std::vector<Variable<double> >::iterator sv_it=ScalarVariables.begin(); sv_it != ScalarVariables.end(); )
	  //     {
	  // 	if( *sv_it == PRESSURE )
	  // 	  rFeatures.mOptions.Set( V_P_LAW );
	  //     }
	  // }
	}

      //...
    }
    
    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
    
    KRATOS_CATCH(" ")
  }



} // Namespace Kratos
