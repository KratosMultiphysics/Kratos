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
#include "custom_laws/linear_elastic_laws/linear_elastic_3D_law.hpp"
#include "custom_utilities/constitutive_law_utilities.hpp"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LinearElastic3DLaw::LinearElastic3DLaw()
    : Elastic3DLaw()
  {
    KRATOS_TRY

    //member variables initialization
    mInitialStrainVector.clear();

    KRATOS_CATCH(" ")
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LinearElastic3DLaw::LinearElastic3DLaw(const LinearElastic3DLaw& rOther)
    : Elastic3DLaw(rOther)
    ,mInitialStrainVector(rOther.mInitialStrainVector)
  {
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer LinearElastic3DLaw::Clone() const
  {
    return ( LinearElastic3DLaw::Pointer(new LinearElastic3DLaw(*this)) );
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LinearElastic3DLaw::~LinearElastic3DLaw()
  {
  }


  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************
  //*************************** SET VALUE: VECTOR **************************************
  //************************************************************************************

  
  void LinearElastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
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


  void  LinearElastic3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
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
      
      this->CalculateLinearElasticMatrix(rConstitutiveMatrix, rMaterialProperties);
      
      this->CalculateStress(rStrainVector, rConstitutiveMatrix, rStressVector);
      
    }
    else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

      Matrix& rConstitutiveMatrix  = rValues.GetConstitutiveMatrix();
      this->CalculateLinearElasticMatrix(rConstitutiveMatrix, rMaterialProperties);
      
    }
    
    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************


  void LinearElastic3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
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
      
      this->CalculateLinearElasticMatrix( rConstitutiveMatrix, rMaterialProperties);
      
      this->CalculateStress( rStrainVector, rConstitutiveMatrix, rStressVector );
      
    }
    else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

      Matrix& rConstitutiveMatrix  = rValues.GetConstitutiveMatrix();
      this->CalculateLinearElasticMatrix(rConstitutiveMatrix, rMaterialProperties);
      
    }

    KRATOS_CATCH(" ")

  }


  //***********************COMPUTE TOTAL STRESS PK2*************************************
  //************************************************************************************


  void LinearElastic3DLaw::CalculateStress(const Vector & rStrainVector,
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


  void LinearElastic3DLaw::CalculateLinearElasticMatrix(Matrix& rConstitutiveMatrix,
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


  void LinearElastic3DLaw::AddInitialStrainVector(Vector & rStrainVector)
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

  void LinearElastic3DLaw::GetLawFeatures(Features& rFeatures)
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



} // Namespace Kratos
