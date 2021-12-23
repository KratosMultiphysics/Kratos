//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   March 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_laws/strain_rate_laws/newtonian_3D_law.hpp"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  NewtonianFluid3DLaw::NewtonianFluid3DLaw()
    : ConstitutiveLaw()
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  NewtonianFluid3DLaw::NewtonianFluid3DLaw(const NewtonianFluid3DLaw& rOther)
    : ConstitutiveLaw(rOther)
  {

  }

  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  NewtonianFluid3DLaw& NewtonianFluid3DLaw::operator=(const NewtonianFluid3DLaw& rOther)
  {
    ConstitutiveLaw::operator=(rOther);
    return *this;
  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer NewtonianFluid3DLaw::Clone() const
  {
    return Kratos::make_shared<NewtonianFluid3DLaw>(*this);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  NewtonianFluid3DLaw::~NewtonianFluid3DLaw()
  {
  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************

  void NewtonianFluid3DLaw::InitializeMaterial( const Properties& rProperties,
                                           const GeometryType& rElementGeometry,
                                           const Vector& rShapeFunctionsValues )
  {
    KRATOS_TRY

    KRATOS_CATCH(" ")
  }

  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void NewtonianFluid3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
  {
    KRATOS_TRY

    //0.- Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    const Flags& rOptions = rValues.GetOptions();


    //2.-Calculate Total kirchhoff stress and  Constitutive Matrix related to Cauchy stress

    const Properties& rProperties  = rValues.GetMaterialProperties();

    // Calculate total Kirchhoff stress

    if( rOptions.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

      Vector& rStrainVector                  = rValues.GetStrainVector();
      Vector& rStressVector                  = rValues.GetStressVector();
      this->CalculateStress(rStressVector, rStrainVector, rProperties);

    }
    if( rOptions.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

      Matrix& rConstitutiveMatrix  = rValues.GetConstitutiveMatrix();
      this->CalculateConstitutiveMatrix(rConstitutiveMatrix, rProperties);

    }

    // std::cout<<" ConstitutiveMatrix "<<rValues.GetConstitutiveMatrix()<<std::endl;
    // std::cout<<" StrainVector "<<rValues.GetStrainVector()<<std::endl;
    // std::cout<<" StressVector "<<rValues.GetStressVector()<<std::endl;


    KRATOS_CATCH(" ")

  }

  //*******************************COMPUTE STRESS VECTOR********************************
  //************************************************************************************


void NewtonianFluid3DLaw::CalculateStress(Vector& rStressVector,
                                       const Vector & rStrainVector,
                                       const Properties& rProperties)
  {
    KRATOS_TRY

    const double& rViscosity = rProperties[DYNAMIC_VISCOSITY];

    const double pressure = (rStrainVector[0]+rStrainVector[1]+rStrainVector[2])/3.0;

    // Cauchy StressVector
    rStressVector[0] = 2.0*rViscosity*(rStrainVector[0] - pressure);
    rStressVector[1] = 2.0*rViscosity*(rStrainVector[1] - pressure);
    rStressVector[2] = 2.0*rViscosity*(rStrainVector[2] - pressure);
    rStressVector[3] = rViscosity*rStrainVector[3];
    rStressVector[4] = rViscosity*rStrainVector[4];
    rStressVector[5] = rViscosity*rStrainVector[5];


    KRATOS_CATCH(" ")

  }


  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************


  void NewtonianFluid3DLaw::CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                   const Properties& rProperties)
  {
    KRATOS_TRY

    // Viscosity
    const double& rViscosity = rProperties[DYNAMIC_VISCOSITY];

    const double diagonal_component = 4.0 * rViscosity / 3.0;
    const double side_component = -0.5 * diagonal_component;

    // 3D linear elastic constitutive matrix
    rConstitutiveMatrix ( 0 , 0 ) = diagonal_component;
    rConstitutiveMatrix ( 1 , 1 ) = diagonal_component;
    rConstitutiveMatrix ( 2 , 2 ) = diagonal_component;

    rConstitutiveMatrix ( 3 , 3 ) = rViscosity;
    rConstitutiveMatrix ( 4 , 4 ) = rViscosity;
    rConstitutiveMatrix ( 5 , 5 ) = rViscosity;

    rConstitutiveMatrix ( 0 , 1 ) = side_component;
    rConstitutiveMatrix ( 1 , 0 ) = side_component;

    rConstitutiveMatrix ( 0 , 2 ) = side_component;
    rConstitutiveMatrix ( 2 , 0 ) = side_component;

    rConstitutiveMatrix ( 1 , 2 ) = side_component;
    rConstitutiveMatrix ( 2 , 1 ) = side_component;

    //initialize to zero other values
    for(unsigned int i=0; i<3; ++i)
    {
      for(unsigned int j=3; j<6; ++j)
      {
        rConstitutiveMatrix ( i , j ) = 0;
        rConstitutiveMatrix ( j , i ) = 0;
      }

    }

    rConstitutiveMatrix ( 3 , 4 ) = 0.0;
    rConstitutiveMatrix ( 3 , 5 ) = 0.0;
    rConstitutiveMatrix ( 4 , 5 ) = 0.0;

    rConstitutiveMatrix ( 4 , 3 ) = 0.0;
    rConstitutiveMatrix ( 5 , 3 ) = 0.0;
    rConstitutiveMatrix ( 5 , 4 ) = 0.0;

    KRATOS_CATCH(" ")
  }



  //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
  //************************************************************************************

  void NewtonianFluid3DLaw::GetLawFeatures(Features& rFeatures)
  {
    KRATOS_TRY

    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Velocity_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    KRATOS_CATCH(" ")
  }

  //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
  //************************************************************************************

  int NewtonianFluid3DLaw::Check(const Properties& rProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    if( rProperties[DYNAMIC_VISCOSITY] <= 0.00 )
      KRATOS_ERROR << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for NewtonianFluid3DLaw: " << rProperties[DYNAMIC_VISCOSITY] << std::endl;

    return 0;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************

  void NewtonianFluid3DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
  {
    KRATOS_TRY

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseCauchy(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    KRATOS_CATCH(" ")
  }


} // Namespace Kratos
