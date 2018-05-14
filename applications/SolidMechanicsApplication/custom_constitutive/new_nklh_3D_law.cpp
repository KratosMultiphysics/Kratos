//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:               DAbadias $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/new_nklh_3D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   NewNKLH3DLaw::NewNKLH3DLaw()
      : LinearElastic3DLaw()
   {
   }

   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   NewNKLH3DLaw::NewNKLH3DLaw(const NewNKLH3DLaw& rOther)
      : LinearElastic3DLaw(rOther)
   {
   }

   //********************************CLONE***********************************************
   //************************************************************************************

   ConstitutiveLaw::Pointer NewNKLH3DLaw::Clone() const
   {
      NewNKLH3DLaw::Pointer p_clone(new NewNKLH3DLaw(*this));
      return p_clone;
   }

   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   NewNKLH3DLaw::~NewNKLH3DLaw()
   {
   }

   // **********************************************************************************
   // construct the vectors
   void NewNKLH3DLaw::InitializeMaterial( const Properties& rMaterialProperties,
         const GeometryType& rElementGeometry,
         const Vector& rShapeFunctionsValues )
   {
      mStressPrevious = ZeroVector(6);
      mStrainPrevious = ZeroVector(6);
      mHistoricalVariablesPrevious = ZeroVector(31);
   }


   //*******************************  CalculateMaterialResp Kirchhoff *******************
   //************************************************************************************
   void NewNKLH3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
   {

      KRATOS_TRY

      //-----------------------------//

      //a.-Check if the constitutive parameters are passed correctly to the law calculation
      //CheckParameters(rValues);

      //b.- Get Values to compute the constitutive law:
      Flags &Options=rValues.GetOptions();

      const Properties& MaterialProperties  = rValues.GetMaterialProperties();

      const Vector& rStrainVector                  = rValues.GetStrainVector();
      Vector& rStressVector                  = rValues.GetStressVector();

      // Convert strain to 3D and compute the incremental strain (3D)
      Vector StrainVector = ZeroVector(6);
      unsigned int thisSize = 3;
      if ( rStressVector.size() == 3)
         thisSize = 2;
      Matrix StrainTensor = ZeroMatrix(thisSize,thisSize);

      StrainTensor = MathUtils<double>::StrainVectorToTensor( rStrainVector);
      StrainTensor = ConvertStrainTensorTo3D( StrainTensor);
      StrainVector = MathUtils<double>::StrainTensorToVector( StrainTensor, 6);
      Vector IncrementalStrain = (StrainVector - mStrainPrevious);


      Matrix ConstitutiveMatrix( StrainVector.size() ,StrainVector.size());
      noalias(ConstitutiveMatrix) = ZeroMatrix( StrainVector.size() ,StrainVector.size());

      Vector StressVector = ZeroVector(6);

      Vector PlasticVariables = mHistoricalVariablesPrevious;


      //-----------------------------//
      //1.- Lame constants
      const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
      const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];
      CalculateLinearElasticStiffnessMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );


      // Constitutive Model
      Vector IncrementalStress = ZeroVector(6);

      IncrementalStress = prod(ConstitutiveMatrix, IncrementalStrain);
      StressVector = mStressPrevious + IncrementalStress;




      //7.-Calculate total Kirchhoff stress

      if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

         Matrix StressTensor = ZeroMatrix(3,3);

         StressTensor = MathUtils<double>::StressVectorToTensor( StressVector);
         rStressVector = MathUtils<double>::StressTensorToVector( StressTensor, rStressVector.size() );
      }


      if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

         Matrix& rConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

         rConstitutiveMatrix = ConvertConstitutiveMatrixAppropiateSize( rConstitutiveMatrix, ConstitutiveMatrix, rStressVector.size() );

      }

      if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) ) {
         mHistoricalVariablesPrevious = PlasticVariables;
         mStrainPrevious = StrainVector;
         mStressPrevious = StressVector;
      }


      /*{
         std::cout << " INPUT //  OUTPUT " << std::endl;
         std::cout<<" Strain "<<rStrainVector<<std::endl;
         if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){
            Vector& rStressVector                  = rValues.GetStressVector();
            std::cout << " stress " << rStressVector << std::endl;
         }
         if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

            Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
         }
      }*/

      KRATOS_CATCH("")
   }

   // ***************************************************************************************
   // convert the constitutive matrix to the appropiate size
   Matrix & NewNKLH3DLaw::ConvertConstitutiveMatrixAppropiateSize( Matrix & rOutput, const Matrix & rInput, const int voigtSize)
   {

      KRATOS_TRY

      if ( voigtSize == 6) // 3D 
      {
         rOutput = rInput;
      } else if ( voigtSize == 3) // Plane strain
      {
         rOutput(0,0) = rInput(0,0);
         rOutput(0,1) = rInput(0,1);
         rOutput(1,0) = rInput(1,0);
         rOutput(1,1) = rInput(1,1);

         rOutput(2,0) = rInput(3,0);
         rOutput(2,1) = rInput(3,1);         
         rOutput(2,2) = rInput(3,3);

         rOutput(0,2) = rInput(0,3);
         rOutput(1,2) = rInput(1,3);

      } else if ( voigtSize == 4) // axisym
      {
         for (unsigned int i = 0; i < 4; ++i)  {
            for (unsigned int j = 0; j < 4; ++j) {
               rOutput( i, j) = rInput( i, j);
            }
         }
      } else {
         KRATOS_ERROR << " this does not have the appropiate sizes in here " << std::endl;
      }

      return rOutput;



      KRATOS_CATCH("")
   }


   // ***********************************************************************************
   // convert a 2D strain Tensor to 3D
   Matrix & NewNKLH3DLaw::ConvertStrainTensorTo3D( Matrix & rStrainTensor)
   {
      KRATOS_TRY

      if ( rStrainTensor.size1() == 3)
         return rStrainTensor;


      Matrix AuxMatrix = ZeroMatrix(2,2);
      AuxMatrix = rStrainTensor;

      
      rStrainTensor = ZeroMatrix(3,3);
      
      for (unsigned int i = 0; i < 2; i++)
         for (unsigned int j = 0; j < 2; j++)
            rStrainTensor(i,j) = AuxMatrix(i,j);

      return rStrainTensor;


      KRATOS_CATCH("")
   }

   //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
   //************************************************************************************


   void NewNKLH3DLaw::CalculateLinearElasticStiffnessMatrix( Matrix& rConstitutiveMatrix,
         const double &rYoungModulus,
         const double &rPoissonCoefficient )
   {
      rConstitutiveMatrix.clear();

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

   }


   //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
   //************************************************************************************

   void NewNKLH3DLaw::GetLawFeatures(Features& rFeatures)
   {
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

   }

   //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
   //************************************************************************************

   bool NewNKLH3DLaw::CheckParameters(Parameters& rValues)
   {
      return rValues.CheckAllParameters();
   }



   int NewNKLH3DLaw::Check(const Properties& rMaterialProperties,
         const GeometryType& rElementGeometry,
         const ProcessInfo& rCurrentProcessInfo)
   {

      if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
         KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ", "" )

            const double& nu = rMaterialProperties[POISSON_RATIO];
      const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

      if(POISSON_RATIO.Key() == 0 || check==true)
         KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO has Key zero invalid value ", "" )


            if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
               KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY has Key zero or invalid value ", "" )


                  return 0;

   }


} // Namespace Kratos
