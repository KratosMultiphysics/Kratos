//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
//#include "custom_elements/large_displacement_element.hpp"
//#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"
//#include "solid_mechanics_application.h"
#include "custom_elements/spatial_lagrangian_U_wP_second_element.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   SpatialLagrangianUwPSecondElement::SpatialLagrangianUwPSecondElement()
      : SpatialLagrangianUwPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPSecondElement::SpatialLagrangianUwPSecondElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : SpatialLagrangianUwPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPSecondElement::SpatialLagrangianUwPSecondElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : SpatialLagrangianUwPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   SpatialLagrangianUwPSecondElement::SpatialLagrangianUwPSecondElement( SpatialLagrangianUwPSecondElement const& rOther)
      :SpatialLagrangianUwPElement(rOther)
       //,mDeterminantF0(rOther.mDeterminantF0)
       //,mDeformationGradientF0(rOther.mDeformationGradientF0)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   SpatialLagrangianUwPSecondElement&  SpatialLagrangianUwPSecondElement::operator=(SpatialLagrangianUwPSecondElement const& rOther)
   {
      SpatialLagrangianUwPElement::operator=(rOther);

      /*    mDeformationGradientF0.clear();
            mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());

            for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
            {
            mDeformationGradientF0[i] = rOther.mDeformationGradientF0[i];
            }

            mDeterminantF0 = rOther.mDeterminantF0;
       */ //MAYBE DONE THERE
      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer SpatialLagrangianUwPSecondElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new SpatialLagrangianUwPSecondElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer SpatialLagrangianUwPSecondElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      SpatialLagrangianUwPSecondElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

      //-----------//

      NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


      if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
         NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

         if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
            KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }

      for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
         NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

      //-----------//

      if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
         NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

      for(unsigned int i=0; i<mDeformationGradientF0.size(); i++)
      {
         NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
      }

      NewElement.mDeterminantF0 = mDeterminantF0;

      return Element::Pointer( new SpatialLagrangianUwPSecondElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPSecondElement::~SpatialLagrangianUwPSecondElement()
   {
   }



   //************************************************************************************
   //************************************************************************************

   int  SpatialLagrangianUwPSecondElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

         int correct = 0;

      correct = LargeDisplacementElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_ERROR( std::logic_error, "constitutive law is not compatible with the U-wP element type ", " SpatialLagrangianUwPElement" )

            //verify that the variables are correctly initialized

            if ( PRESSURE.Key() == 0 )
               KRATOS_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

                  return correct;

      KRATOS_CATCH( "" );
   }


   void SpatialLagrangianUwPSecondElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {

      /*if ( this->Id() == 1) {
        std::cout <<  " ID1 rHS RHS rhs R: detF0 " << rVariables.detF0 << " detF " << rVariables.detF << std::endl;
        }*/
      // EM FALTA UN CACHO
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contribution of the internal and external forces
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
      SpatialLagrangianUwPElement::CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

      // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
      SpatialLagrangianUwPElement::CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
      this->CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

      rVariables.detF0 *= rVariables.detF;
      DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
      SpatialLagrangianUwPElement::CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rIntegrationWeight);

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;
      //KRATOS_WATCH( rRightHandSideVector )
   }



   void SpatialLagrangianUwPSecondElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
   {


      /* std::cout << " I HAVE IT ALL WRONG: detF0 " << 1.0 - rVariables.detF0 << " detF " << 1.0 - rVariables.detF << std::endl;
         std::cout << " PART 2 F  " << rVariables.F << std::endl;
         std::cout << " PART 3 F0 " << rVariables.F0 << std::endl;
         std::cout << std::endl;
       */


      // EM FALTA UN CACHO
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      //contributions of the stiffness matrix calculated on the reference configuration
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      // operation performed: add Km to the rLefsHandSideMatrix

      //respect to the current configuration n+1
      SpatialLagrangianUwPElement::CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kg to the rLefsHandSideMatrix
      SpatialLagrangianUwPElement::CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kup to the rLefsHandSideMatrix
      SpatialLagrangianUwPElement::CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      // operation performed: add Kpu to the rLefsHandSideMatrix
      this->CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      rVariables.detF0 *= rVariables.detF;
      DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      // operation performed: add Kpp to the rLefsHandSideMatrix
      SpatialLagrangianUwPElement::CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      // operation performed: add Kpp Stab to the rLefsHandSideMatrix
      SpatialLagrangianUwPElement::CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;




      //KRATOS_WATCH( rLeftHandSideMatrix )
   }




   void SpatialLagrangianUwPSecondElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         double ScalingConstant = GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;


      Matrix K = ZeroMatrix(dimension);
      Vector b = ZeroVector(dimension);
      double WaterDensity = GetProperties()[DENSITY_WATER];
      b(dimension-1) = -10.0*WaterDensity;
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = Permeability;
      // WHY??

      double detFend = rVariables.detF0 * rVariables.detF; 

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         rRightHandSideVector[indexp] -= rVariables.N[i] * (detFend - rVariables.detF0) / pow( detFend, 2.0) * rIntegrationWeight * ScalingConstant;

         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            // consistent=1;
            // if(i==j)
            //     consistent=2;

            const double& CurrentPressure   = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE );
            const double& PreviousPressure = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE , 1);
            double DeltaPressure = CurrentPressure - PreviousPressure;

            //const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[j].FastGetSolutionStepValue( DISPLACEMENT );
            //const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[j].FastGetSolutionStepValue( DISPLACEMENT , 1 );
            //array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;


            if (false)
               rRightHandSideVector[indexp] += (1.0/WaterBulk) * rVariables.N[i] * rVariables.N[j] * DeltaPressure * rIntegrationWeight * ScalingConstant / detFend;



            for ( unsigned int p = 0; p < dimension; ++p )
            {

               //rRightHandSideVector[indexp] -= rVariables.N[i]*DeltaDisplacement[p] * rVariables.DN_DX(j,p) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
               for ( unsigned int q = 0; q < dimension; ++q )
               {

                  rRightHandSideVector[indexp] += DeltaTime * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q)* CurrentPressure *rIntegrationWeight * ScalingConstant / detFend;

                  if (j == 0) {
                     // Gravity term of the DARCY FLOW.
                     rRightHandSideVector[indexp] += DeltaTime * rVariables.DN_DX(i,p) * K(p,q)*b(q)*rIntegrationWeight* ScalingConstant / detFend;
                  }

               }


            }

         }

         indexp += (dimension + 1);

      }


      // std::cout<<std::endl;
      // std::cout<<" auxiliar " <<auxiliar<<" F0 "<<rVariables.detF0<<std::endl;
      // std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )

   }



   //************************************************************************************
   //************************************************************************************


   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPSecondElement::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY
         double ScalingConstant = GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      Matrix K = ZeroMatrix(dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = Permeability;

      MatrixType Kh=rLeftHandSideMatrix;

      //std::cout << " JUST TO CHEK: F0 " << rVariables.detF0 << " detF " << rVariables.detF << std::endl;
      //std::cout << " JUST TO CHEK: F0 " << (1.0 - rVariables.detF0) << " detF " << ( 1.0 -  rVariables.detF ) << std::endl;

      double detFend = rVariables.detF0 * rVariables.detF; 
      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexp = dimension;

      //double auxiliar = (rVariables.detF0*rVariables.detF0 + 1)/(rVariables.detF0*rVariables.detF0); //(J²-1)
      //double auxiliar = (1.0-std::log(rVariables.detF0))/(rVariables.detF0*rVariables.detF0);   //(ln(J))
      for ( unsigned int i = 0; i < number_of_nodes; ++i )
      {
         for ( unsigned int j = 0; j < number_of_nodes; ++j )
         {
            int indexup = dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexp, indexup+k) += rVariables.detF0 * rVariables.N[i] * rVariables.DN_DX( j , k ) * rIntegrationWeight * ScalingConstant / pow(detFend, 2.0);

            }
         }
         indexp += (dimension + 1);
      }



      // std::cout<<std::endl;
      // std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;


      KRATOS_CATCH( "" )
   }



   //*********************** SAVE *******************************************************
   //************************************************************************************

   void SpatialLagrangianUwPSecondElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   //*************************************  LOAD ****************************************
   //************************************************************************************

   void SpatialLagrangianUwPSecondElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }

}
