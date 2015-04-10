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
#include "custom_elements/spatial_lagrangian_U_wP_FIC_element.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   SpatialLagrangianUwPFICElement::SpatialLagrangianUwPFICElement()
      : SpatialLagrangianUwPStabElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPFICElement::SpatialLagrangianUwPFICElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : SpatialLagrangianUwPStabElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPFICElement::SpatialLagrangianUwPFICElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : SpatialLagrangianUwPStabElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   SpatialLagrangianUwPFICElement::SpatialLagrangianUwPFICElement( SpatialLagrangianUwPFICElement const& rOther)
      :SpatialLagrangianUwPStabElement(rOther)
       //,mDeterminantF0(rOther.mDeterminantF0)
       //,mDeformationGradientF0(rOther.mDeformationGradientF0)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   SpatialLagrangianUwPFICElement&  SpatialLagrangianUwPFICElement::operator=(SpatialLagrangianUwPFICElement const& rOther)
   {
      SpatialLagrangianUwPStabElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer SpatialLagrangianUwPFICElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new SpatialLagrangianUwPFICElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer SpatialLagrangianUwPFICElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      SpatialLagrangianUwPFICElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new SpatialLagrangianUwPFICElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPFICElement::~SpatialLagrangianUwPFICElement()
   {
   }



   //************************************************************************************
   //************************************************************************************

   int  SpatialLagrangianUwPFICElement::Check( const ProcessInfo& rCurrentProcessInfo )
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





   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPFICElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant ; //= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
      double Permeability; double WaterBulk; double DeltaTime;
      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;

      Vector he = ZeroVector(dimension);
      double Factor = GetProperties()[STABILIZATION];
      for (unsigned int i = 0; i < dimension; i++)
      {
         he(i) = Factor*sqrt( 4.0*rIntegrationWeight/ sqrt(3.0) ); 
      }


      double Caux; 
      ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;
      GetValueOnIntegrationPoints( M_MODULUS, Mmodulus, CurrentProcessInfo);
      Caux = 1.0/Mmodulus[0];

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            const double& CurrentPressure   = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE );
            const double& PreviousPressure = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE , 1);
            double DeltaPressure = CurrentPressure - PreviousPressure;

            for ( unsigned int k = 0; k < dimension; ++k )
            {

               rRightHandSideVector[indexp] += rVariables.N[i]* Caux*(1.0/2.0) * he(k) *  rVariables.DN_DX(j,k) * DeltaPressure * rIntegrationWeight * ScalingConstant / rVariables.detF0;

            }

         }

         indexp += (dimension + 1);

      }

      KRATOS_CATCH( "" )
   }



   void SpatialLagrangianUwPFICElement::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant ; //= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
      double Permeability; double WaterBulk; double DeltaTime;
      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      Vector he = ZeroVector(dimension);
      double Factor = GetProperties()[STABILIZATION];
      for (unsigned int i = 0; i < dimension; i++)
      {
         he(i) = Factor*sqrt( 4.0*rIntegrationWeight/ sqrt(3.0) ); 
      }

      double Caux; 
      ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;
      GetValueOnIntegrationPoints( M_MODULUS, Mmodulus, CurrentProcessInfo);
      Caux = 1.0/Mmodulus[0];
      MatrixType Kh=rLeftHandSideMatrix;

      unsigned int indexpi = dimension;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            for ( unsigned int k = 0; k < dimension; ++k ) {
               rLeftHandSideMatrix(indexpi, indexpj) -= rVariables.N[i] * Caux*(1.0/2.0)*he(k) * rVariables.DN_DX(j,k) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
            }

            indexpj += (dimension + 1);
         }

         indexpi += (dimension + 1);
      }

      KRATOS_CATCH( "" )
   }



   //************************************************************************************
   //************************************************************************************

   void SpatialLagrangianUwPFICElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void SpatialLagrangianUwPFICElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }


   //****************************************************************************************
   //****************************************************************************************
   //   RESTOS DEL NAUFRAGIO RESTOS DEL NAUFRAGIO RESTOS DEL NAUFRAGIO RESTOS DEL NAUFRAGIO
   //****************************************************************************************
   //****************************************************************************************


}
