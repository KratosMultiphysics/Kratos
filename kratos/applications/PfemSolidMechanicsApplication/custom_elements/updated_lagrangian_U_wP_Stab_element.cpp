//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_U_wP_Stab_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   UpdatedLagrangianUwPStabElement::UpdatedLagrangianUwPStabElement()
      : UpdatedLagrangianUwPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPStabElement::UpdatedLagrangianUwPStabElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUwPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPStabElement::UpdatedLagrangianUwPStabElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUwPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUwPStabElement::UpdatedLagrangianUwPStabElement( UpdatedLagrangianUwPStabElement const& rOther)
      :UpdatedLagrangianUwPElement(rOther)
       //,mDeterminantF0(rOther.mDeterminantF0)
       //,mDeformationGradientF0(rOther.mDeformationGradientF0)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUwPStabElement&  UpdatedLagrangianUwPStabElement::operator=(UpdatedLagrangianUwPStabElement const& rOther)
   {
      UpdatedLagrangianUwPElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUwPStabElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUwPStabElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUwPStabElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUwPStabElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

      //-----------//

      NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


      if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
         NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

         if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
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

      return Element::Pointer( new UpdatedLagrangianUwPStabElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPStabElement::~UpdatedLagrangianUwPStabElement()
   {
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUwPStabElement::InitializeGeneralVariables (GeneralVariables & rVariables, 
								     const ProcessInfo& rCurrentProcessInfo)
   {
     KRATOS_TRY

     UpdatedLagrangianUwPElement::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

     //stabilization factor
     double StabilizationFactor = 1.0;
     if( GetProperties().Has(STABILIZATION_FACTOR) ){
       StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
     }
     else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR) ){
       StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR];
     }
     GetProperties().SetValue(STABILIZATION_FACTOR, StabilizationFactor);

     KRATOS_CATCH( "" )
   }

   //********* STABILIZATION TERM *******************************************************
   //****************** Fluid Pressure Laplacian ****************************************

   void UpdatedLagrangianUwPStabElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant ;
      double Permeability; double WaterBulk; double DeltaTime;
      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      double StabilizationAlpha, Caux, StabilizationFactor; 

      ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;
      GetValueOnIntegrationPoints(M_MODULUS, Mmodulus, CurrentProcessInfo);
      Caux = 1.0/Mmodulus[0];

      double he;
      he = GetElementSize( rVariables.DN_DX);
      StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];

      StabilizationAlpha = pow(he, 2.0) * Caux / ( 6.0) - DeltaTime*Permeability; 
      StabilizationAlpha *= StabilizationFactor;

      if (StabilizationAlpha < 0.0)
      {
         return;
      }

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;


      Matrix K = ZeroMatrix(dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = 1.0;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {


            const double& CurrentPressure   = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE );
            const double& PreviousPressure = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE , 1);
            double DeltaPressure = CurrentPressure - PreviousPressure;

            for ( unsigned int p = 0; p < dimension; ++p )
            {

               for ( unsigned int q = 0; q < dimension; ++q )
               {

                  rRightHandSideVector[indexp] += StabilizationAlpha * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q)* DeltaPressure *rIntegrationWeight * ScalingConstant / rVariables.detF0;

               }


            }

         }

         indexp += (dimension + 1);

      }


      // std::cout<<std::endl;
      // std::cout<<" auxiliar " <<auxiliar<<" F0 "<<rVariables.detF0<<std::endl;
      // std::cout<<" Fpres STABILI "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )

   }


   //********* TANGENT MATRIX  STABILIZATION TERM ***************************************
   //****************** Fluid Pressure Laplacian ****************************************

   void UpdatedLagrangianUwPStabElement::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;
      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      double StabilizationAlpha, Caux, StabilizationFactor; 

      ProcessInfo CurrentProcessInfo;
      std::vector<double> Mmodulus;
      GetValueOnIntegrationPoints(M_MODULUS, Mmodulus, CurrentProcessInfo);
      Caux = 1.0/Mmodulus[0];

      double he;
      he = GetElementSize( rVariables.DN_DX);
    	StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];

      StabilizationAlpha = pow(he, 2.0) * Caux / ( 6.0) - DeltaTime*Permeability; 
      StabilizationAlpha *= StabilizationFactor;

      if (StabilizationAlpha < 0.0)
      {
         return;
      }

      Matrix K = ZeroMatrix(dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = 1.0;

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            for ( unsigned int p = 0; p < dimension; ++p ) {
               for ( unsigned int q = 0; q < dimension; ++q ) {
                  rLeftHandSideMatrix(indexpi, indexpj) -= StabilizationAlpha * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q) * rIntegrationWeight* ScalingConstant / rVariables.detF0;
               }
            }

            indexpj += (dimension + 1);
         }

         indexpi += (dimension + 1);
      }

      // std::cout<<std::endl;
      // std::cout<<" Kpp STABILI "<< (rLeftHandSideMatrix-Kh) <<std::endl;
      //std::cout<<" Kpp "<< (rLeftHandSideMatrix-Kh) / ScalingConstant / rIntegrationWeight* rVariables.detF0  * WaterBulk <<std::endl;


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUwPStabElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void UpdatedLagrangianUwPStabElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }

}
