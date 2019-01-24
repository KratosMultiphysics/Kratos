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
#include "custom_elements/axisym_updated_lagrangian_U_wP_Stab_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"
#include "custom_utilities/axisym_water_pressure_utilities.hpp"
namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   AxisymUpdatedLagrangianUwPStabElement::AxisymUpdatedLagrangianUwPStabElement()
      : AxisymUpdatedLagrangianUwPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUwPStabElement::AxisymUpdatedLagrangianUwPStabElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : AxisymUpdatedLagrangianUwPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUwPStabElement::AxisymUpdatedLagrangianUwPStabElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : AxisymUpdatedLagrangianUwPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUwPStabElement::AxisymUpdatedLagrangianUwPStabElement( AxisymUpdatedLagrangianUwPStabElement const& rOther)
      :AxisymUpdatedLagrangianUwPElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   AxisymUpdatedLagrangianUwPStabElement&  AxisymUpdatedLagrangianUwPStabElement::operator=(AxisymUpdatedLagrangianUwPStabElement const& rOther)
   {
      AxisymUpdatedLagrangianUwPElement::operator=(rOther);

      mTimeStep = rOther.mTimeStep;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUwPStabElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new AxisymUpdatedLagrangianUwPStabElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUwPStabElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      AxisymUpdatedLagrangianUwPStabElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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
      NewElement.mTimeStep = mTimeStep; 

      NewElement.SetData(this->GetData());
      NewElement.SetFlags(this->GetFlags());

      return Element::Pointer( new AxisymUpdatedLagrangianUwPStabElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUwPStabElement::~AxisymUpdatedLagrangianUwPStabElement()
   {
   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUwPStabElement::InitializeElementData (ElementDataType & rVariables, 
         const ProcessInfo& rCurrentProcessInfo)
   {
      KRATOS_TRY

      AxisymUpdatedLagrangianUwPElement::InitializeElementData(rVariables,rCurrentProcessInfo);

      //stabilization factor
      double StabilizationFactor = 1.0;
      if( GetProperties().Has(STABILIZATION_FACTOR_WP) ){
         StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_WP];
      }
      else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR_WP) ){
         StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR_WP];
      }
      GetProperties().SetValue(STABILIZATION_FACTOR_WP, StabilizationFactor);

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUwPStabElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      KRATOS_TRY

      // define some variables
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;


      // ComputeBaseClass LHS
      double IntegrationWeight = rIntegrationWeight; 
      AxisymUpdatedLagrangianUwPElement::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight);

      IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      // Add stabilization
      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
      AxisymWaterPressureUtilities WaterUtility; 
      int number_of_variables = dimension + 1; 
      ProcessInfo SomeProcessInfo; 
      std::vector< double> Mmodulus;
      GetValueOnIntegrationPoints( M_MODULUS, Mmodulus, SomeProcessInfo);
      double ConstrainedModulus = Mmodulus[0];
      if ( ConstrainedModulus < 1e-5)
      {
         const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
         const double& nu    = GetProperties()[POISSON_RATIO];
         ConstrainedModulus =  YoungModulus * ( 1.0-nu)/(1.0+nu) / (1.0-2.0*nu);
      }
      // 1. Create (make pointers) variables
      WaterPressureUtilities::HydroMechanicalVariables HMVariables(GetGeometry(), GetProperties() );

      HMVariables.SetBMatrix( rVariables.B);
      HMVariables.SetShapeFunctionsDerivatives( rVariables.DN_DX);
      HMVariables.SetShapeFunctions( rVariables.N);

      HMVariables.DeltaTime = mTimeStep;
      HMVariables.detF0 = rVariables.detF0;
      HMVariables.CurrentRadius = rVariables.CurrentRadius;
      HMVariables.ConstrainedModulus = ConstrainedModulus;
      HMVariables.number_of_variables = number_of_variables;

      rLeftHandSideMatrix = WaterUtility.CalculateAndAddStabilizationLHS( HMVariables, rLeftHandSideMatrix, IntegrationWeight);

      rVariables.detF = DeterminantF; 
      rVariables.detF0 /= rVariables.detF; 

      KRATOS_CATCH("")
   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUwPStabElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {
      KRATOS_TRY

      // define some variables
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      rVariables.detF0 *= rVariables.detF;
      double DeterminantF = rVariables.detF;
      rVariables.detF = 1.0;

      // Compute Base Class RHS
      double IntegrationWeight = rIntegrationWeight; 
      AxisymUpdatedLagrangianUwPElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, rIntegrationWeight);

      // Add stabilization
      Vector& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
      AxisymWaterPressureUtilities WaterUtility; 
      int number_of_variables = dimension + 1; 
      IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];
      ProcessInfo SomeProcessInfo; 
      std::vector< double> Mmodulus;
      GetValueOnIntegrationPoints( M_MODULUS, Mmodulus, SomeProcessInfo);
      double ConstrainedModulus = Mmodulus[0];
      if ( ConstrainedModulus < 1e-5)
      {
         const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
         const double& nu    = GetProperties()[POISSON_RATIO];
         ConstrainedModulus =  YoungModulus * ( 1.0-nu)/(1.0+nu) / (1.0-2.0*nu);
      }

      // 1. Create (make pointers) variables
      WaterPressureUtilities::HydroMechanicalVariables HMVariables(GetGeometry(), GetProperties() );

      HMVariables.SetBMatrix( rVariables.B);
      HMVariables.SetShapeFunctionsDerivatives( rVariables.DN_DX);
      HMVariables.SetShapeFunctions( rVariables.N);

      HMVariables.DeltaTime = mTimeStep;
      HMVariables.detF0 = rVariables.detF0;
      HMVariables.CurrentRadius = rVariables.CurrentRadius;
      HMVariables.ConstrainedModulus = ConstrainedModulus;
      HMVariables.number_of_variables = number_of_variables;

      rRightHandSideVector = WaterUtility.CalculateAndAddStabilization( HMVariables, rRightHandSideVector, IntegrationWeight);


      rVariables.detF = DeterminantF;
      rVariables.detF0 /= rVariables.detF;

      KRATOS_CATCH("")
   }



   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUwPStabElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUwPElement )
   }

   void AxisymUpdatedLagrangianUwPStabElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AxisymUpdatedLagrangianUwPElement )
   }

}  // END KRATOS NAMESPACE
