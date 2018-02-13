//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                     PNavas $
//   Last modified by:    $Co-Author:               LMonforte $
//   Date:                $Date:                 October 2017 $
//   Revision:            $Revision:                     -0.1 $
//
//
//   Implementation of the Fluid Saturated porous media in a u-w Formulation

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_J_W_wP_stab_element.hpp"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   UpdatedLagrangianUJWwPStabElement::UpdatedLagrangianUJWwPStabElement()
      : UpdatedLagrangianUJWwPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJWwPStabElement::UpdatedLagrangianUJWwPStabElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUJWwPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJWwPStabElement::UpdatedLagrangianUJWwPStabElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUJWwPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUJWwPStabElement::UpdatedLagrangianUJWwPStabElement( UpdatedLagrangianUJWwPStabElement const& rOther)
      :UpdatedLagrangianUJWwPElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUJWwPStabElement&  UpdatedLagrangianUJWwPStabElement::operator=(UpdatedLagrangianUJWwPStabElement const& rOther)
   {
      UpdatedLagrangianUJWwPElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJWwPStabElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUJWwPStabElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJWwPStabElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUJWwPStabElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      NewElement.SetData(this->GetData());
      NewElement.SetFlags(this->GetFlags());

      return Element::Pointer( new UpdatedLagrangianUJWwPStabElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJWwPStabElement::~UpdatedLagrangianUJWwPStabElement()
   {
   }

   // *********************************************************************************
   //         Calculate the Mass matrix part due to the stabilization
   void UpdatedLagrangianUJWwPStabElement::CalculateAndAddMassStabilizationMatrix( MatrixType & rMassMatrix, ElementVariables & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dofs_per_node = 2*dimension + 2;

      double density_mixture0 = GetProperties()[DENSITY];
      double WaterDensity =GetProperties().GetValue(DENSITY_WATER);
      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);

      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0; 
      double density_solid = (density_mixture0 - porosity0*WaterDensity) / ( 1.0 - porosity0);
      double CurrentDensity = ( 1.0 - porosity) * density_solid + porosity * WaterDensity;

      const double &  rWaterBulk = GetProperties()[WATER_BULK_MODULUS];

      const double & rStabFactor = GetProperties().GetValue( STABILIZATION_FACTOR_WP);

      Matrix SmallMatrix = ZeroMatrix( number_of_nodes, number_of_nodes);

      if ( ( fabs(rStabFactor) > 1.0e-9) && dimension==2)  {

         // assume that K is a unit tensor by k
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               SmallMatrix(i,j) -= rVariables.N[i] * rVariables.N[j] * rStabFactor * CurrentDensity * rIntegrationWeight / rWaterBulk;
            }
         } 
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               rMassMatrix( (i+1)*dofs_per_node-1, (j+1)*dofs_per_node-1) += SmallMatrix(i,j);
            }
         }
      }

   KRATOS_CATCH("")

   }

   // *********************************************************************************
   //         Calculate the Damping matrix part due to the stabilization
   void UpdatedLagrangianUJWwPStabElement::CalculateAndAddDampingStabilizationMatrix( MatrixType & rDampingMatrix, ElementVariables & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      KRATOS_CATCH("")
   }

   //************************************************************************************
   //         Matrix due to the the water pressure stabilization contribution to the internal forces   
   void UpdatedLagrangianUJWwPStabElement::CalculateAndAddKPPStab( MatrixType & rLeftHandSide, ElementVariables & rVariables, double & rIntegrationWeight)
   {

      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = 2*dimension + 2;

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

      const double & rWaterBulk = GetProperties()[WATER_BULK_MODULUS];
      double ConstrainedModulusMixture = 1 + ConstrainedModulus / rWaterBulk;

      //double StabFactor = CalculateStabilizationFactor( Variables, StabFactor); 
      double StabFactor = GetProperties().GetValue( STABILIZATION_FACTOR_WP);

      Matrix SmallMatrix = ZeroMatrix( number_of_nodes, number_of_nodes);

      if ( ( fabs(StabFactor) > 1.0e-9) && dimension==2)  {

         // assume that K is a unit tensor by k
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               for (unsigned p = 0; p < dimension; p++) {
                  SmallMatrix(i,j) -= rVariables.DN_DX(i,p) * rVariables.DN_DX(j,p) * StabFactor * ConstrainedModulusMixture * rIntegrationWeight;
               }
            }
         } 
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               rLeftHandSide( (i+1)*dofs_per_node-1, (j+1)*dofs_per_node-1) += SmallMatrix(i,j);
            }
         }
      }

      KRATOS_CATCH("")
   }
 
   // ****************************************************************************
   //     Right hand side part due to the stabilization technique
   void UpdatedLagrangianUJWwPStabElement::CalculateAndAddStabilizationRHS(VectorType& rRightHandSideVector,
         ElementVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = 2*dimension + 2;

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

      const double  rWaterBulk = GetProperties()[WATER_BULK_MODULUS];
      double ConstrainedModulusMixture = 1 + ConstrainedModulus / rWaterBulk;

      //double StabFactor = CalculateStabilizationFactor( Variables, StabFactor); 
      double StabFactor = GetProperties().GetValue( STABILIZATION_FACTOR_WP); 

      Vector SmallRHS = ZeroVector(number_of_nodes);

      if ( ( fabs(StabFactor) > 1.0e-9) && dimension==2)  {
         
         for (unsigned int i = 0; i < number_of_nodes; i++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
               for (unsigned p = 0; p < dimension; p++) {
		  const double & WaterPressure = GetGeometry()[i].FastGetSolutionStepValue(WATER_PRESSURE);
                  SmallRHS(i) += rVariables.DN_DX(i,p) * rVariables.DN_DX(j,p) * StabFactor * ConstrainedModulusMixture * WaterPressure * rIntegrationWeight;
               }
            }
         }
      	for (unsigned int i = 0; i < number_of_nodes; i++) {
            rRightHandSideVector((i+1)*dofs_per_node-1 ) += SmallRHS(i);
      	}
      } 

      KRATOS_CATCH("")

   }


   void UpdatedLagrangianUJWwPStabElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianUJWwPElement )
   }

   void UpdatedLagrangianUJWwPStabElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianUJWwPElement )
   }

}  // END KRATOS NAMESPACE
