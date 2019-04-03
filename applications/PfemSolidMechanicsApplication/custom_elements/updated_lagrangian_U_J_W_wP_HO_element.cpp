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
#include "custom_elements/updated_lagrangian_U_J_W_wP_HO_element.hpp"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   UpdatedLagrangianUJWwPHOElement::UpdatedLagrangianUJWwPHOElement()
      : UpdatedLagrangianUJWwPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJWwPHOElement::UpdatedLagrangianUJWwPHOElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUJWwPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJWwPHOElement::UpdatedLagrangianUJWwPHOElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUJWwPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUJWwPHOElement::UpdatedLagrangianUJWwPHOElement( UpdatedLagrangianUJWwPHOElement const& rOther)
      :UpdatedLagrangianUJWwPElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUJWwPHOElement&  UpdatedLagrangianUJWwPHOElement::operator=(UpdatedLagrangianUJWwPHOElement const& rOther)
   {
      UpdatedLagrangianUJWwPElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJWwPHOElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUJWwPHOElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUJWwPHOElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUJWwPHOElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new UpdatedLagrangianUJWwPHOElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUJWwPHOElement::~UpdatedLagrangianUJWwPHOElement()
   {
   }

   // *********************************************************************************
   //         Calculate the Damping matrix part due to the high order terms
   void UpdatedLagrangianUJWwPHOElement::CalculateAndAddHighOrderDampingMatrix( MatrixType & rDampingMatrix, ElementDataType & rVariables, double & rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dofs_per_node = 2*dimension + 2;

      double WaterDensity =GetProperties().GetValue(DENSITY_WATER);
      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);

      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0; 

      Matrix ShapeFunction = ZeroMatrix(dimension,dimension*number_of_nodes);
      Matrix Divergence  = ZeroMatrix(1,dimension*number_of_nodes);
      Matrix NodalWaterVelocities = ZeroMatrix(dimension*number_of_nodes,1); 
      Matrix WaterVelocity = ZeroMatrix(dimension,1);

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < dimension; j++) {
            ShapeFunction(j,i*dimension+j)=rVariables.N[i];
            Divergence(0,i*dimension+j)=rVariables.DN_DX(i,j);
         }
      }

      Vector WaterVelocityAux = ZeroVector(3);
      for (unsigned int i = 0; i < number_of_nodes; i++){
         const array_1d<double, 3> &  rWaterVelocity = GetGeometry()[i].FastGetSolutionStepValue( WATER_VELOCITY );
         WaterVelocityAux += rVariables.N[i] * rWaterVelocity;
      }

      for (unsigned int j = 0; j < dimension; j++)
         WaterVelocity(j,0)=WaterVelocityAux(j);

      Matrix SmallMatrix = prod( Matrix(prod( trans(ShapeFunction),WaterVelocity)), Divergence);
      SmallMatrix *= WaterDensity*rIntegrationWeight/porosity;

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            for ( unsigned int iDim = 0; iDim < dimension; iDim++ ){
               for ( unsigned int jDim = 0; jDim < dimension; jDim++ ){
                  rDampingMatrix( i*dofs_per_node+iDim, j*dofs_per_node+jDim) += porosity*SmallMatrix(i*dimension+iDim,j*dimension+jDim);
                  rDampingMatrix( i*dofs_per_node+(dimension+1)+iDim, j*dofs_per_node+jDim) += SmallMatrix(i*dimension+iDim,j*dimension+jDim);
                  rDampingMatrix( i*dofs_per_node+iDim, j*dofs_per_node+(dimension+1)+jDim) += SmallMatrix(i*dimension+iDim,j*dimension+jDim);
                  rDampingMatrix( i*dofs_per_node+(dimension+1)+iDim, j*dofs_per_node+(dimension+1)+jDim) += SmallMatrix(i*dimension+iDim,j*dimension+jDim)/porosity;
               }
            }
         }
      }

      KRATOS_CATCH("")

   }

   //************************************************************************************
   //         Matrix due to the the water pressure stabilization contribution to the internal forces   
   void UpdatedLagrangianUJWwPHOElement::CalculateAndAddHighOrderKPP( MatrixType & rLeftHandSide, ElementDataType & rVariables, double & rIntegrationWeight)
   {

      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = 2*dimension + 2;

      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);
      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0; 
      const double & rWaterBulk = GetProperties()[WATER_BULK_MODULUS];

      Matrix ShapeFunction = ZeroMatrix(1,number_of_nodes);
      Matrix Gradient  = ZeroMatrix(dimension,number_of_nodes);
      Matrix NodalWaterVelocities = ZeroMatrix(dimension*number_of_nodes,1); 
      Matrix WaterVelocity = ZeroMatrix(1,dimension);

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         ShapeFunction(0,i)=rVariables.N[i];
         for (unsigned int j = 0; j < dimension; j++)
            Gradient(j,i)=rVariables.DN_DX(i,j);
      }

      Vector WaterVelocityAux = ZeroVector(3);
      for (unsigned int i = 0; i < number_of_nodes; i++){
         const array_1d<double, 3> &  rWaterVelocity = GetGeometry()[i].FastGetSolutionStepValue( WATER_VELOCITY );
         WaterVelocityAux += rVariables.N[i] * rWaterVelocity;
      }

      for (unsigned int j = 0; j < dimension; j++)
         WaterVelocity(0,j)=WaterVelocityAux(j);

      Matrix SmallMatrix = prod( Matrix(prod( trans(ShapeFunction),WaterVelocity)), Gradient);
      SmallMatrix *= rIntegrationWeight/porosity/rWaterBulk;    

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            rLeftHandSide( (i+1)*dofs_per_node-1, (j+1)*dofs_per_node-1) += SmallMatrix(i,j);
         }
      }

      KRATOS_CATCH("")
   }

   // ****************************************************************************
   //     Right hand side part due to the stabilization technique
   void UpdatedLagrangianUJWwPHOElement::CalculateAndAddHighOrderRHS(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int dofs_per_node = 2*dimension + 2;

      double porosity0 = GetProperties().GetValue( INITIAL_POROSITY);
      double porosity = 1.0 - (1.0-porosity0) / rVariables.detF0; 
      const double & rWaterBulk = GetProperties()[WATER_BULK_MODULUS];

      Matrix ShapeFunction = ZeroMatrix(1,number_of_nodes);
      Matrix Gradient  = ZeroMatrix(dimension,number_of_nodes);
      Matrix NodalWaterVelocities = ZeroMatrix(dimension*number_of_nodes,1); 
      Matrix WaterVelocity = ZeroMatrix(1,dimension);

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         ShapeFunction(0,i)=rVariables.N[i];
         for (unsigned int j = 0; j < dimension; j++)
            Gradient(j,i)=rVariables.DN_DX(i,j);
      }

      Vector WaterVelocityAux = ZeroVector(3);
      Vector WaterPressure = ZeroVector(number_of_nodes);
      for (unsigned int i = 0; i < number_of_nodes; i++){
         const array_1d<double, 3> &  rWaterVelocity = GetGeometry()[i].FastGetSolutionStepValue( WATER_VELOCITY );
         const double & rWaterPressure = GetGeometry()[i].FastGetSolutionStepValue(WATER_PRESSURE);
         WaterVelocityAux += rVariables.N[i] * rWaterVelocity;
         WaterPressure(i) = rWaterPressure;
      }

      for (unsigned int j = 0; j < dimension; j++)
         WaterVelocity(0,j)=WaterVelocityAux(j);

      Vector SmallRHS = prod(Matrix(prod( Matrix(prod( trans(ShapeFunction),WaterVelocity)), Gradient)),WaterPressure);
      SmallRHS *= rIntegrationWeight/porosity/rWaterBulk;

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         rRightHandSideVector((i+1)*dofs_per_node-1 ) -= SmallRHS(i);
      }

      KRATOS_CATCH("")

   }


   void UpdatedLagrangianUJWwPHOElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianUJWwPElement )
   }

   void UpdatedLagrangianUJWwPHOElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianUJWwPElement )
   }

}  // END KRATOS NAMESPACE
