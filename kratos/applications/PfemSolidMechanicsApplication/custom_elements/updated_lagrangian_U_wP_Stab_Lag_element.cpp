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
#include "custom_elements/updated_lagrangian_U_wP_Stab_Lag_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application_variables.h"

// TO BE REMOVED. The differences (ie: constant lagrangian permeability) are now commented in the base element

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   UpdatedLagrangianUwPStabLagElement::UpdatedLagrangianUwPStabLagElement()
      : UpdatedLagrangianUwPStabElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPStabLagElement::UpdatedLagrangianUwPStabLagElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUwPStabElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPStabLagElement::UpdatedLagrangianUwPStabLagElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUwPStabElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   UpdatedLagrangianUwPStabLagElement::UpdatedLagrangianUwPStabLagElement( UpdatedLagrangianUwPStabLagElement const& rOther)
      :UpdatedLagrangianUwPStabElement(rOther)
       //,mDeterminantF0(rOther.mDeterminantF0)
       //,mDeformationGradientF0(rOther.mDeformationGradientF0)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   UpdatedLagrangianUwPStabLagElement&  UpdatedLagrangianUwPStabLagElement::operator=(UpdatedLagrangianUwPStabLagElement const& rOther)
   {
      UpdatedLagrangianUwPStabElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUwPStabLagElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new UpdatedLagrangianUwPStabLagElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer UpdatedLagrangianUwPStabLagElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      UpdatedLagrangianUwPStabLagElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new UpdatedLagrangianUwPStabLagElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   UpdatedLagrangianUwPStabLagElement::~UpdatedLagrangianUwPStabLagElement()
   {
   }



   // function to get the Eulerian permeability assuming a constant eulerian permeability
   //**************************************************************************************
   void UpdatedLagrangianUwPStabLagElement::GetPermeabilityTensor( const double& rPermeability, const Matrix& rF, Matrix& rPermeabilityTensor)
   {

      unsigned int thisSize = rF.size1();
      rPermeabilityTensor = ZeroMatrix( thisSize);
      for (unsigned int i = 0; i < thisSize; ++i)
      {
         rPermeabilityTensor(i,i) = rPermeability;
      }

      rPermeabilityTensor = prod(rPermeabilityTensor, trans( rF) );
      rPermeabilityTensor = prod( rF, rPermeabilityTensor);
      rPermeabilityTensor /= MathUtils<double>::Det(rF);

   }

   // Tangents to the Permeability Tensor with respect to the displacement
   // ********************************************************************
   double UpdatedLagrangianUwPStabLagElement::GetPermeabilityLDTerm( const Matrix& rPermeability, const Matrix& rF, const int i, const int j, const int k, const int l)
   {

      unsigned int thisSize = rF.size1();
      Matrix Identity = ZeroMatrix(thisSize);
      for (unsigned int h = 0; h < thisSize; h++)
         Identity(h,h) = 1.0;

      double LDTerm;
      LDTerm = Identity(i,k)*rPermeability(j,l) + Identity(i,l)*rPermeability(j,k);
      LDTerm = -rPermeability(i,j)*Identity(k,l);
      return LDTerm;

   }

   //************************************************************************************
   //************************************************************************************

   void UpdatedLagrangianUwPStabLagElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
      rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void UpdatedLagrangianUwPStabLagElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
      rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }







} //namespace kratos
