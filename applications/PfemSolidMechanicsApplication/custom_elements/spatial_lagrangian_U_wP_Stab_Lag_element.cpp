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
#include "custom_elements/spatial_lagrangian_U_wP_Stab_Lag_element.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   SpatialLagrangianUwPStabLagElement::SpatialLagrangianUwPStabLagElement()
      : SpatialLagrangianUwPStabElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPStabLagElement::SpatialLagrangianUwPStabLagElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : SpatialLagrangianUwPStabElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPStabLagElement::SpatialLagrangianUwPStabLagElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : SpatialLagrangianUwPStabElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   SpatialLagrangianUwPStabLagElement::SpatialLagrangianUwPStabLagElement( SpatialLagrangianUwPStabLagElement const& rOther)
      :SpatialLagrangianUwPStabElement(rOther)
       //,mDeterminantF0(rOther.mDeterminantF0)
       //,mDeformationGradientF0(rOther.mDeformationGradientF0)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   SpatialLagrangianUwPStabLagElement&  SpatialLagrangianUwPStabLagElement::operator=(SpatialLagrangianUwPStabLagElement const& rOther)
   {
      SpatialLagrangianUwPStabElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer SpatialLagrangianUwPStabLagElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new SpatialLagrangianUwPStabLagElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer SpatialLagrangianUwPStabLagElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      SpatialLagrangianUwPStabLagElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new SpatialLagrangianUwPStabLagElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   SpatialLagrangianUwPStabLagElement::~SpatialLagrangianUwPStabLagElement()
   {
   }



   // function to get the Eulerian permeability assuming a constant eulerian permeability
   //**************************************************************************************
   void SpatialLagrangianUwPStabLagElement::GetPermeabilityTensor( const double& rPermeability, const Matrix& rF, Matrix& rPermeabilityTensor)
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
   double SpatialLagrangianUwPStabLagElement::GetPermeabilityLDTerm( const Matrix& rPermeability, const Matrix& rF, const int i, const int j, const int k, const int l)
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

   void SpatialLagrangianUwPStabLagElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
      rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void SpatialLagrangianUwPStabLagElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
      rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }







} //namespace kratos
