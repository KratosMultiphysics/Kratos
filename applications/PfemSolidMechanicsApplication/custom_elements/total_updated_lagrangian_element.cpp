//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/total_updated_lagrangian_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TotalUpdatedLagrangianElement::TotalUpdatedLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangianElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TotalUpdatedLagrangianElement::TotalUpdatedLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangianElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TotalUpdatedLagrangianElement::TotalUpdatedLagrangianElement( TotalUpdatedLagrangianElement const& rOther)
    :UpdatedLagrangianElement(rOther)
{
}



//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

TotalUpdatedLagrangianElement&  TotalUpdatedLagrangianElement::operator=(TotalUpdatedLagrangianElement const& rOther)
{
    UpdatedLagrangianElement::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer TotalUpdatedLagrangianElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new TotalUpdatedLagrangianElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer TotalUpdatedLagrangianElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    TotalUpdatedLagrangianElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );


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

    return Element::Pointer( new TotalUpdatedLagrangianElement(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TotalUpdatedLagrangianElement::~TotalUpdatedLagrangianElement()
{
}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


void TotalUpdatedLagrangianElement::InitializeElementData (ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    LargeDisplacementElement::InitializeElementData(rVariables,rCurrentProcessInfo);

    //Calculate Delta Position
    ElementUtilities::CalculateDeltaPosition(rVariables.DeltaPosition,this->GetGeometry());

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

}



//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void TotalUpdatedLagrangianElement::CalculateKinematics(ElementDataType& rVariables,
        const double& rPointNumber)

{
    KRATOS_TRY

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Compute cartesian derivatives [dN/dx_n]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    //
    //
    //Deformation Gradient F [dx_n+1/dx_n] to be updated
    noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //
    //
    //

    //
    //

    //Determinant of the Deformation Gradient F0
    //
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void TotalUpdatedLagrangianElement::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 0 );
            rB( 1, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 1 ) + rF( 0, 1 ) * rDN_DX( i, 0 );
            rB( 2, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 1 ) + rF( 1, 1 ) * rDN_DX( i, 0 );

        }
    }
    else if( dimension == 3 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 2 ) = rF( 2, 0 ) * rDN_DX( i, 0 );
            rB( 1, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 2 ) = rF( 2, 1 ) * rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rF( 0, 2 ) * rDN_DX( i, 2 );
            rB( 2, index + 1 ) = rF( 1, 2 ) * rDN_DX( i, 2 );
            rB( 2, index + 2 ) = rF( 2, 2 ) * rDN_DX( i, 2 );
            rB( 3, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 1 ) + rF( 0, 1 ) * rDN_DX( i, 0 );
            rB( 3, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 1 ) + rF( 1, 1 ) * rDN_DX( i, 0 );
            rB( 3, index + 2 ) = rF( 2, 0 ) * rDN_DX( i, 1 ) + rF( 2, 1 ) * rDN_DX( i, 0 );
            rB( 4, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 2 ) + rF( 0, 2 ) * rDN_DX( i, 1 );
            rB( 4, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 2 ) + rF( 1, 2 ) * rDN_DX( i, 1 );
            rB( 4, index + 2 ) = rF( 2, 1 ) * rDN_DX( i, 2 ) + rF( 2, 2 ) * rDN_DX( i, 1 );
            rB( 5, index + 0 ) = rF( 0, 2 ) * rDN_DX( i, 0 ) + rF( 0, 0 ) * rDN_DX( i, 2 );
            rB( 5, index + 1 ) = rF( 1, 2 ) * rDN_DX( i, 0 ) + rF( 1, 0 ) * rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rF( 2, 2 ) * rDN_DX( i, 0 ) + rF( 2, 0 ) * rDN_DX( i, 2 );

        }
    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void TotalUpdatedLagrangianElement::TransformElementData(ElementDataType& rVariables, const double& rPointNumber)
{
  // pull_back the stresses to last_known configuration
  mConstitutiveLawVector[rPointNumber]->TransformStresses(rVariables.StressVector, rVariables.F, rVariables.detF, ConstitutiveLaw::StressMeasure_Cauchy, ConstitutiveLaw::StressMeasure_PK2);

  // pull_back the constitutive tensor to last_known configuration
  mConstitutiveLawVector[rPointNumber]->PullBackConstitutiveMatrix(rVariables.ConstitutiveMatrix, rVariables.F);
}

//************************************************************************************
//************************************************************************************

void TotalUpdatedLagrangianElement::GetHistoricalVariables( ElementDataType& rVariables, const double& rPointNumber )
{
    LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

    //Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& TotalUpdatedLagrangianElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void TotalUpdatedLagrangianElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianElement )

}

void TotalUpdatedLagrangianElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianElement )

}


} // Namespace Kratos
