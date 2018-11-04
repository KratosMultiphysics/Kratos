//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/updated_lagrangian_U_P_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUPElement::UpdatedLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementUPElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUPElement::UpdatedLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementUPElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianUPElement::UpdatedLagrangianUPElement( UpdatedLagrangianUPElement const& rOther)
    :LargeDisplacementUPElement(rOther)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianUPElement&  UpdatedLagrangianUPElement::operator=(UpdatedLagrangianUPElement const& rOther)
{
    LargeDisplacementUPElement::operator=(rOther);

    mDeformationGradientF0.clear();
    mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
    {
        mDeformationGradientF0[i] = rOther.mDeformationGradientF0[i];
    }

    mDeterminantF0 = rOther.mDeterminantF0;

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared< UpdatedLagrangianUPElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianUPElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

    return Kratos::make_shared< UpdatedLagrangianUPElement >(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUPElement::~UpdatedLagrangianUPElement()
{
}

//************************************************************************************
//************************************************************************************


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void UpdatedLagrangianUPElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();


    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
      {
	mDeterminantF0[PointNumber] = rValues[PointNumber];
	mConstitutiveLawVector[PointNumber]->SetValue(rVariable, rValues[PointNumber], rCurrentProcessInfo);
      }

  }
  else{

    LargeDisplacementUPElement::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }


}


//************************************************************************************
//************************************************************************************

//**********************************GET DOUBLE VALUE**********************************
//************************************************************************************


void UpdatedLagrangianUPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
      rValues.resize( integration_points_number );

    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
      {
	rValues[PointNumber] = mDeterminantF0[PointNumber];
      }

  }
  else{

    LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }

}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPElement::Initialize()
{
    KRATOS_TRY

    LargeDisplacementElement::Initialize();

    SizeType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    //Resize historic deformation gradient
    if ( mDeformationGradientF0.size() != integration_points_number )
        mDeformationGradientF0.resize( integration_points_number );

    if ( mDeterminantF0.size() != integration_points_number )
        mDeterminantF0.resize( integration_points_number, false );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
    {
        mDeterminantF0[PointNumber] = 1;
        mDeformationGradientF0[PointNumber] = identity_matrix<double> (dimension);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPElement::InitializeElementData (ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    LargeDisplacementUPElement::InitializeElementData(rVariables,rCurrentProcessInfo);

    //Calculate Delta Position
    ElementUtilities::CalculateDeltaPosition(rVariables.DeltaPosition,this->GetGeometry());

    //set variables including all integration points values

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

}


////************************************************************************************
////************************************************************************************

void UpdatedLagrangianUPElement::FinalizeStepVariables( ElementDataType & rVariables, const double& rPointNumber )
{
    //update internal (historical) variables
    mDeterminantF0[rPointNumber] = rVariables.detF * rVariables.detF0;
    noalias(mDeformationGradientF0[rPointNumber]) = prod(rVariables.F, rVariables.F0);
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
{

    //contributions to stiffness matrix calculated on the reference config
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    LargeDisplacementUPElement::CalculateAndAddLHS( rLocalSystem, rVariables, rIntegrationWeight );

    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    //KRATOS_WATCH( rLeftHandSideMatrix )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    //contribution to external forces
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    LargeDisplacementUPElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, rIntegrationWeight );

    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    //KRATOS_WATCH( rRightHandSideVector )
}



//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void UpdatedLagrangianUPElement::CalculateKinematics(ElementDataType& rVariables,
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
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Current Deformation Gradient [dx_n+1/dx_n]
    //this->CalculateDeformationGradient (rVariables.F, rVariables.DN_DX, rVariables.DeltaPosition);

    //Deformation Gradient F [dx_n+1/dx_n] to be updated
    noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    noalias(rVariables.DN_DX) = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Determinant of the Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];


    //Compute the deformation matrix B
    const GeometryType& rGeometry = GetGeometry();
    ElementUtilities::CalculateLinearDeformationMatrix(rVariables.B,rGeometry,rVariables.DN_DX);

    KRATOS_CATCH( "" )
}


//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void UpdatedLagrangianUPElement::CalculateDeformationGradient(Matrix& rF,
                                                              const Matrix& rDN_DX,
                                                              const Matrix& rDeltaPosition)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    rF = identity_matrix<double> ( dimension );

    if( dimension == 2 )
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
        }

    }
    else if( dimension == 3)
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {

            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 0 , 2 ) += rDeltaPosition(i,0)*rDN_DX ( i , 2 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
            rF ( 1 , 2 ) += rDeltaPosition(i,1)*rDN_DX ( i , 2 );
            rF ( 2 , 0 ) += rDeltaPosition(i,2)*rDN_DX ( i , 0 );
            rF ( 2 , 1 ) += rDeltaPosition(i,2)*rDN_DX ( i , 1 );
            rF ( 2 , 2 ) += rDeltaPosition(i,2)*rDN_DX ( i , 2 );
        }

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPElement::GetHistoricalVariables( ElementDataType& rVariables, const double& rPointNumber )
{
    LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

    //Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    noalias(rVariables.F0) = mDeformationGradientF0[rPointNumber];
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianUPElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************


void UpdatedLagrangianUPElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementUPElement )
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
}

void UpdatedLagrangianUPElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementUPElement )
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}


} // Namespace Kratos
