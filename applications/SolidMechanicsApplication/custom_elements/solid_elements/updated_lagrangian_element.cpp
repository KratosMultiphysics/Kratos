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
#include "custom_elements/solid_elements/updated_lagrangian_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianElement::UpdatedLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianElement::UpdatedLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianElement::UpdatedLagrangianElement( UpdatedLagrangianElement const& rOther)
    :LargeDisplacementElement(rOther)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianElement&  UpdatedLagrangianElement::operator=(UpdatedLagrangianElement const& rOther)
{
    LargeDisplacementElement::operator=(rOther);

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

Element::Pointer UpdatedLagrangianElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared< UpdatedLagrangianElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianElement NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
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

    return Kratos::make_shared< UpdatedLagrangianElement >(NewElement);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianElement::~UpdatedLagrangianElement()
{
}



//************************************************************************************
//************************************************************************************


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void UpdatedLagrangianElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
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

    LargeDisplacementElement::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }


}


//************************************************************************************
//************************************************************************************

//**********************************GET DOUBLE VALUE**********************************
//************************************************************************************


void UpdatedLagrangianElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
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


void UpdatedLagrangianElement::Initialize()
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

void UpdatedLagrangianElement::InitializeElementData (ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    LargeDisplacementElement::InitializeElementData(rVariables,rCurrentProcessInfo);

    //Calculate Delta Position
    rVariables.DeltaPosition = this->CalculateDeltaPosition(rVariables.DeltaPosition);

    //set variables including all integration points values

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianElement::FinalizeStepVariables( ElementDataType & rVariables, const double& rPointNumber )
{
    //update internal (historical) variables
    mDeterminantF0[rPointNumber] = rVariables.detF * rVariables.detF0;
    noalias(mDeformationGradientF0[rPointNumber]) = prod(rVariables.F, rVariables.F0);
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void UpdatedLagrangianElement::CalculateKinematics(ElementDataType& rVariables,
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

    //
    //
    //Deformation Gradient F [dx_n+1/dx_n] to be updated
    noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    noalias(rVariables.DN_DX) = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Determinant of the Deformation Gradient F0
    //
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    const GeometryType& rGeometry = GetGeometry();
    ElementUtilities::CalculateLinearDeformationMatrix(rVariables.B,rGeometry,rVariables.DN_DX);

    KRATOS_CATCH( "" )
}

//*********************************COMPUTE KINETICS***********************************
//************************************************************************************

void UpdatedLagrangianElement::CalculateKinetics(ElementDataType& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    //TotalDeltaPosition must not be used in this element as mDeterminantF0 and mDeformationGradientF0 are stored for reduced order
    //however then the storage of variables in the full integration order quadrature must be considered

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    rVariables.DeltaPosition = this->CalculateTotalDeltaPosition(rVariables.DeltaPosition);

    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod, rVariables.DeltaPosition );

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], InvJ, rVariables.detJ);

    //Calculating the cartesian derivatives [dN/dx_n] = [dN/d£][d£/dx_0]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Deformation Gradient F [dx_n+1/dx_0] = [dx_n+1/d£] [d£/dx_0]
    noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //Determinant of the Deformation Gradient F0
    // (in this element F = F0, then F0 is set to the identity for coherence in the constitutive law)
    rVariables.detF0 = 1;
    rVariables.F0    = identity_matrix<double> ( dimension );

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    KRATOS_CATCH( "" )
}


//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void UpdatedLagrangianElement::CalculateDeformationGradient(Matrix& rF,
                                                            const Matrix& rDN_DX,
                                                            const Matrix& rDeltaPosition)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
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

      KRATOS_ERROR << " something is wrong with the dimension when computing F " << std::endl;

    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************


void UpdatedLagrangianElement::GetHistoricalVariables( ElementDataType& rVariables, const double& rPointNumber )
{
    LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

    //Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    noalias(rVariables.F0) = mDeformationGradientF0[rPointNumber];
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int UpdatedLagrangianElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementElement::Check(rCurrentProcessInfo);

    // Check compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
      if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
	correct_strain_measure = true;
    }

    if( correct_strain_measure == false )
      KRATOS_ERROR <<  "Large Displacement element with no Deformation Gradient strain measure" << std::endl;


    // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
    for(SizeType i=0; i<this->GetGeometry().size(); ++i)
      {
	// Nodal data
	Node<3> &rNode = this->GetGeometry()[i];
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rNode);
	//KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rNode);

	// Nodal dofs
	KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rNode);
	KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rNode);
	if( rCurrentProcessInfo[SPACE_DIMENSION] == 3)
	  KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rNode);
      }

    return ErrorCode;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************


void UpdatedLagrangianElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);

}

void UpdatedLagrangianElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);

}


} // Namespace Kratos
