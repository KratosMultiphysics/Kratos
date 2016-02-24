

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_updated_lagrangian_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymUpdatedLagrangianElement::AxisymUpdatedLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymUpdatedLagrangianElement::AxisymUpdatedLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementElement( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    //mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    //mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymUpdatedLagrangianElement::AxisymUpdatedLagrangianElement( AxisymUpdatedLagrangianElement const& rOther)
    :LargeDisplacementElement(rOther)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

AxisymUpdatedLagrangianElement&  AxisymUpdatedLagrangianElement::operator=(AxisymUpdatedLagrangianElement const& rOther)
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

Element::Pointer AxisymUpdatedLagrangianElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new AxisymUpdatedLagrangianElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer AxisymUpdatedLagrangianElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    AxisymUpdatedLagrangianElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

        
    return Element::Pointer( new AxisymUpdatedLagrangianElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

AxisymUpdatedLagrangianElement::~AxisymUpdatedLagrangianElement()
{
}

//************************************************************************************
//************************************************************************************


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
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


void AxisymUpdatedLagrangianElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
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

void AxisymUpdatedLagrangianElement::Initialize()
{
    KRATOS_TRY

    LargeDisplacementElement::Initialize();

    SizeType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    //Resize historic deformation gradient
    if ( mDeformationGradientF0.size() != integration_points_number )
        mDeformationGradientF0.resize( integration_points_number );

    if ( mDeterminantF0.size() != integration_points_number )
        mDeterminantF0.resize( integration_points_number, false );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
    {
        mDeterminantF0[PointNumber] = 1;
        mDeformationGradientF0[PointNumber] = identity_matrix<double> (3);
    }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();

    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    rVariables.B.resize( 4 , number_of_nodes * 2 );

    rVariables.F.resize( 3, 3 );

    rVariables.F0.resize( 3, 3 );

    rVariables.FT.resize( 3, 3 );

    rVariables.ConstitutiveMatrix.resize( 4, 4 );

    rVariables.StrainVector.resize( 4 );

    rVariables.StressVector.resize( 4 );

    rVariables.DN_DX.resize( number_of_nodes, 2 );

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );


    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    
}

////************************************************************************************
////************************************************************************************

void AxisymUpdatedLagrangianElement::FinalizeStepVariables( GeneralVariables & rVariables, const double& rPointNumber )
{ 
    //update internal (historical) variables
    mDeterminantF0[rPointNumber]         = rVariables.detF * rVariables.detF0;
    mDeformationGradientF0[rPointNumber] = prod(rVariables.F, rVariables.F0);
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

    //contributions to stiffness matrix calculated on the reference config

    LargeDisplacementElement::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )
}


//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

    //contribution to external forces

    LargeDisplacementElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, IntegrationWeight );

    //KRATOS_WATCH( rRightHandSideVector )
}



//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& AxisymUpdatedLagrangianElement::CalculateTotalMass( double& rTotalMass, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Compute the Volume Change acumulated:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    rTotalMass = 0;
    //reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	//compute element kinematics
	this->CalculateKinematics(Variables,PointNumber);
	
	//getting informations for integration
        double IntegrationWeight = Variables.detJ * integration_points[PointNumber].Weight();

	//compute point volume change
	double PointVolumeChange = 0;
	PointVolumeChange = this->CalculateVolumeChange( PointVolumeChange, Variables );
	
	rTotalMass += PointVolumeChange * GetProperties()[DENSITY] * 2.0 * 3.141592654 * Variables.CurrentRadius * IntegrationWeight;

      }

    return rTotalMass;

    KRATOS_CATCH( "" )
}



//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void AxisymUpdatedLagrangianElement::CalculateKinematics(GeneralVariables& rVariables,
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
    rVariables.N=row( Ncontainer, rPointNumber);

    //Calculate IntegrationPoint radius
    CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius, rVariables.N);

    //Current Deformation Gradient [dx_n+1/dx_n]
    CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition, rVariables.CurrentRadius, rVariables.ReferenceRadius);

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Determinant of the Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    //Compute the deformation matrix B
    CalculateDeformationMatrix(rVariables.B, rVariables.DN_DX, rVariables.N, rVariables.CurrentRadius);


    KRATOS_CATCH( "" )
}



//*************************COMPUTE AXYSIMMETRIC RADIUS********************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::CalculateRadius(double & rCurrentRadius,
        double & rReferenceRadius,
        const Vector& rN)


{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rCurrentRadius=0;
    rReferenceRadius=0;

    if ( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
	    array_1d<double, 3 > & CurrentPosition      = GetGeometry()[i].Coordinates();
	    array_1d<double, 3 > ReferencePosition      = CurrentPosition - DeltaDisplacement;

            rCurrentRadius   += CurrentPosition[0]*rN[i];
            rReferenceRadius += ReferencePosition[0]*rN[i];
            //std::cout<<" node "<<i<<" -> DeltaDisplacement : "<<DeltaDisplacement<<std::endl;
        }
    }


    if ( dimension == 3 )
    {
        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::CalculateDeformationGradient(const Matrix& rDN_DX,
        Matrix&  rF,
        Matrix&  rDeltaPosition,
        double & rCurrentRadius,
        double & rReferenceRadius)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rF = identity_matrix<double> ( 3 );

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
        }

        rF ( 2 , 2 ) = rCurrentRadius/rReferenceRadius;
    }
    else if( dimension == 3)
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************


void AxisymUpdatedLagrangianElement::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rDN_DX,
        Vector& rN,
        double & rCurrentRadius)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rB.clear(); //set all components to zero

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rN[i]/rCurrentRadius;
            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

        }

    }
    else if( dimension == 3 )
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::CalculateGreenLagrangeStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    //Right Cauchy-Green Calculation
    Matrix C ( 3, 3 );
    noalias( C ) = prod( trans( rF ), rF );

    if( dimension == 2 )
    {

        //Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 4 ) rStrainVector.resize( 4, false );

        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        rStrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

        rStrainVector[3] = C( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::CalculateAlmansiStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( 3, 3 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    if( dimension == 2 )
    {

        //Almansi Strain Calculation
        if ( rStrainVector.size() != 4 ) rStrainVector.resize( 4, false );

        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = 0.5 * ( 1.00 - InverseLeftCauchyGreen( 2, 2 ) );

        rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }


    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::CalculateAndAddKuug(MatrixType& rK,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();

    //Matrix Kh = rK;

    // axisymmetric geometric matrix

    double alpha1 = 0;
    double alpha2 = 0;
    double alpha3 = 0;

    unsigned int indexi = 0;
    unsigned int indexj = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        indexj =0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            alpha1 = rVariables.DN_DX(j,0) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[0] + rVariables.DN_DX(i,1) * rVariables.StressVector[3] );
            alpha2 = rVariables.DN_DX(j,1) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[3] + rVariables.DN_DX(i,1) * rVariables.StressVector[1] );
            alpha3 = rVariables.N[i] * rVariables.N[j] * rVariables.StressVector[2] * (1.0/rVariables.CurrentRadius*rVariables.CurrentRadius);

            rK(indexi,indexj)     += (alpha1 + alpha2 + alpha3) * rIntegrationWeight ;
            rK(indexi+1,indexj+1) += (alpha1 + alpha2) * rIntegrationWeight ;

            indexj+=2;
        }

        indexi+=2;

    }

    //std::cout<<std::endl;
    //std::cout<<" Kgeo "<<rK-Kh<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianElement::GetHistoricalVariables( GeneralVariables& rVariables, const double& rPointNumber )
{
    LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

    //Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    rVariables.CurrentRadius = rVariables.ReferenceRadius;
}


//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& AxisymUpdatedLagrangianElement::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY
      
    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void AxisymUpdatedLagrangianElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
}

void AxisymUpdatedLagrangianElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}


} // Namespace Kratos


