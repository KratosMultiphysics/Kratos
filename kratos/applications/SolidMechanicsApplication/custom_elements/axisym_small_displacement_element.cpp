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
#include "custom_elements/axisym_small_displacement_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymSmallDisplacementElement::AxisymSmallDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : SmallDisplacementElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymSmallDisplacementElement::AxisymSmallDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : SmallDisplacementElement( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    //mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    //mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymSmallDisplacementElement::AxisymSmallDisplacementElement( AxisymSmallDisplacementElement const& rOther)
    :SmallDisplacementElement(rOther)
{
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer AxisymSmallDisplacementElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new AxisymSmallDisplacementElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer AxisymSmallDisplacementElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    AxisymSmallDisplacementElement NewElement ( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
      }
    
       
    return Element::Pointer( new AxisymSmallDisplacementElement(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

AxisymSmallDisplacementElement::~AxisymSmallDisplacementElement()
{
}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


//************************************************************************************
//************************************************************************************

void AxisymSmallDisplacementElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();

    rVariables.detF  = 1;

    rVariables.B.resize( 4 , number_of_nodes * 2 );

    rVariables.F.resize( 3, 3 );

    rVariables.ConstitutiveMatrix.resize( 4, 4 );

    rVariables.StrainVector.resize( 4 );

    rVariables.StressVector.resize( 4 );

    rVariables.DN_DX.resize( number_of_nodes, 2 );

    //needed parameters for consistency with the general constitutive law: small displacements
    rVariables.detF  = 1;
    rVariables.F     = identity_matrix<double>(3);

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );


}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void AxisymSmallDisplacementElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.Radius / GetProperties()[THICKNESS];

    //contributions to stiffness matrix calculated on the reference config
    SmallDisplacementElement::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )
}


//************************************************************************************
//************************************************************************************

void AxisymSmallDisplacementElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.Radius / GetProperties()[THICKNESS];

    //contribution to external forces
    SmallDisplacementElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, IntegrationWeight );

    //KRATOS_WATCH( rRightHandSideVector )
}


//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& AxisymSmallDisplacementElement::CalculateTotalMass( double& rTotalMass )
{
    KRATOS_TRY

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( GeometryData::GI_GAUSS_1 );
    const unsigned int PointNumber = 0;
    Vector N = row(Ncontainer , PointNumber);

    double Radius = 0;
    CalculateRadius (Radius, N);

    rTotalMass = GetGeometry().DomainSize() * GetProperties()[DENSITY] * 2.0 * 3.141592654 * Radius;

    return rTotalMass;

    KRATOS_CATCH( "" )
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void AxisymSmallDisplacementElement::CalculateKinematics(GeneralVariables& rVariables,
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
    CalculateRadius(rVariables.Radius, rVariables.N);

    //Displacement Gradient H [dU/dx_n]
    CalculateDisplacementGradient(rVariables.H, rVariables.DN_DX, rVariables.N, rVariables.Radius);

    //Compute the deformation matrix B
    CalculateDeformationMatrix(rVariables.B, rVariables.DN_DX, rVariables.N, rVariables.Radius);

    //Compute infinitessimal strain
    this->CalculateInfinitesimalStrain(rVariables.H, rVariables.StrainVector);


    KRATOS_CATCH( "" )
}


//*************************COMPUTE AXYSIMMETRIC RADIUS********************************
//************************************************************************************

void AxisymSmallDisplacementElement::CalculateRadius(double & rRadius,
        const Vector& rN)


{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rRadius=0;

    if ( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            array_1d<double, 3 > & ReferencePosition = GetGeometry()[i].Coordinates();

            rRadius   += ReferencePosition[0]*rN[i];
        }
    }


    if ( dimension == 3 )
    {
        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE DISPLACEMENT GRADIENT******************************
//************************************************************************************

void AxisymSmallDisplacementElement::CalculateDisplacementGradient(Matrix& rH,
        const Matrix& rDN_DX,
        const Vector & rN,
        const double & rRadius)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rH = zero_matrix<double> ( 3 );

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {

            array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            rH ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            rH ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            rH ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            rH ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
            rH ( 2 , 2 ) += rN[i]*Displacement[0]/rRadius;

        }


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


void AxisymSmallDisplacementElement::CalculateDeformationMatrix(Matrix& rB,
        const Matrix& rDN_DX,
        const Vector& rN,
        const double & rRadius)
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
            rB( 2, index + 0 ) = rN[i]/rRadius;
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

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void AxisymSmallDisplacementElement::CalculateInfinitesimalStrain(const Matrix& rH,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {

        //Infinitesimal Strain Calculation
        if ( rStrainVector.size() != 4 ) rStrainVector.resize( 4, false );

        rStrainVector[0] = rH( 0, 0 );

        rStrainVector[1] = rH( 1, 1 );

        rStrainVector[2] = rH( 2, 2 );

        rStrainVector[3] = (rH( 0, 1 ) + rH( 1, 0 )); // xy

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
void AxisymSmallDisplacementElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacementElement )
}

void AxisymSmallDisplacementElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacementElement )
}


} // Namespace Kratos


