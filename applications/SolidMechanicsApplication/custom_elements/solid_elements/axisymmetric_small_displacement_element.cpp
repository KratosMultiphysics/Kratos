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
#include "custom_elements/solid_elements/axisymmetric_small_displacement_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymmetricSmallDisplacementElement::AxisymmetricSmallDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : SmallDisplacementElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymmetricSmallDisplacementElement::AxisymmetricSmallDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : SmallDisplacementElement( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    //mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    //mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymmetricSmallDisplacementElement::AxisymmetricSmallDisplacementElement( AxisymmetricSmallDisplacementElement const& rOther)
    :SmallDisplacementElement(rOther)
{
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer AxisymmetricSmallDisplacementElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new AxisymmetricSmallDisplacementElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer AxisymmetricSmallDisplacementElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    AxisymmetricSmallDisplacementElement NewElement ( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
      }
    
    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer( new AxisymmetricSmallDisplacementElement(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

AxisymmetricSmallDisplacementElement::~AxisymmetricSmallDisplacementElement()
{
}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


//************************************************************************************
//************************************************************************************

void AxisymmetricSmallDisplacementElement::InitializeElementData (ElementDataPointerType & pVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType& dimension       = this->Dimension();
    const unsigned int voigt_size      = 4;
    
    pVariables->Initialize(voigt_size,dimension,number_of_nodes);

    //needed parameters for consistency with the general constitutive law: small displacements
    pVariables->detF  = 1;
    pVariables->F     = identity_matrix<double>(3);

    //set variables including all integration points values

    //reading shape functions
    pVariables->SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    pVariables->SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    pVariables->j = GetGeometry().Jacobian( pVariables->j, mThisIntegrationMethod );

    //Calculate Delta Position
    pVariables->DeltaPosition = this->CalculateTotalDeltaPosition(pVariables->DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    pVariables->J = GetGeometry().Jacobian( pVariables->J, mThisIntegrationMethod, pVariables->DeltaPosition );


}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void AxisymmetricSmallDisplacementElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataPointerType& pVariables, double& rIntegrationWeight)
{
  
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * pVariables->ReferenceRadius;
    if ( this->GetProperties().Has( THICKNESS ) )
      IntegrationWeight /= GetProperties()[THICKNESS];
  
    //contributions to stiffness matrix calculated on the reference config
    SmallDisplacementElement::CalculateAndAddLHS( rLocalSystem, pVariables, IntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )
}


//************************************************************************************
//************************************************************************************

void AxisymmetricSmallDisplacementElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataPointerType& pVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * pVariables->ReferenceRadius;
    if ( this->GetProperties().Has( THICKNESS ) )
      IntegrationWeight /= GetProperties()[THICKNESS];
  
    //contribution to external forces
    SmallDisplacementElement::CalculateAndAddRHS( rLocalSystem, pVariables, rVolumeForce, IntegrationWeight );

    //KRATOS_WATCH( rRightHandSideVector )
}


//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& AxisymmetricSmallDisplacementElement::CalculateTotalMass( double& rTotalMass, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Compute the Volume Change acumulated:
    ElementDataPointerType Variables(make_unique<ElementDataType>());
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	//compute element kinematics
	this->CalculateKinematics(Variables,PointNumber);
	
	//getting informations for integration
        double IntegrationWeight = Variables->detJ * integration_points[PointNumber].Weight() * 2.0 * 3.141592654 * Variables->ReferenceRadius;

	//compute point volume changes	
	rTotalMass += GetProperties()[DENSITY] * IntegrationWeight;
      }


    return rTotalMass;

    KRATOS_CATCH( "" )
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void AxisymmetricSmallDisplacementElement::CalculateKinematics(ElementDataPointerType& pVariables,
        const double& rPointNumber)

{
    KRATOS_TRY
      
    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = pVariables->GetShapeFunctionsGradients();
    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = pVariables->GetShapeFunctions();

    //Parent to reference configuration
    pVariables->StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( pVariables->J[rPointNumber], InvJ, pVariables->detJ);

    //Compute cartesian derivatives [dN/dx_n]
    noalias( pVariables->DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Set Shape Functions Values for this integration point
    noalias(pVariables->N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Calculate IntegrationPoint radius
    CalculateRadius(pVariables->ReferenceRadius, pVariables->N);

    //Displacement Gradient H [dU/dx_n]
    CalculateDisplacementGradient(pVariables->H, pVariables->DN_DX, pVariables->N, pVariables->ReferenceRadius);

    //Compute the deformation matrix B
    CalculateDeformationMatrix(pVariables->B, pVariables->DN_DX, pVariables->N, pVariables->ReferenceRadius);

    //Compute infinitessimal strain
    this->CalculateInfinitesimalStrain(pVariables->H, pVariables->StrainVector);


    KRATOS_CATCH( "" )
}


//*************************COMPUTE AXYSIMMETRIC RADIUS********************************
//************************************************************************************

void AxisymmetricSmallDisplacementElement::CalculateRadius(double & rRadius,
        const Vector& rN)


{

    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();

    const SizeType& dimension = this->Dimension();

    rRadius=0;

    if ( dimension == 2 )
    {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
            //array_1d<double, 3 > & ReferencePosition = GetGeometry()[i].Coordinates();
            //rRadius   += ReferencePosition[0]*rN[i];

           rRadius += rN[i] * GetGeometry()[i].X0(); 
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

void AxisymmetricSmallDisplacementElement::CalculateDisplacementGradient(Matrix& rH,
        const Matrix& rDN_DX,
        const Vector & rN,
        const double & rRadius)

{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();

    const SizeType& dimension = this->Dimension();

    rH = zero_matrix<double> ( 3 );

    if( dimension == 2 )
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
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


void AxisymmetricSmallDisplacementElement::CalculateDeformationMatrix(Matrix& rB,
        const Matrix& rDN_DX,
        const Vector& rN,
        const double & rRadius)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType& dimension = this->Dimension();

    rB.clear(); //set all components to zero

    if( dimension == 2 )
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
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

void AxisymmetricSmallDisplacementElement::CalculateInfinitesimalStrain(const Matrix& rH,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const SizeType& dimension = this->Dimension();

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
  
int AxisymmetricSmallDisplacementElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = SmallDisplacementElement::Check(rCurrentProcessInfo);     

    return ErrorCode;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void AxisymmetricSmallDisplacementElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacementElement )
}

void AxisymmetricSmallDisplacementElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacementElement )
}


} // Namespace Kratos


