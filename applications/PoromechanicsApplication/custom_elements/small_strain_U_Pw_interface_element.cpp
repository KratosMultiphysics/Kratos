//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

/* Project includes */
#include "includes/define.h"
#include "custom_elements/small_strain_U_Pw_interface_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
SmallStrainUPwInterfaceElement::SmallStrainUPwInterfaceElement() : Element() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SmallStrainUPwInterfaceElement::SmallStrainUPwInterfaceElement( IndexType NewId, GeometryType::Pointer pGeometry ) : Element( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SmallStrainUPwInterfaceElement::SmallStrainUPwInterfaceElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : Element( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
SmallStrainUPwInterfaceElement::~SmallStrainUPwInterfaceElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Element::Pointer SmallStrainUPwInterfaceElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallStrainUPwInterfaceElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

int  SmallStrainUPwInterfaceElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    //verify that the variables are correctly initialized

    //Solid variables
    if ( DISPLACEMENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )

    if ( VELOCITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" )

    if ( ACCELERATION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" )

    if ( DENSITY_SOLID.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY_SOLID has Key zero! (check if the application is correctly registered", "" )

    //Fluid variables
    if ( WATER_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "WATER_PRESSURE has Key zero! (check if the application is correctly registered", "" )

    if ( DERIVATIVE_WATER_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DERIVATIVE_WATER_PRESSURE has Key zero! (check if the application is correctly registered", "" )

    if ( DENSITY_WATER.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY_WATER has Key zero! (check if the application is correctly registered", "" )

    //Interface variables
    if ( FRACTURE_APERTURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "FRACTURE_APERTURE has Key zero! (check if the application is correctly registered", "" )

    //verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )


        if ( this->GetGeometry()[i].SolutionStepsDataHas( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable WATER_PRESSURE on node ", this->GetGeometry()[i].Id() )

        if ( this->GetGeometry()[i].HasDofFor( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable WATER_PRESSURE on node ", GetGeometry()[i].Id() )
    }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW_POINTER ) == false )
        KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() )

	//verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW_POINTER )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
	    if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
		    correct_strain_measure = true;
    }

    if( correct_strain_measure == false )
	    KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the element type ", " StrainMeasure_Infinitesimal " );

    //verify that the constitutive law has the correct dimension
    if ( dimension == 2 )
    {
        if ( this->GetProperties().Has( THICKNESS ) == false )
            KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() )

        if ( THICKNESS.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered)", "" )
    }

	this->GetProperties().GetValue( CONSTITUTIVE_LAW_POINTER )->Check( this->GetProperties(), this->GetGeometry(), rCurrentProcessInfo );

    return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::Initialize()
{
    KRATOS_TRY

	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    if ( GetProperties()[CONSTITUTIVE_LAW_POINTER] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW_POINTER]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element with ID ", this->Id() )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    rElementalDofList.resize( 0 );
    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

        rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ));
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dimension + 1);

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != MatSize )
        rLeftHandSideMatrix.resize( MatSize, MatSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != MatSize )
        rRightHandSideVector.resize( MatSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( MatSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int MatSize = number_of_nodes * (dimension + 1);

    //Resetting mass matrix
    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );
    noalias( rMassMatrix ) = ZeroMatrix( MatSize, MatSize );

    //Defining shape functions and the determinant of the jacobian at all integration points
    Matrix Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    //TODO: s'ha de calcular la matriu de rotacio per obtenir els desplaçaments normals de la junta (utilitat que rebi la posicio dels nodes?)
    Vector detJcontainer;
    GetGeometry().DeterminantOfJacobian(detJcontainer,mThisIntegrationMethod);

    //Loop over integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    double IntegrationCoefficient;
    double Porosity = GetProperties()[POROSITY];
    double Density = Porosity*GetProperties()[DENSITY_WATER] + (1-Porosity)*GetProperties()[DENSITY_SOLID];
    unsigned int node;
    Matrix Nut = ZeroMatrix( dimension + 1 , number_of_nodes * (dimension + 1) );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //Setting the shape function matrix
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = (dimension+1)*i;
            Nut(0,node)=Ncontainer(PointNumber,i);
            Nut(1,node+1)=Ncontainer(PointNumber,i);
            if(dimension==3)
                Nut(2,node+2)=Ncontainer(PointNumber,i);
        }

        //calculating weighting coefficient for integration
        this->CalculateIntegrationCoefficient( IntegrationCoefficient, detJcontainer[PointNumber], integration_points[PointNumber].Weight() );

        //Adding contribution to Mass matrix
        noalias(rMassMatrix) += Density*prod(trans(Nut),Nut)*IntegrationCoefficient;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

} // Namespace Kratos
