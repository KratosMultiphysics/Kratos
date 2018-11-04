//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"

// Application includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"

namespace Kratos
{

// Default Constructor
SmallStrainUPwDiffOrderElement::SmallStrainUPwDiffOrderElement() : Element() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SmallStrainUPwDiffOrderElement::SmallStrainUPwDiffOrderElement( IndexType NewId, GeometryType::Pointer pGeometry ) : Element( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SmallStrainUPwDiffOrderElement::SmallStrainUPwDiffOrderElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : Element( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
SmallStrainUPwDiffOrderElement::~SmallStrainUPwDiffOrderElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Element::Pointer SmallStrainUPwDiffOrderElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallStrainUPwDiffOrderElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

int  SmallStrainUPwDiffOrderElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    unsigned int dimension = rGeom.WorkingSpaceDimension();

    //verify that the variables are correctly initialized

    //Solid variables
    if ( DISPLACEMENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT Key is 0. Check if all applications were correctly registered.", "" )

    if ( VELOCITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY Key is 0. Check if all applications were correctly registered.", "" )

    if ( ACCELERATION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION Key is 0. Check if all applications were correctly registered.", "" )

    if ( DENSITY_SOLID.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY_SOLID Key is 0. Check if all applications were correctly registered.", "" )

    //Fluid variables
    if ( WATER_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "WATER_PRESSURE Key is 0. Check if all applications were correctly registered.", "" )

    if ( DT_WATER_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DT_WATER_PRESSURE Key is 0. Check if all applications were correctly registered.", "" )

    if ( DENSITY_WATER.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY_WATER Key is 0. Check if all applications were correctly registered.", "" )

    //verify that the dofs exist
    for ( unsigned int i = 0; i < rGeom.size(); i++ )
    {
        if ( rGeom[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", rGeom[i].Id() )

        if ( rGeom[i].HasDofFor( DISPLACEMENT_X ) == false || rGeom[i].HasDofFor( DISPLACEMENT_Y ) == false || rGeom[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", rGeom[i].Id() )


        if ( rGeom[i].SolutionStepsDataHas( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable WATER_PRESSURE on node ", rGeom[i].Id() )

        if ( rGeom[i].HasDofFor( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable WATER_PRESSURE on node ", rGeom[i].Id() )
    }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() )

	//verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

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
            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS Key is 0. Check if all applications were correctly registered.)", "" )
    }

	this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), rGeom, rCurrentProcessInfo );

    return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::Initialize()
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
	const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), rGeom,row( rGeom.ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element with ID ", this->Id() )


    const SizeType NumUNodes = rGeom.PointsNumber();

    switch(NumUNodes)
    {
        case 6: //2D T6P3
            mpPressureGeometry = GeometryType::Pointer( new Triangle2D3< Node<3> >(rGeom(0), rGeom(1), rGeom(2)) );
            break;
        case 8: //2D Q8P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral2D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 9: //2D Q9P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral2D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 10: //3D T10P4
            mpPressureGeometry = GeometryType::Pointer( new Tetrahedra3D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 20: //3D H20P8
            mpPressureGeometry = GeometryType::Pointer( new Hexahedra3D8< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7)) );
            break;
        case 27: //3D H27P8
            mpPressureGeometry = GeometryType::Pointer( new Hexahedra3D8< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7)) );
            break;
        default:
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected geometry type for different order interpolation element","");
            break;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    if(rElementalDofList.size() != ElementSize)
        rElementalDofList.resize(ElementSize);

    SizeType Index = 0;
/*
    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        if(Dim > 2)
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( WATER_PRESSURE );
    }

    for(SizeType i=NumPNodes; i<NumUNodes; i++)
    {
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        if(Dim > 2)
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
    }
*/

    for(SizeType i = 0; i < NumUNodes; i++)
    {
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        if(Dim > 2)
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
    }

    for(SizeType i=0; i<NumPNodes; i++)
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( WATER_PRESSURE );

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != ElementSize )
        rLeftHandSideMatrix.resize( ElementSize, ElementSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( ElementSize, ElementSize );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != ElementSize )
        rRightHandSideVector.resize( ElementSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ElementSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    KRATOS_THROW_ERROR(std::logic_error,"SmallStrainUPwDiffOrderElement::CalculateLeftHandSide not implemented","");
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    //Resetting the RHS
    if ( rRightHandSideVector.size() != ElementSize )
        rRightHandSideVector.resize( ElementSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ElementSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType BlockElementSize = NumUNodes * Dim;
    const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints( mThisIntegrationMethod );
    const SizeType NumGPoints = integration_points.size();
    
    Matrix M = ZeroMatrix(BlockElementSize, BlockElementSize);

    //Defining shape functions and the determinant of the jacobian at all integration points
    Matrix Nucontainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
    Vector detJcontainer = ZeroVector(NumGPoints);
    rGeom.DeterminantOfJacobian(detJcontainer,mThisIntegrationMethod);

    //Loop over integration points
    double IntegrationCoefficient;
    double Porosity = GetProperties()[POROSITY];
    double Density = Porosity*GetProperties()[DENSITY_WATER] + (1.0-Porosity)*GetProperties()[DENSITY_SOLID];
    SizeType Index = 0;
    Matrix Nu = ZeroMatrix( Dim , NumUNodes * Dim );

    for ( SizeType PointNumber = 0; PointNumber < NumGPoints; PointNumber++ )
    {
        //Setting the shape function matrix
        for ( SizeType i = 0; i < NumUNodes; i++ )
        {
            Nu(0,Index++)=Nucontainer(PointNumber,i);
            Nu(1,Index++)=Nucontainer(PointNumber,i);
            if(Dim>2)
                Nu(2,Index++)=Nucontainer(PointNumber,i);
        }

        //calculating weighting coefficient for integration
        this->CalculateIntegrationCoefficient( IntegrationCoefficient, detJcontainer[PointNumber], integration_points[PointNumber].Weight() );

        //Adding contribution to Mass matrix
        noalias(M) += Density*prod(trans(Nu),Nu)*IntegrationCoefficient;
    }

    //Distribute mass block matrix into the elemental matrix
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    if ( rMassMatrix.size1() != ElementSize )
        rMassMatrix.resize( ElementSize, ElementSize, false );
    noalias( rMassMatrix ) = ZeroMatrix( ElementSize, ElementSize );

    SizeType Index_i, Index_j;

    for(SizeType i = 0; i < NumUNodes; i++)
    {
        Index_i = i * Dim;

        for(SizeType j = 0; j < NumUNodes; j++)
        {
            Index_j = j * Dim;

            rMassMatrix(Index_i,Index_j)     += M(Index_i,Index_j);
            rMassMatrix(Index_i,Index_j+1)   += M(Index_i,Index_j+1);
            rMassMatrix(Index_i+1,Index_j)   += M(Index_i+1,Index_j);
            rMassMatrix(Index_i+1,Index_j+1) += M(Index_i+1,Index_j+1);
            if(Dim > 2)
            {
                rMassMatrix(Index_i,Index_j+2)   += M(Index_i,Index_j+2);
                rMassMatrix(Index_i+1,Index_j+2) += M(Index_i+1,Index_j+2);
                rMassMatrix(Index_i+2,Index_j)   += M(Index_i+2,Index_j);
                rMassMatrix(Index_i+2,Index_j+1) += M(Index_i+2,Index_j+1);
                rMassMatrix(Index_i+2,Index_j+2) += M(Index_i+2,Index_j+2);
            }
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    if ( rResult.size() != ElementSize )
        rResult.resize( ElementSize, false );

    SizeType Index = 0;

    for ( SizeType i = 0; i < NumUNodes; i++ )
    {
        rResult[Index++] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[Index++] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        if(Dim > 2)
            rResult[Index++] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    for ( SizeType i = 0; i < NumPNodes; i++ )
        rResult[Index++] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    if ( rValues.size() != ElementSize )
        rValues.resize( ElementSize, false );

    SizeType Index = 0;

    for ( SizeType i = 0; i < NumUNodes; i++ )
    {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
        if ( Dim > 2 )
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
    }

    for ( SizeType i = 0; i < NumPNodes; i++ )
        rValues[Index++] = 0.0;
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    //Definition of variables
    ElementalVariables Variables;
    this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        //compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables,PointNumber);

        //set gauss points variables to constitutivelaw parameters
        this->SetElementalVariables(Variables,ConstitutiveParameters);

        //compute constitutive tensor and/or stresses
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
    }


    //Assign pressure values to the intermediate nodes for post-processing
    GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();

    switch (NumUNodes)
    {
        case 6: //2D T6P3
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[3],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[4],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[5],WATER_PRESSURE, 0.5 * (p2 + p0) );
            break;
        }
        case 8: //2D Q8P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],WATER_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[7],WATER_PRESSURE, 0.5 * (p3 + p0) );
            break;
        }
        case 9: //2D Q9P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],WATER_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[7],WATER_PRESSURE, 0.5 * (p3 + p0) );
            ThreadSafeNodeWrite(rGeom[8],WATER_PRESSURE, 0.25 * (p0 + p1 + p2 + p3) );
            break;
        }
        case 10: //3D T10P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],WATER_PRESSURE, 0.5 * (p2 + p0) );
            ThreadSafeNodeWrite(rGeom[7],WATER_PRESSURE, 0.5 * (p0 + p3) );
            ThreadSafeNodeWrite(rGeom[8],WATER_PRESSURE, 0.5 * (p1 + p3) );
            ThreadSafeNodeWrite(rGeom[9],WATER_PRESSURE, 0.5 * (p2 + p3) );
            break;
        }
        case 20: //3D H20P8
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p4 = rGeom[4].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p5 = rGeom[5].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p6 = rGeom[6].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p7 = rGeom[7].FastGetSolutionStepValue(WATER_PRESSURE);
            // edges -- bottom
            ThreadSafeNodeWrite(rGeom[8],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[9],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[10],WATER_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[11],WATER_PRESSURE, 0.5 * (p3 + p0) );
            // edges -- middle
            ThreadSafeNodeWrite(rGeom[12],WATER_PRESSURE, 0.5 * (p4 + p0) );
            ThreadSafeNodeWrite(rGeom[13],WATER_PRESSURE, 0.5 * (p5 + p1) );
            ThreadSafeNodeWrite(rGeom[14],WATER_PRESSURE, 0.5 * (p6 + p2) );
            ThreadSafeNodeWrite(rGeom[15],WATER_PRESSURE, 0.5 * (p7 + p3) );
            // edges -- top
            ThreadSafeNodeWrite(rGeom[16],WATER_PRESSURE, 0.5 * (p4 + p5) );
            ThreadSafeNodeWrite(rGeom[17],WATER_PRESSURE, 0.5 * (p5 + p6) );
            ThreadSafeNodeWrite(rGeom[18],WATER_PRESSURE, 0.5 * (p6 + p7) );
            ThreadSafeNodeWrite(rGeom[19],WATER_PRESSURE, 0.5 * (p7 + p0) );
            break;
        }
        case 27: //3D H27P8
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p4 = rGeom[4].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p5 = rGeom[5].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p6 = rGeom[6].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p7 = rGeom[7].FastGetSolutionStepValue(WATER_PRESSURE);
            // edges -- bottom
            ThreadSafeNodeWrite(rGeom[8],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[9],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[10],WATER_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[11],WATER_PRESSURE, 0.5 * (p3 + p0) );
            // edges -- middle
            ThreadSafeNodeWrite(rGeom[12],WATER_PRESSURE, 0.5 * (p4 + p0) );
            ThreadSafeNodeWrite(rGeom[13],WATER_PRESSURE, 0.5 * (p5 + p1) );
            ThreadSafeNodeWrite(rGeom[14],WATER_PRESSURE, 0.5 * (p6 + p2) );
            ThreadSafeNodeWrite(rGeom[15],WATER_PRESSURE, 0.5 * (p7 + p3) );
            // edges -- top
            ThreadSafeNodeWrite(rGeom[16],WATER_PRESSURE, 0.5 * (p4 + p5) );
            ThreadSafeNodeWrite(rGeom[17],WATER_PRESSURE, 0.5 * (p5 + p6) );
            ThreadSafeNodeWrite(rGeom[18],WATER_PRESSURE, 0.5 * (p6 + p7) );
            ThreadSafeNodeWrite(rGeom[19],WATER_PRESSURE, 0.5 * (p7 + p0) );
            // face centers
            ThreadSafeNodeWrite(rGeom[20],WATER_PRESSURE, 0.25 * (p0 + p1 + p2 + p3) );
            ThreadSafeNodeWrite(rGeom[21],WATER_PRESSURE, 0.25 * (p0 + p1 + p4 + p5) );
            ThreadSafeNodeWrite(rGeom[22],WATER_PRESSURE, 0.25 * (p1 + p2 + p5 + p6) );
            ThreadSafeNodeWrite(rGeom[23],WATER_PRESSURE, 0.25 * (p2 + p3 + p6 + p7) );
            ThreadSafeNodeWrite(rGeom[24],WATER_PRESSURE, 0.25 * (p3 + p0 + p7 + p4) );
            ThreadSafeNodeWrite(rGeom[25],WATER_PRESSURE, 0.25 * (p4 + p5 + p6 + p7) );
            // element center
            ThreadSafeNodeWrite(rGeom[26],WATER_PRESSURE, 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7) );
            break;
        }
        default:
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected geometry type for different order interpolation element","");
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable,std::vector<Vector>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,std::vector<Matrix>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == VON_MISES_STRESS )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    else
    {
        const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

        if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

        for ( unsigned int i = 0; i < integration_points_number; i++ )
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,std::vector<Vector>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    else if ( rVariable == FLUID_FLUX_VECTOR )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    else
    {
        const unsigned int& integration_points_number = mConstitutiveLawVector.size();

        if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

        for ( unsigned int i = 0;  i < integration_points_number; i++ )
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,std::vector<Matrix>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    else
    {
        const unsigned int& integration_points_number = mConstitutiveLawVector.size();

        if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

        for ( unsigned int i = 0;  i < integration_points_number; i++ )
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,std::vector<ConstitutiveLaw::Pointer>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if(rVariable == CONSTITUTIVE_LAW)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for(unsigned int i=0; i<rValues.size(); i++)
            rValues[i] = mConstitutiveLawVector[i];
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int& integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number, false );

    if ( rVariable == VON_MISES_STRESS )
    {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

        //Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            //set gauss points variables to constitutivelaw parameters
            this->SetElementalVariables(Variables,ConstitutiveParameters);

            //compute constitutive tensor and/or stresses
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            rOutput[PointNumber] =  ElementUtilities::CalculateVonMises(Variables.StressVector);
        }
    }
    else
    {
        for ( unsigned int i = 0; i < integration_points_number; i++ )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int& integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_VECTOR )
    {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

        //Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            //set gauss points variables to constitutivelaw parameters
            this->SetElementalVariables(Variables,ConstitutiveParameters);

            //compute constitutive tensor and/or stresses
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            if ( rOutput[PointNumber].size() != Variables.StressVector.size() )
                rOutput[PointNumber].resize( Variables.StressVector.size(), false );

            rOutput[PointNumber] = Variables.StressVector;
        }
    }
    else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
    {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            if ( rOutput[PointNumber].size() != Variables.StrainVector.size() )
                rOutput[PointNumber].resize( Variables.StrainVector.size(), false );

            rOutput[PointNumber] = Variables.StrainVector;
        }
    }
    else if ( rVariable == FLUID_FLUX_VECTOR )
    {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            //Compute FluidFlux vector q [m/s]
            const SizeType Dim = rGeom.WorkingSpaceDimension();
            const SizeType NumUNodes = rGeom.PointsNumber();
            Vector BodyAcceleration = ZeroVector(Dim);
            SizeType Index = 0;
            for(SizeType i = 0; i < NumUNodes; i++)
            {
                BodyAcceleration[0] += Variables.Nu[i]*Variables.BodyAcceleration[Index++];
                BodyAcceleration[1] += Variables.Nu[i]*Variables.BodyAcceleration[Index++];
                if (Dim>2)
                    BodyAcceleration[2] += Variables.Nu[i]*Variables.BodyAcceleration[Index++];
            }
            
            Vector GradPressureTerm(Dim);
            noalias(GradPressureTerm) = prod(trans(Variables.GradNpT),Variables.PressureVector);
            noalias(GradPressureTerm) -= GetProperties()[DENSITY_WATER]*BodyAcceleration;

            Vector AuxFluidFlux = ZeroVector(Dim);
            AuxFluidFlux = - 1.0/Variables.DynamicViscosity * prod(Variables.IntrinsicPermeability, GradPressureTerm );

            Vector FluidFlux = ZeroVector(3);
            FluidFlux[0] = AuxFluidFlux[0];
            FluidFlux[1] = AuxFluidFlux[1];
            if(Dim>2)
                FluidFlux[2] = AuxFluidFlux[2];

            if ( rOutput[PointNumber].size() != 3 )
                rOutput[PointNumber].resize( 3, false );

            rOutput[PointNumber] = FluidFlux;
        }
    }
    else
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable , rOutput[i] );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int& integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    const unsigned int dimension       = rGeom.WorkingSpaceDimension();

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_TENSOR )
    {
        std::vector<Vector> StressVector;

        this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);
        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
    {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            rOutput[PointNumber] = MathUtils<double>::StrainVectorToTensor(StrainVector[PointNumber]);
        }
    }
    else
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable , rOutput[i] );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
                                bool CalculateLHSMatrixFlag, bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    
    //Definition of variables
    ElementalVariables Variables;
    this->InitializeElementalVariables(Variables,rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,GetProperties(),rCurrentProcessInfo);
    if(CalculateLHSMatrixFlag)
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if(CalculateResidualVectorFlag)
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

    //Loop over integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints( mThisIntegrationMethod );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables,PointNumber);

        //set gauss points variables to constitutivelaw parameters
        this->SetElementalVariables(Variables,ConstitutiveParameters);

        //compute constitutive tensor and/or stresses
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //calculating weighting coefficient for integration
        this->CalculateIntegrationCoefficient( Variables.IntegrationCoefficient, Variables.detJuContainer[PointNumber], integration_points[PointNumber].Weight() );

        //Contributions to the left hand side
        if ( CalculateLHSMatrixFlag )
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if ( CalculateResidualVectorFlag )
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::InitializeElementalVariables (ElementalVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    
    //Variables at all integration points
    (rVariables.NuContainer).resize(NumGPoints,NumUNodes,false);
    rVariables.NuContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
    
    (rVariables.NpContainer).resize(NumGPoints,NumPNodes,false);
    rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

    (rVariables.Nu).resize(NumUNodes,false);
    (rVariables.Np).resize(NumPNodes,false);
    
    (rVariables.DNu_DXContainer).resize(NumGPoints,false);
    for(SizeType i = 0; i<NumGPoints; i++)
        ((rVariables.DNu_DXContainer)[i]).resize(NumUNodes,Dim,false);
    (rVariables.DNu_DX).resize(NumUNodes,Dim,false);
    (rVariables.detJuContainer).resize(NumGPoints,false);
    rGeom.ShapeFunctionsIntegrationPointsGradients(rVariables.DNu_DXContainer,rVariables.detJuContainer,mThisIntegrationMethod);

    (rVariables.DNp_DXContainer).resize(NumGPoints,false);
    for(SizeType i = 0; i<NumGPoints; i++)
        ((rVariables.DNp_DXContainer)[i]).resize(NumPNodes,Dim,false);
    (rVariables.GradNpT).resize(NumPNodes,Dim,false);
    Vector detJpContainer = ZeroVector(NumGPoints);
    mpPressureGeometry->ShapeFunctionsIntegrationPointsGradients(rVariables.DNp_DXContainer,detJpContainer,mThisIntegrationMethod);
    
    //Variables computed at each integration point    
    unsigned int voigtsize  = 3;
    if( Dim == 3 ) voigtsize  = 6;
    (rVariables.B).resize(voigtsize, NumUNodes * Dim, false);
    noalias(rVariables.B) = ZeroMatrix( voigtsize, NumUNodes * Dim );
    (rVariables.StrainVector).resize(voigtsize,false);
    (rVariables.ConstitutiveMatrix).resize(voigtsize, voigtsize, false);
    (rVariables.StressVector).resize(voigtsize,false);
    
    //Needed parameters for consistency with the general constitutive law
    rVariables.detF  = 1.0;
    (rVariables.F).resize(Dim, Dim, false);
    noalias(rVariables.F) = identity_matrix<double>(Dim);

    //Nodal variables
    this->InitializeNodalVariables(rVariables);

    //Properties variables
    this->InitializeProperties(rVariables);

    //ProcessInfo variables
    rVariables.NewmarkCoefficient1 = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.NewmarkCoefficient2 = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::InitializeNodalVariables (ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    SizeType Local_i;
    Vector BodyAccelerationAux    = ZeroVector(3);
    (rVariables.BodyAcceleration).resize(NumUNodes * Dim,false);
    (rVariables.DisplacementVector).resize(NumUNodes * Dim,false);
    (rVariables.VelocityVector).resize(NumUNodes * Dim,false);
    for(SizeType i=0; i<NumUNodes; i++)
    {
        Local_i = i * Dim;
        BodyAccelerationAux = rGeom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        rVariables.BodyAcceleration[Local_i]   = BodyAccelerationAux[0];
        rVariables.DisplacementVector[Local_i] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
        rVariables.VelocityVector[Local_i]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);

        rVariables.BodyAcceleration[Local_i+1]   = BodyAccelerationAux[1];
        rVariables.DisplacementVector[Local_i+1] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        rVariables.VelocityVector[Local_i+1]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);

        if(Dim >2)
        {
            rVariables.BodyAcceleration[Local_i+2]   = BodyAccelerationAux[2];
            rVariables.DisplacementVector[Local_i+2] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
            rVariables.VelocityVector[Local_i+2]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z);
        }
    }

    (rVariables.PressureVector).resize(NumPNodes,false);
    (rVariables.PressureDtVector).resize(NumPNodes,false);
    for(SizeType i=0; i<NumPNodes; i++)
    {
        rVariables.PressureVector[i]   = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.PressureDtVector[i] = rGeom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::InitializeProperties (ElementalVariables& rVariables)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double BulkModulus = GetProperties()[YOUNG_MODULUS]/(3.0*(1.0-2.0*GetProperties()[POISSON_RATIO]));
    double BulkModulusSolid = GetProperties()[BULK_MODULUS_SOLID];
    rVariables.BiotCoefficient = 1.0-BulkModulus/BulkModulusSolid;
    double Porosity = GetProperties()[POROSITY];
    rVariables.BiotModulusInverse = (rVariables.BiotCoefficient-Porosity)/BulkModulusSolid + Porosity/GetProperties()[BULK_MODULUS_FLUID];
    rVariables.DynamicViscosity = GetProperties()[DYNAMIC_VISCOSITY];
    //Setting the intrinsic permeability matrix
    (rVariables.IntrinsicPermeability).resize(dimension,dimension,false);
    rVariables.IntrinsicPermeability(0,0) = GetProperties()[PERMEABILITY_XX];
    rVariables.IntrinsicPermeability(1,1) = GetProperties()[PERMEABILITY_YY];
    rVariables.IntrinsicPermeability(0,1) = GetProperties()[PERMEABILITY_XY];
    rVariables.IntrinsicPermeability(1,0) = rVariables.IntrinsicPermeability(0,1);
    if(dimension==3)
    {
        rVariables.IntrinsicPermeability(2,2) = GetProperties()[PERMEABILITY_ZZ];
        rVariables.IntrinsicPermeability(2,0) = GetProperties()[PERMEABILITY_ZX];
        rVariables.IntrinsicPermeability(1,2) = GetProperties()[PERMEABILITY_YZ];
        rVariables.IntrinsicPermeability(0,2) = rVariables.IntrinsicPermeability(2,0);
        rVariables.IntrinsicPermeability(2,1) = rVariables.IntrinsicPermeability(1,2);
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateKinematics(ElementalVariables& rVariables, unsigned int PointNumber)

{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();

    //Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Nu) = row(rVariables.NuContainer, PointNumber);
    noalias(rVariables.Np) = row(rVariables.NpContainer, PointNumber);

    noalias(rVariables.DNu_DX) = rVariables.DNu_DXContainer[PointNumber];
    noalias(rVariables.GradNpT) = rVariables.DNp_DXContainer[PointNumber];

    SizeType node;

    //Compute the deformation matrix B
    if( Dim == 2 )
    {
        for ( SizeType i = 0; i < NumUNodes; i++ )
        {
            node = 2 * i;

            rVariables.B( 0, node + 0 ) = rVariables.DNu_DX( i, 0 );
            rVariables.B( 1, node + 1 ) = rVariables.DNu_DX( i, 1 );
            rVariables.B( 2, node + 0 ) = rVariables.DNu_DX( i, 1 );
            rVariables.B( 2, node + 1 ) = rVariables.DNu_DX( i, 0 );
        }
    }
    else
    {
        for ( SizeType i = 0; i < NumUNodes; i++ )
        {
            node = 3 * i;

            rVariables.B( 0, node + 0 ) = rVariables.DNu_DX( i, 0 );
            rVariables.B( 1, node + 1 ) = rVariables.DNu_DX( i, 1 );
            rVariables.B( 2, node + 2 ) = rVariables.DNu_DX( i, 2 );
            rVariables.B( 3, node + 0 ) = rVariables.DNu_DX( i, 1 );
            rVariables.B( 3, node + 1 ) = rVariables.DNu_DX( i, 0 );
            rVariables.B( 4, node + 1 ) = rVariables.DNu_DX( i, 2 );
            rVariables.B( 4, node + 2 ) = rVariables.DNu_DX( i, 1 );
            rVariables.B( 5, node + 0 ) = rVariables.DNu_DX( i, 2 );
            rVariables.B( 5, node + 2 ) = rVariables.DNu_DX( i, 0 );
        }
    }

    //Compute infinitessimal strain
    rVariables.StrainVector = prod(rVariables.B,rVariables.DisplacementVector);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::SetElementalVariables(ElementalVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters)
{
    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rConstitutiveParameters.SetStressVector(rVariables.StressVector);

    //Needed parameters for consistency with the general constitutive law
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.DNu_DX);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Nu);

    rConstitutiveParameters.SetDeterminantF(rVariables.detF);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, double detJu, double weight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rIntegrationCoefficient = detJu * weight;

    if( dimension == 2 )
        rIntegrationCoefficient *= GetProperties()[THICKNESS];
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
}


//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    Matrix StiffnessMatrix = prod(trans(rVariables.B), Matrix(prod(rVariables.ConstitutiveMatrix, rVariables.B)))*rVariables.IntegrationCoefficient;

    //Distribute stiffness block matrix into the elemental matrix
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    SizeType Index_i, Index_j;

    for(SizeType i = 0; i < NumUNodes; i++)
    {
        Index_i = i * Dim;

        for(SizeType j = 0; j < NumUNodes; j++)
        {
            Index_j = j * Dim;

            rLeftHandSideMatrix(Index_i,Index_j)     += StiffnessMatrix(Index_i,Index_j);
            rLeftHandSideMatrix(Index_i,Index_j+1)   += StiffnessMatrix(Index_i,Index_j+1);
            rLeftHandSideMatrix(Index_i+1,Index_j)   += StiffnessMatrix(Index_i+1,Index_j);
            rLeftHandSideMatrix(Index_i+1,Index_j+1) += StiffnessMatrix(Index_i+1,Index_j+1);
            if(Dim > 2)
            {
                rLeftHandSideMatrix(Index_i,Index_j+2)   += StiffnessMatrix(Index_i,Index_j+2);
                rLeftHandSideMatrix(Index_i+1,Index_j+2) += StiffnessMatrix(Index_i+1,Index_j+2);
                rLeftHandSideMatrix(Index_i+2,Index_j)   += StiffnessMatrix(Index_i+2,Index_j);
                rLeftHandSideMatrix(Index_i+2,Index_j+1) += StiffnessMatrix(Index_i+2,Index_j+1);
                rLeftHandSideMatrix(Index_i+2,Index_j+2) += StiffnessMatrix(Index_i+2,Index_j+2);
            }
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    unsigned int voigtsize  = 3;
    if( Dim == 3 ) voigtsize  = 6;

    Vector VoigtVector = ZeroVector(voigtsize);
    VoigtVector[0] = 1.0;
    VoigtVector[1] = 1.0;
    if(Dim == 3) VoigtVector[2] = 1.0;

    Matrix CouplingMatrix = rVariables.BiotCoefficient*prod(trans(rVariables.B),Matrix(outer_prod(VoigtVector,rVariables.Np)))*rVariables.IntegrationCoefficient;

    //Distribute coupling block matrix into the elemental matrix
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    SizeType Index_i;

    for(SizeType i = 0; i<NumUNodes; i++)
    {
        Index_i = i * Dim;

        for(SizeType j = 0; j<NumPNodes; j++)
        {
            rLeftHandSideMatrix(Index_i,NumUNodes*Dim+j) -= CouplingMatrix(Index_i,j);
            rLeftHandSideMatrix(Index_i+1,NumUNodes*Dim+j) -= CouplingMatrix(Index_i+1,j);
            if(Dim > 2)
                rLeftHandSideMatrix(Index_i+2,NumUNodes*Dim+j) -= CouplingMatrix(Index_i+2,j);
        }
    }

    Matrix CouplingMatrixT = rVariables.NewmarkCoefficient1*trans(CouplingMatrix);

    //Distribute transposed coupling block matrix into the elemental matrix
    SizeType Index_j;
    for(SizeType i = 0; i<NumPNodes; i++)
    {
        for(SizeType j = 0; j<NumUNodes; j++)
        {
            Index_j = j * Dim;

            rLeftHandSideMatrix(NumUNodes*Dim+i,Index_j) += CouplingMatrixT(i,Index_j);
            rLeftHandSideMatrix(NumUNodes*Dim+i,Index_j+1) += CouplingMatrixT(i,Index_j+1);
            if(Dim > 2)
                rLeftHandSideMatrix(NumUNodes*Dim+i,Index_j+2) += CouplingMatrixT(i,Index_j+2);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    Matrix CompressibilityMatrix = rVariables.NewmarkCoefficient2*rVariables.BiotModulusInverse*outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;

    //Distribute compressibility block matrix into the elemental matrix
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for(SizeType i = 0; i < NumPNodes; i++)
    {
        for(SizeType j=0; j < NumPNodes; j++)
        {
            rLeftHandSideMatrix(NumUNodes*Dim+i,NumUNodes*Dim+j) += CompressibilityMatrix(i,j);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    Matrix PermeabilityMatrix = 1.0/rVariables.DynamicViscosity*
                                prod(rVariables.GradNpT,Matrix(prod(rVariables.IntrinsicPermeability,trans(rVariables.GradNpT))))*
                                rVariables.IntegrationCoefficient;

    //Distribute permeability block matrix into the elemental matrix
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for(SizeType i = 0; i < NumPNodes; i++)
    {
        for(SizeType j=0; j < NumPNodes; j++)
        {
            rLeftHandSideMatrix(NumUNodes*Dim+i,NumUNodes*Dim+j) += PermeabilityMatrix(i,j);
        }
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Vector StiffnessForce = prod(trans(rVariables.B), rVariables.StressVector)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block vector into the elemental vector
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();

    SizeType Index;

    for(SizeType i = 0; i < NumUNodes; i++)
    {
        Index = i * Dim;

        rRightHandSideVector[Index]   -= StiffnessForce[Index];
        rRightHandSideVector[Index+1] -= StiffnessForce[Index+1];
        if(Dim > 2)
            rRightHandSideVector[Index+2] -= StiffnessForce[Index+2];
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();

    double Porosity = GetProperties()[POROSITY];
    double Density = Porosity*GetProperties()[DENSITY_WATER] + (1.0-Porosity)*GetProperties()[DENSITY_SOLID];

    Vector BodyAcceleration = ZeroVector(Dim);
    SizeType Index = 0;
    for(SizeType i = 0; i < NumUNodes; i++)
    {
        BodyAcceleration[0] += rVariables.Nu[i]*rVariables.BodyAcceleration[Index++];
        BodyAcceleration[1] += rVariables.Nu[i]*rVariables.BodyAcceleration[Index++];
        if (Dim>2)
            BodyAcceleration[2] += rVariables.Nu[i]*rVariables.BodyAcceleration[Index++];
    }

    for(SizeType i=0; i < NumUNodes; i++)
    {
        Index = i * Dim;

        rRightHandSideVector[Index] += rVariables.Nu[i] * Density * BodyAcceleration[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Index+1] += rVariables.Nu[i] * Density * BodyAcceleration[1] * rVariables.IntegrationCoefficient;
        if(Dim > 2)
            rRightHandSideVector[Index+2] += rVariables.Nu[i] * Density * BodyAcceleration[2] * rVariables.IntegrationCoefficient;
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    unsigned int voigtsize  = 3;
    if( Dim == 3 ) voigtsize  = 6;

    Vector VoigtVector = ZeroVector(voigtsize);
    VoigtVector[0] = 1.0;
    VoigtVector[1] = 1.0;
    if(Dim == 3) VoigtVector[2] = 1.0;

    Matrix CouplingMatrix = rVariables.BiotCoefficient*prod(trans(rVariables.B),Matrix(outer_prod(VoigtVector,rVariables.Np)))*rVariables.IntegrationCoefficient;

    Vector CouplingForce = prod(CouplingMatrix,rVariables.PressureVector);

    //Distribute coupling block vector 1 into the elemental vector
    const SizeType NumUNodes = rGeom.PointsNumber();

    SizeType Index;

    for(SizeType i = 0; i<NumUNodes; i++)
    {
        Index = i * Dim;

        rRightHandSideVector[Index] += CouplingForce[Index];
        rRightHandSideVector[Index+1] += CouplingForce[Index+1];
        if(Dim > 2)
            rRightHandSideVector[Index+2] += CouplingForce[Index+2];
    }

    Vector CouplingFlow = prod(trans(CouplingMatrix),rVariables.VelocityVector);

    //Distribute coupling block vector 2 into the elemental vector
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for(SizeType i = 0; i<NumPNodes; i++)
    {
        rRightHandSideVector[NumUNodes*Dim+i] -= CouplingFlow[i];
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Matrix CompressibilityMatrix = rVariables.BiotModulusInverse*outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;

    Vector CompressibilityFlow = prod(CompressibilityMatrix,rVariables.PressureDtVector);

    //Distribute compressibility block vector into the elemental vector
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rRightHandSideVector[NumUNodes*Dim+i] -= CompressibilityFlow[i];
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Matrix PermeabilityMatrix = 1.0/rVariables.DynamicViscosity*prod(rVariables.GradNpT,Matrix(prod(rVariables.IntrinsicPermeability,trans(rVariables.GradNpT))))*rVariables.IntegrationCoefficient;

    Vector PermeabilityFlow = prod(PermeabilityMatrix,rVariables.PressureVector);

    //Distribute permeability block vector into the elemental vector
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rRightHandSideVector[NumUNodes*Dim+i] -= PermeabilityFlow[i];
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwDiffOrderElement::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Matrix GradNpTPerm = 1.0/rVariables.DynamicViscosity*GetProperties()[DENSITY_WATER]*
                         prod(rVariables.GradNpT,rVariables.IntrinsicPermeability)*rVariables.IntegrationCoefficient;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    Vector BodyAcceleration = ZeroVector(Dim);
    SizeType Index = 0;
    for(SizeType i = 0; i < NumUNodes; i++)
    {
        BodyAcceleration[0] += rVariables.Nu[i]*rVariables.BodyAcceleration[Index++];
        BodyAcceleration[1] += rVariables.Nu[i]*rVariables.BodyAcceleration[Index++];
        if (Dim>2)
            BodyAcceleration[2] += rVariables.Nu[i]*rVariables.BodyAcceleration[Index++];
    }

    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rRightHandSideVector[NumUNodes*Dim+i] += inner_prod(row(GradNpTPerm,i),BodyAcceleration);
    }
}

} // Namespace Kratos
