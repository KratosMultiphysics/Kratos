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
#include "custom_elements/one-phase_flow/small_strain_U_Pl_diff_order_element.hpp"

namespace Kratos
{

// Default Constructor
SmallStrainUPlDiffOrderElement::SmallStrainUPlDiffOrderElement() : Element() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SmallStrainUPlDiffOrderElement::SmallStrainUPlDiffOrderElement( IndexType NewId, GeometryType::Pointer pGeometry ) : Element( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SmallStrainUPlDiffOrderElement::SmallStrainUPlDiffOrderElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : Element( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = this->GetIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
SmallStrainUPlDiffOrderElement::~SmallStrainUPlDiffOrderElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Element::Pointer SmallStrainUPlDiffOrderElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallStrainUPlDiffOrderElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

int  SmallStrainUPlDiffOrderElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
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

    //Liquid variables
    if ( LIQUID_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "LIQUID_PRESSURE Key is 0. Check if all applications were correctly registered.", "" )

    if ( DT_LIQUID_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DT_LIQUID_PRESSURE Key is 0. Check if all applications were correctly registered.", "" )

    if ( DENSITY_LIQUID.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY_LIQUID Key is 0. Check if all applications were correctly registered.", "" )

    //verify that the dofs exist
    for ( unsigned int i = 0; i < rGeom.size(); i++ )
    {
        if ( rGeom[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", rGeom[i].Id() )

        if ( rGeom[i].HasDofFor( DISPLACEMENT_X ) == false || rGeom[i].HasDofFor( DISPLACEMENT_Y ) == false || rGeom[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", rGeom[i].Id() )


        if ( rGeom[i].SolutionStepsDataHas( LIQUID_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable LIQUID_PRESSURE on node ", rGeom[i].Id() )

        if ( rGeom[i].HasDofFor( LIQUID_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable LIQUID_PRESSURE on node ", rGeom[i].Id() )
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

void SmallStrainUPlDiffOrderElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& rGeom = GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != NumGPoints )
        mConstitutiveLawVector.resize( NumGPoints );

    //Imposed Z strain vector initialisation
    if ( mImposedZStrainVector.size() != NumGPoints )
        mImposedZStrainVector.resize( NumGPoints );

    if (Prop[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] =Prop[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( Prop, rGeom,row( rGeom.ShapeFunctionsValues( mThisIntegrationMethod ), i ) );

            mImposedZStrainVector[i] = 0.0;
        }
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element with ID ", this->Id() )


    const SizeType NumUNodes = rGeom.PointsNumber();

    switch(NumUNodes)
    {
        case 6: //2D T6P3
            mpPressureGeometry = GeometryType::Pointer( new Triangle2D3< Node >(rGeom(0), rGeom(1), rGeom(2)) );
            break;
        case 8: //2D Q8P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral2D4< Node >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 9: //2D Q9P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral2D4< Node >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 10: //3D T10P4
            mpPressureGeometry = GeometryType::Pointer( new Tetrahedra3D4< Node >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 20: //3D H20P8
            mpPressureGeometry = GeometryType::Pointer( new Hexahedra3D8< Node >(rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7)) );
            break;
        case 27: //3D H27P8
            mpPressureGeometry = GeometryType::Pointer( new Hexahedra3D8< Node >(rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7)) );
            break;
        default:
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected geometry type for different order interpolation element","");
            break;
    }

    // Initializing the intrinsic permeability matrix from the properties
    const SizeType Dim = rGeom.WorkingSpaceDimension();

    PoroElementUtilities::CalculatePermeabilityMatrix(mIntrinsicPermeability,Prop,Dim);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo ) const
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
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( LIQUID_PRESSURE );
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
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof( LIQUID_PRESSURE );

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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

void SmallStrainUPlDiffOrderElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR(std::logic_error,"SmallStrainUPlDiffOrderElement::CalculateLeftHandSide not implemented","");

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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

void SmallStrainUPlDiffOrderElement::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
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
    double Density = Porosity*GetProperties()[DENSITY_LIQUID] + (1.0-Porosity)*GetProperties()[DENSITY_SOLID];
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

void SmallStrainUPlDiffOrderElement::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
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
        rResult[Index++] = GetGeometry()[i].GetDof( LIQUID_PRESSURE ).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::GetSecondDerivativesVector( Vector& rValues, int Step ) const
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

void SmallStrainUPlDiffOrderElement::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    //Definition of variables
    ElementalVariables Variables;
    this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

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
            const double p0 = rGeom[0].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(LIQUID_PRESSURE);
            ThreadSafeNodeWrite(rGeom[3],LIQUID_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[4],LIQUID_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[5],LIQUID_PRESSURE, 0.5 * (p2 + p0) );
            break;
        }
        case 8: //2D Q8P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(LIQUID_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],LIQUID_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],LIQUID_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],LIQUID_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[7],LIQUID_PRESSURE, 0.5 * (p3 + p0) );
            break;
        }
        case 9: //2D Q9P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(LIQUID_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],LIQUID_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],LIQUID_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],LIQUID_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[7],LIQUID_PRESSURE, 0.5 * (p3 + p0) );
            ThreadSafeNodeWrite(rGeom[8],LIQUID_PRESSURE, 0.25 * (p0 + p1 + p2 + p3) );
            break;
        }
        case 10: //3D T10P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(LIQUID_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],LIQUID_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],LIQUID_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],LIQUID_PRESSURE, 0.5 * (p2 + p0) );
            ThreadSafeNodeWrite(rGeom[7],LIQUID_PRESSURE, 0.5 * (p0 + p3) );
            ThreadSafeNodeWrite(rGeom[8],LIQUID_PRESSURE, 0.5 * (p1 + p3) );
            ThreadSafeNodeWrite(rGeom[9],LIQUID_PRESSURE, 0.5 * (p2 + p3) );
            break;
        }
        case 20: //3D H20P8
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p4 = rGeom[4].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p5 = rGeom[5].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p6 = rGeom[6].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p7 = rGeom[7].FastGetSolutionStepValue(LIQUID_PRESSURE);
            // edges -- bottom
            ThreadSafeNodeWrite(rGeom[8],LIQUID_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[9],LIQUID_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[10],LIQUID_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[11],LIQUID_PRESSURE, 0.5 * (p3 + p0) );
            // edges -- middle
            ThreadSafeNodeWrite(rGeom[12],LIQUID_PRESSURE, 0.5 * (p4 + p0) );
            ThreadSafeNodeWrite(rGeom[13],LIQUID_PRESSURE, 0.5 * (p5 + p1) );
            ThreadSafeNodeWrite(rGeom[14],LIQUID_PRESSURE, 0.5 * (p6 + p2) );
            ThreadSafeNodeWrite(rGeom[15],LIQUID_PRESSURE, 0.5 * (p7 + p3) );
            // edges -- top
            ThreadSafeNodeWrite(rGeom[16],LIQUID_PRESSURE, 0.5 * (p4 + p5) );
            ThreadSafeNodeWrite(rGeom[17],LIQUID_PRESSURE, 0.5 * (p5 + p6) );
            ThreadSafeNodeWrite(rGeom[18],LIQUID_PRESSURE, 0.5 * (p6 + p7) );
            ThreadSafeNodeWrite(rGeom[19],LIQUID_PRESSURE, 0.5 * (p7 + p0) );
            break;
        }
        case 27: //3D H27P8
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p4 = rGeom[4].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p5 = rGeom[5].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p6 = rGeom[6].FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double p7 = rGeom[7].FastGetSolutionStepValue(LIQUID_PRESSURE);
            // edges -- bottom
            ThreadSafeNodeWrite(rGeom[8],LIQUID_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[9],LIQUID_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[10],LIQUID_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[11],LIQUID_PRESSURE, 0.5 * (p3 + p0) );
            // edges -- middle
            ThreadSafeNodeWrite(rGeom[12],LIQUID_PRESSURE, 0.5 * (p4 + p0) );
            ThreadSafeNodeWrite(rGeom[13],LIQUID_PRESSURE, 0.5 * (p5 + p1) );
            ThreadSafeNodeWrite(rGeom[14],LIQUID_PRESSURE, 0.5 * (p6 + p2) );
            ThreadSafeNodeWrite(rGeom[15],LIQUID_PRESSURE, 0.5 * (p7 + p3) );
            // edges -- top
            ThreadSafeNodeWrite(rGeom[16],LIQUID_PRESSURE, 0.5 * (p4 + p5) );
            ThreadSafeNodeWrite(rGeom[17],LIQUID_PRESSURE, 0.5 * (p5 + p6) );
            ThreadSafeNodeWrite(rGeom[18],LIQUID_PRESSURE, 0.5 * (p6 + p7) );
            ThreadSafeNodeWrite(rGeom[19],LIQUID_PRESSURE, 0.5 * (p7 + p0) );
            // face centers
            ThreadSafeNodeWrite(rGeom[20],LIQUID_PRESSURE, 0.25 * (p0 + p1 + p2 + p3) );
            ThreadSafeNodeWrite(rGeom[21],LIQUID_PRESSURE, 0.25 * (p0 + p1 + p4 + p5) );
            ThreadSafeNodeWrite(rGeom[22],LIQUID_PRESSURE, 0.25 * (p1 + p2 + p5 + p6) );
            ThreadSafeNodeWrite(rGeom[23],LIQUID_PRESSURE, 0.25 * (p2 + p3 + p6 + p7) );
            ThreadSafeNodeWrite(rGeom[24],LIQUID_PRESSURE, 0.25 * (p3 + p0 + p7 + p4) );
            ThreadSafeNodeWrite(rGeom[25],LIQUID_PRESSURE, 0.25 * (p4 + p5 + p6 + p7) );
            // element center
            ThreadSafeNodeWrite(rGeom[26],LIQUID_PRESSURE, 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7) );
            break;
        }
        default:
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected geometry type for different order interpolation element","");
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::SetValuesOnIntegrationPoints( const Variable<double>& rVariable,const std::vector<double>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if (rVariable == IMPOSED_Z_STRAIN_VALUE) {
        for ( IndexType PointNumber = 0; PointNumber < mImposedZStrainVector.size(); ++PointNumber ) {
            mImposedZStrainVector[PointNumber] = rValues[PointNumber];
        }

    } else {
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable,const std::vector<Vector>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable,const std::vector<Matrix>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if (rVariable == PERMEABILITY_MATRIX) {
        // Permeability is set only on the element, not on every GP
        noalias(mIntrinsicPermeability) = rValues[0];

    } else {
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number, false );

    if ( rVariable == VON_MISES_STRESS ) {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

        //Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

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
    } else {
        for ( unsigned int i = 0; i < integration_points_number; i++ ) {
            rOutput[i] = 0.0;
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_VECTOR ) {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

        //Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

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
    } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ) {
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
    } else {
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ ) {
            if ( rOutput[i].size() != dimension )
                rOutput[i].resize( dimension, false );
            noalias(rOutput[i]) = ZeroVector(dimension);
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable , rOutput[i] );
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == LIQUID_FLUX_VECTOR ) {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            //Compute LiquidFlux vector q [m/s]
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
            noalias(GradPressureTerm) -= GetProperties()[DENSITY_LIQUID]*BodyAcceleration;

            Vector AuxLiquidFlux = ZeroVector(Dim);
            AuxLiquidFlux = - 1.0/Variables.DynamicViscosity * prod(mIntrinsicPermeability, GradPressureTerm );

            array_1d<double,3> LiquidFlux = ZeroVector(3);
            LiquidFlux[0] = AuxLiquidFlux[0];
            LiquidFlux[1] = AuxLiquidFlux[1];
            if(Dim>2)
                LiquidFlux[2] = AuxLiquidFlux[2];

            rOutput[PointNumber] = LiquidFlux;
        }
    } else if ( rVariable == LIQUID_PRESSURE_GRADIENT ) {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            //Compute LiquidFlux vector q [m/s]
            const SizeType Dim = rGeom.WorkingSpaceDimension();

            Vector GradPressure(Dim);
            noalias(GradPressure) = prod(trans(Variables.GradNpT),Variables.PressureVector);

            array_1d<double,3> GradPressureVector = ZeroVector(3);
            GradPressureVector[0] = GradPressure[0];
            GradPressureVector[1] = GradPressure[1];
            if(Dim>2)
                GradPressureVector[2] = GradPressure[2];

            rOutput[PointNumber] = GradPressureVector;
        }
    } else {
        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); i++ )
        {
            noalias(rOutput[i]) = ZeroVector(3);
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    const unsigned int dimension       = rGeom.WorkingSpaceDimension();
    const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == EFFECTIVE_STRESS_TENSOR ) {
        std::vector<Vector> StressVector;

        this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != cl_dimension )
                rOutput[PointNumber].resize( cl_dimension, cl_dimension, false );

            noalias(rOutput[PointNumber]) = MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);
        }
    } else if ( rVariable == TOTAL_STRESS_TENSOR ) {
        std::vector<Vector> StressVector;
        this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );

        double Pressure;
        const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
        Vector PressureVector(NumPNodes);
        for(SizeType i=0; i<NumPNodes; i++)
        {
            PressureVector[i] = rGeom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
        }
        const double BiotCoefficient = this->GetProperties()[BIOT_COEFFICIENT];
        Matrix NpContainer(integration_points_number,NumPNodes);
        noalias(NpContainer) = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

        const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
        Vector VoigtVector(strain_size);
        noalias(VoigtVector) = ZeroVector(strain_size);
        if(cl_dimension == 3) {
            VoigtVector[0] = 1.0;
            VoigtVector[1] = 1.0;
            VoigtVector[2] = 1.0;
        } else {
            VoigtVector[0] = 1.0;
            VoigtVector[1] = 1.0;
        }

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            Pressure = 0.0;
            for(unsigned int i = 0; i < NumPNodes; i++)
            {
                Pressure += NpContainer(PointNumber,i)*PressureVector[i];
            }

            noalias(StressVector[PointNumber]) += -BiotCoefficient*Pressure*VoigtVector;

            if ( rOutput[PointNumber].size2() != cl_dimension )
                rOutput[PointNumber].resize( cl_dimension, cl_dimension, false );

            noalias(rOutput[PointNumber]) = MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);
        }
    } else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ) {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != cl_dimension )
                rOutput[PointNumber].resize( cl_dimension, cl_dimension, false );

            noalias(rOutput[PointNumber]) = MathUtils<double>::StrainVectorToTensor(StrainVector[PointNumber]);
        }
    } else if(rVariable == PERMEABILITY_MATRIX) {
        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            noalias(rOutput[PointNumber]) = mIntrinsicPermeability;
        }
    } else {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ ){
            if ( rOutput[i].size2() != dimension )
                rOutput[i].resize( dimension, dimension, false );
            noalias(rOutput[i]) = ZeroMatrix(dimension, dimension);
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable , rOutput[i] );
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,std::vector<ConstitutiveLaw::Pointer>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const unsigned int integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (unsigned int point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo,
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
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if(CalculateResidualVectorFlag)
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

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

void SmallStrainUPlDiffOrderElement::InitializeElementalVariables (ElementalVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

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
    (rVariables.B).resize(strain_size, NumUNodes * Dim, false);
    noalias(rVariables.B) = ZeroMatrix( strain_size, NumUNodes * Dim );
    (rVariables.StrainVector).resize(strain_size,false);
    (rVariables.ConstitutiveMatrix).resize(strain_size, strain_size, false);
    (rVariables.StressVector).resize(strain_size,false);

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
    rVariables.NewmarkCoefficient2 = rCurrentProcessInfo[DT_LIQUID_PRESSURE_COEFFICIENT];
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::InitializeNodalVariables (ElementalVariables& rVariables)
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
        rVariables.PressureVector[i]   = rGeom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
        rVariables.PressureDtVector[i] = rGeom[i].FastGetSolutionStepValue(DT_LIQUID_PRESSURE);
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::InitializeProperties (ElementalVariables& rVariables)
{
    double BulkModulusSolid = GetProperties()[BULK_MODULUS_SOLID];
    rVariables.BiotCoefficient = GetProperties()[BIOT_COEFFICIENT];
    double Porosity = GetProperties()[POROSITY];
    rVariables.BiotModulusInverse = (rVariables.BiotCoefficient-Porosity)/BulkModulusSolid + Porosity/GetProperties()[BULK_MODULUS_LIQUID];
    rVariables.DynamicViscosity = GetProperties()[DYNAMIC_VISCOSITY_LIQUID];
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateKinematics(ElementalVariables& rVariables, unsigned int PointNumber)

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

    const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();
    // 2.5D element (2D Geometry with 3D ConstitutiveLaw)
    if (cl_dimension > Dim) {

        // StrainVector must have the shape of a 3D element
        rVariables.StrainVector[3] = rVariables.StrainVector[2];
        rVariables.StrainVector[2] = mImposedZStrainVector[PointNumber];

        // B matrix must have the shape of a 3D element
        for ( SizeType i = 0; i < NumUNodes; i++ )
        {
            node = 2 * i;

            rVariables.B( 3, node + 0 ) = rVariables.B( 2, node + 0 );
            rVariables.B( 3, node + 1 ) = rVariables.B( 2, node + 1 );
            rVariables.B( 2, node + 0 ) = 0.0;
            rVariables.B( 2, node + 1 ) = 0.0;
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::SetElementalVariables(ElementalVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters)
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

void SmallStrainUPlDiffOrderElement::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, double detJu, double weight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rIntegrationCoefficient = detJu * weight;

    if( dimension == 2 )
        rIntegrationCoefficient *= GetProperties()[THICKNESS];
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
}


//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
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

void SmallStrainUPlDiffOrderElement::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();

    Vector VoigtVector = ZeroVector(strain_size);
    VoigtVector[0] = 1.0;
    VoigtVector[1] = 1.0;
    if(cl_dimension == 3) VoigtVector[2] = 1.0;

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

void SmallStrainUPlDiffOrderElement::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
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

void SmallStrainUPlDiffOrderElement::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    Matrix PermeabilityMatrix = 1.0/rVariables.DynamicViscosity*
                                prod(rVariables.GradNpT,Matrix(prod(mIntrinsicPermeability,trans(rVariables.GradNpT))))*
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

void SmallStrainUPlDiffOrderElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPlDiffOrderElement::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
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

void SmallStrainUPlDiffOrderElement::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();

    double Porosity = GetProperties()[POROSITY];
    double Density = Porosity*GetProperties()[DENSITY_LIQUID] + (1.0-Porosity)*GetProperties()[DENSITY_SOLID];

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

void SmallStrainUPlDiffOrderElement::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();

    Vector VoigtVector = ZeroVector(strain_size);
    VoigtVector[0] = 1.0;
    VoigtVector[1] = 1.0;
    if(cl_dimension == 3) VoigtVector[2] = 1.0;

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

void SmallStrainUPlDiffOrderElement::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
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

void SmallStrainUPlDiffOrderElement::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Matrix PermeabilityMatrix = 1.0/rVariables.DynamicViscosity*prod(rVariables.GradNpT,Matrix(prod(mIntrinsicPermeability,trans(rVariables.GradNpT))))*rVariables.IntegrationCoefficient;

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

void SmallStrainUPlDiffOrderElement::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Matrix GradNpTPerm = 1.0/rVariables.DynamicViscosity*GetProperties()[DENSITY_LIQUID]*
                         prod(rVariables.GradNpT,mIntrinsicPermeability)*rVariables.IntegrationCoefficient;

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
