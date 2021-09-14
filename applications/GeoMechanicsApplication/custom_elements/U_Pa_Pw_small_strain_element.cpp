//
//   Project Name:        Kratos GeoMechanics
//   Author:              Hoang-Giang Bui, Felix Nagel
//   Date:                31 Mar 2020
//   Revision:            1.1
//
//

// System includes
#include <typeinfo>


// External includes


// Project includes
#include "custom_retention/retention_law_factory.h"
#include "custom_elements/U_Pa_Pw_small_strain_element.h"
#include "geo_mechanics_application_variables.h"

//#define ENABLE_DEBUG_CONSTITUTIVE_LAW
#define CHECK_NAN
// #define DEBUG_AIR

namespace Kratos
{

UPaPwSmallStrainElement::UPaPwSmallStrainElement( IndexType NewId,
        GeometryType::Pointer pGeometry )
    : BaseType( NewId, pGeometry )
{
    this->InitializeSubGeometry();
}

UPaPwSmallStrainElement::UPaPwSmallStrainElement( IndexType NewId,
        GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseType( NewId, pGeometry, pProperties )
{
    this->InitializeSubGeometry();
}

//************************************************************************************
//************************************************************************************

Element::Pointer UPaPwSmallStrainElement::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UPaPwSmallStrainElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer UPaPwSmallStrainElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UPaPwSmallStrainElement( NewId, pGeom, pProperties ) );
}

//************************************************************************************
//************************************************************************************

UPaPwSmallStrainElement::~UPaPwSmallStrainElement()
{
}

//************************************************************************************
//************************************************************************************
void UPaPwSmallStrainElement::InitializeSubGeometry()
{
    // #ifdef ENABLE_BEZIER_GEOMETRY
    // // account for Bezier geometry
    // if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier2D
    //   || GetGeometry().GetGeometryType() == GeometryData::Kratos_Bezier3D )
    // {
    //     mpSubGeometry = this->pGetGeometry();
    //     mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    //     mIsStabilized = false;
    //     return;
    // }
    // #endif

    if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D27 )
    {
        mpSubGeometry = GeometryType::Pointer( new Hexahedra3D8 <Node<3> >(
                                 GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ),
                                 GetGeometry()( 4 ), GetGeometry()( 5 ), GetGeometry()( 6 ), GetGeometry()( 7 ) ) );
        mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        mIsStabilized = false;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D20 ) //remarks: this element does not work correctly with the 20 nodes discretisation. See the cube consolidation test
    {
        mpSubGeometry = GeometryType::Pointer( new Hexahedra3D8 <Node<3> >(
                                 GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ),
                                 GetGeometry()( 4 ), GetGeometry()( 5 ), GetGeometry()( 6 ), GetGeometry()( 7 ) ) );
        mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        mIsStabilized = false;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Tetrahedra3D10 )
    {
        mpSubGeometry = GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ) ) );
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        mIsStabilized = false;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Prism3D15 )
    {
        mpSubGeometry = GeometryType::Pointer( new Prism3D6 <Node<3> >( GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ), GetGeometry()( 4 ), GetGeometry()( 5 ) ) );
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        mIsStabilized = false;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
           || GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8 )
    {
        mpSubGeometry = GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >(
                                 GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ) ) );
        mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        mIsStabilized = false;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Triangle2D6 )
    {
        mpSubGeometry = GeometryType::Pointer( new Triangle2D3 <Node<3> >(
                                 GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ) ) );
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        mIsStabilized = false;
    }
    /**** low order element, stabilisation activated ****/
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Hexahedra3D8 )
    {
        mpSubGeometry = this->pGetGeometry();
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2; //remarks: GI_GAUSS_2 gives better result than GI_GAUSS_3 for tunnel problem
        mIsStabilized = true;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Tetrahedra3D4 )
    {
        mpSubGeometry = this->pGetGeometry();
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        mIsStabilized = true;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4 )
    {
        mpSubGeometry = this->pGetGeometry();
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        mIsStabilized = true;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Triangle2D3 )
    {
        mpSubGeometry = this->pGetGeometry();
        mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        mIsStabilized = true;
    }
    else if ( GetGeometry().GetGeometryType() == GeometryData::Kratos_Prism3D6 )
    {
        mpSubGeometry = this->pGetGeometry();
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        mIsStabilized = true;
    }
    else
    {
        KRATOS_THROW_ERROR( std::logic_error, "This element matches only with a quadratic hexahedra (8, 20 or 27), tetrahedra (4, 10) or prism (15) geometry" , *this );
    }
}

//************************************************************************************
//************************************************************************************
void UPaPwSmallStrainElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if (rCurrentProcessInfo[RESET_CONFIGURATION] == 0)
    {
        // integration rule
        if(this->Has( INTEGRATION_ORDER ))
        {
            if(this->GetValue(INTEGRATION_ORDER) == 1)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 2)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 3)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 4)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
            }
            else if(this->GetValue(INTEGRATION_ORDER) == 5)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "SoilsElement does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
        }
        else if(GetProperties().Has( INTEGRATION_ORDER ))
        {
            if(GetProperties()[INTEGRATION_ORDER] == 1)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 2)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 3)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 4)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 5)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "SoilsElement does not support for integration rule", GetProperties()[INTEGRATION_ORDER])
        }
        else
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

        // number of integration points used, mThisIntegrationMethod refers to the
        // integration method defined in the constructor
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        // constitutive Law initialization
        mConstitutiveLawVector.resize( integration_points.size() );
        mRetentionLawVector.resize( integration_points.size() );
        mReferencePressures.resize( integration_points.size() );
        this->InitializeMaterial(rCurrentProcessInfo);

        // initialize zero displacement
        mInitialDisp.resize( GetGeometry().size(), dim, false );
        noalias(mInitialDisp) = ZeroMatrix(GetGeometry().size(), dim);

        // initializing the Jacobian in the reference configuration
        GeometryType::JacobiansType J0;
        Matrix DeltaPosition(GetGeometry().size(), 3);

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        // calculating the domain size
        double TotalDomainInitialSize = 0.0;
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            //getting informations for integration
            double IntegrationWeight = integration_points[PointNumber].Weight();
            //calculating the total area/volume
            TotalDomainInitialSize += MathUtils<double>::Det(J0[PointNumber]) * IntegrationWeight;
        }
        // extern Variable<double> GEOMETRICAL_DOMAIN_SIZE;
        // this->SetValue(GEOMETRICAL_DOMAIN_SIZE, TotalDomainInitialSize);

        if (TotalDomainInitialSize < 0)
        {
            std::cout << "error on element -> " << this->Id() << std::endl;
            std::cout << "TotalDomainInitialSize = " << TotalDomainInitialSize << std::endl;
            std::cout << "Properties " << GetProperties().Id() << ": " << GetProperties() << std::endl;
            std::cout << "Element nodes:" << std::endl;
            for (std::size_t i = 0; i < GetGeometry().size(); ++i)
                std::cout << " " << GetGeometry()[i].Id() << ": " << GetGeometry()[i].GetInitialPosition() << std::endl;
            std::cout << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Domain size can not be less than 0. Please check Jacobian.", __FUNCTION__);
        }
    }
    else if (rCurrentProcessInfo[RESET_CONFIGURATION] == 1)
    {
        // setup initial displacement for stress-free activation of elements
        if (mInitialDisp.size1() != GetGeometry().size() || mInitialDisp.size2() != dim)
            mInitialDisp.resize( GetGeometry().size(), dim, false );
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            for ( unsigned int i = 0; i < 3; ++i )
                mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];
    }

    this->InitializeSubGeometry();

    KRATOS_CATCH( "" )
}

void UPaPwSmallStrainElement::ResetConstitutiveLaw()
{
    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
    {
        Vector dummy;
        // dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), PointNumber );
        mConstitutiveLawVector[PointNumber]->ResetMaterial( GetProperties(), GetGeometry(),  dummy);
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif
}

void UPaPwSmallStrainElement::InitializeSolutionStep( const ProcessInfo& CurrentProcessInfo )
{
    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
    {
        Vector dummy;
        // dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), PointNumber );
        mConstitutiveLawVector[Point]->InitializeSolutionStep( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif
}

void UPaPwSmallStrainElement::InitializeNonLinearIteration( const ProcessInfo& CurrentProcessInfo )
{
    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
    {
        Vector dummy;
        // dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), PointNumber );
        mConstitutiveLawVector[Point]->InitializeNonLinearIteration( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif
}

void UPaPwSmallStrainElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

void UPaPwSmallStrainElement::CalculateRightHandSide( VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
}

void UPaPwSmallStrainElement::FinalizeNonLinearIteration( const ProcessInfo& CurrentProcessInfo )
{
    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
    {
        Vector dummy;
        // dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), PointNumber );
        mConstitutiveLawVector[Point]->FinalizeNonLinearIteration( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif
}

void UPaPwSmallStrainElement::FinalizeSolutionStep( const ProcessInfo& CurrentProcessInfo )
{
    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
    {
        Vector dummy;
        // dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), PointNumber );
        mConstitutiveLawVector[Point]->FinalizeSolutionStep( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif

    const double TOL = 1.0e-10;

    // check the porosity
    std::vector<double> porosities;
    this->CalculateOnIntegrationPoints(POROSITY, porosities, CurrentProcessInfo);
    for(unsigned int Point = 0; Point < porosities.size(); ++Point)
    {
        if(porosities[Point] < -TOL || porosities[Point] > 1.0+TOL)
        {
            std::cout << "porosities at element " << Id() << ":" << std::endl;
            for (std::size_t i = 0; i < porosities.size(); ++i)
                std::cout << " " << porosities[i] << std::endl;
            std::cout << std::endl;

            std::cout << "displacement at element " << Id() << ":" << std::endl;
            for (std::size_t i = 0; i < GetGeometry().size(); ++i)
                std::cout << " " << GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT) << std::endl;
            std::cout << std::endl;

            std::stringstream ss;
            ss << "The porosity is not in range [0.0, 1.0] at element " << Id() << ", point " << Point
               << ". It is " << porosities[Point];

            if (GetValue(USE_DISTRIBUTED_PROPERTIES) == false)
                ss << ". FIX_POROSITY = " << GetProperties()[FIX_POROSITY];
            else
                ss << ". FIX_POROSITY = " << GetValue(FIX_POROSITY);

            ss << ". Properties Id = " << GetProperties().Id() << "." << std::endl;
            ss << GetProperties();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), __FUNCTION__)
        }
    }

    // check the saturation
    std::vector<double> saturations;
    this->CalculateOnIntegrationPoints(SATURATION, saturations, CurrentProcessInfo);
    for(unsigned int Point = 0; Point < saturations.size(); ++Point)
    {
        if(saturations[Point] < -TOL || saturations[Point] > 1.0+TOL)
        {
            std::stringstream ss;
            ss << "The saturation is not in range [0.0, 1.0] at element " << Id()
               << ", point " << Point << ". It is " << saturations[Point];
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), __FUNCTION__)
        }
    }
}

void UPaPwSmallStrainElement::InitializeMaterial(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
    }

    for ( unsigned int i = 0; i < mRetentionLawVector.size(); i++ )
    {
        mRetentionLawVector[i] = RetentionLawFactory::Clone(GetProperties());
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo);
    }

    for ( unsigned int i = 0; i < mRetentionLawVector.size(); ++i )
    {
        mRetentionLawVector[i]-> InitializeMaterial(  GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    GetGeometry().Clean();
    #endif

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UPaPwSmallStrainElement::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
{
    unsigned int dim_press = 2;//two pressure dofs
    unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
    unsigned int MatSize = mpSubGeometry->size() * dim_press + GetGeometry().size() * dim_disp;

    if( rResult.size() != MatSize )
        rResult.resize( MatSize, false );

    unsigned int cnt = 0;

    for ( unsigned int i = 0; i < mpSubGeometry->size(); i++ )
    {
        rResult[cnt++] = (*mpSubGeometry)[i].GetDof( WATER_PRESSURE ).EquationId();
        rResult[cnt++] = (*mpSubGeometry)[i].GetDof( AIR_PRESSURE ).EquationId();
    }

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rResult[cnt++] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[cnt++] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[cnt++] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }
}

void UPaPwSmallStrainElement::GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
{
    unsigned int dim_press = 2;//two pressure dofs
    unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
    unsigned int MatSize = mpSubGeometry->size() * dim_press + GetGeometry().size() * dim_disp;

    if( ElementalDofList.size() != MatSize )
        ElementalDofList.resize( MatSize );

    unsigned int cnt = 0;

    for ( unsigned int i = 0; i < mpSubGeometry->size(); i++ )
    {
        ElementalDofList[cnt++] = (*mpSubGeometry)[i].pGetDof( WATER_PRESSURE );
        ElementalDofList[cnt++] = (*mpSubGeometry)[i].pGetDof( AIR_PRESSURE );
    }

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        ElementalDofList[cnt++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        ElementalDofList[cnt++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        ElementalDofList[cnt++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
    }
}

void UPaPwSmallStrainElement::GetValuesVector( Vector& values, int Step ) const
{
    unsigned int dim_press = 2;//one pressure dofs
    unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
    unsigned int MatSize = mpSubGeometry->size() * dim_press + GetGeometry().size() * dim_disp;

    if ( values.size() != MatSize )
        values.resize( MatSize, false );

    unsigned int cnt = 0;

    for ( unsigned int i = 0; i < mpSubGeometry->size(); i++ )
    {
        values( cnt++ ) = (*mpSubGeometry)[i].GetSolutionStepValue( WATER_PRESSURE, Step );
        values( cnt++ ) = (*mpSubGeometry)[i].GetSolutionStepValue( AIR_PRESSURE, Step );
    }

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
        values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}

void UPaPwSmallStrainElement::GetFirstDerivativesVector( Vector& values, int Step ) const
{
    unsigned int dim_press = 1;//one pressure dofs
    unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
    unsigned int MatSize = mpSubGeometry->size() * dim_press + GetGeometry().size() * dim_disp;

    if ( values.size() != MatSize )
        values.resize( MatSize, false );

    unsigned int cnt = 0;

    for ( unsigned int i = 0; i < mpSubGeometry->size(); i++ )
    {
        values( cnt++ ) = (*mpSubGeometry)[i].GetSolutionStepValue( WATER_PRESSURE_DT, Step );
        values( cnt++ ) = (*mpSubGeometry)[i].GetSolutionStepValue( AIR_PRESSURE_DT, Step );
    }

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
        if (dim_disp == 3)
            values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
    }
}

void UPaPwSmallStrainElement::GetSecondDerivativesVector( Vector& values, int Step ) const
{
    unsigned int dim_press = 1;//one pressure dofs
    unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
    unsigned int MatSize = mpSubGeometry->size() * dim_press + GetGeometry().size() * dim_disp;

    if ( values.size() != MatSize )
        values.resize( MatSize, false );

    unsigned int cnt = 0;

    for ( unsigned int i = 0; i < mpSubGeometry->size(); i++ )
    {
        values( cnt++ ) = (*mpSubGeometry)[i].GetSolutionStepValue( WATER_PRESSURE_ACCELERATION, Step );
        values( cnt++ ) = (*mpSubGeometry)[i].GetSolutionStepValue( AIR_PRESSURE_ACCELERATION, Step );
    }

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
        if (dim_disp == 3)
            values( cnt++ ) = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }
}

void UPaPwSmallStrainElement::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int number_of_nodes_disp = GetGeometry().size();

    unsigned int number_of_nodes_press = mpSubGeometry->size();

    unsigned int number_of_nodes_air = mpSubGeometry->size();

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int strain_size = dim * (dim + 1) / 2;

    unsigned int MatSizeU = number_of_nodes_disp * dim;

    unsigned int MatSizeW = number_of_nodes_press;

    unsigned int MatSizeA = number_of_nodes_air;

    unsigned int MatSize = MatSizeU + MatSizeW + MatSizeA;

    //resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize || rLeftHandSideMatrix.size2() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    #ifdef ENABLE_BEZIER_GEOMETRY
    // initialize the geometry
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    #ifdef ENABLE_BEZIER_GEOMETRY
    mpSubGeometry->Initialize(integration_points); // this is to make sure that the Bezier geometry of the pressure will have the same number of integration points as the geometry for displacement
    #endif

    const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =
        GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    const GeometryType::ShapeFunctionsGradientsType& DN_De_Pressure =
        mpSubGeometry->ShapeFunctionsLocalGradients( mThisIntegrationMethod ); // in the case of Bezier geometry, mThisIntegrationMethod does not play a role

    const Matrix& Ncontainer_Displacement = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    const Matrix& Ncontainer_Pressure = mpSubGeometry->ShapeFunctionsValues( mThisIntegrationMethod ); // in the case of Bezier geometry, mThisIntegrationMethod does not play a role

    double Weight;

    double capillaryPressure;

    double waterPressure;

    double airPressure;

    double waterPressure_Dt;

    double airPressure_Dt;

    double porosity;

    double Dn_DdivU;

    double density;

    double density_soil;

    double permeability_water;

    double density_water;

    double saturation;

    double DS_Dpc;

    double Dpc_Dt;

    double divU_Dt;

    double Drho_DdivU;

    double Drho_Dpw;

    double Drho_Dpa;

    Vector grad_water( dim );

    Vector flow_water( dim );

    double permeability_air;

    double density_air;

    double bulk_air;

    Vector grad_air( dim );

    Vector flow_air( dim );

    double D2S_Dpc2;

    Vector DflowWater_Dpw( dim );

    double DflowWater_Dgradpw;

    Vector DflowWater_Dpa( dim );

    Vector DflowAir_Dpw( dim );

    double DflowAir_Dgradpa;

    Vector DflowAir_Dpa( dim );

    Matrix InvJ(dim, dim);

    double DetJ = 0.0;

    Matrix Help_K_UU;

    Matrix Help_K_UW;

    Matrix Help_K_UA;

    Matrix Help_K_WU;

    Matrix Help_K_WW;

    Matrix Help_K_WA;

    Matrix Help_K_AU;

    Matrix Help_K_AW;

    Matrix Help_K_AA;

    Vector Help_R_U;

    Vector Help_R_W;

    Vector Help_R_A;

    Matrix DN_DX_PRESS( number_of_nodes_press, dim );

    Vector N_DISP( number_of_nodes_disp );

    Vector N_PRESS( number_of_nodes_press );

    Matrix TanC_W( dim, dim );

    Matrix TanC_A( dim, dim );

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        Help_K_UU.resize( MatSizeU, MatSizeU, false );
        Help_K_UW.resize( MatSizeU, MatSizeW, false );
        Help_K_UA.resize( MatSizeU, MatSizeA, false );

        Help_K_WU.resize( MatSizeW, MatSizeU, false );
        Help_K_WW.resize( MatSizeW, MatSizeW, false );
        Help_K_WA.resize( MatSizeA, MatSizeW, false );

        Help_K_AU.resize( MatSizeA, MatSizeU, false );
        Help_K_AW.resize( MatSizeA, MatSizeW, false );
        Help_K_AA.resize( MatSizeA, MatSizeA, false );
    }

    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        Help_R_U.resize( MatSizeU, false );
        Help_R_W.resize( MatSizeW, false );
        Help_R_A.resize( MatSizeA, false );
    }

    //Initialize local variables
    Matrix B( strain_size, MatSizeU );

    Matrix TanC_U( strain_size, strain_size );

    Vector StrainVector( strain_size );

    Vector StressVector( strain_size );

    Matrix DN_DX_DISP( number_of_nodes_disp, dim );

    Matrix CurrentDisp( number_of_nodes_disp, dim );

    //Current displacements
    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
    {
        CurrentDisp( node, 0 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_X );
        CurrentDisp( node, 1 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_Y );
        CurrentDisp( node, 2 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_Z );
    }

    //initializing the Jacobian in the reference configuration
    GeometryType::JacobiansType J0;
    Matrix DeltaPosition(GetGeometry().size(), 3);

    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
    {
        noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
    }

    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

    // if (Id() == 1)
    // {
    //     for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
    //     {
    //         std::cout << "node " << GetGeometry()[node].Id() << " displacement:"
    //                   << " " << GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_X)
    //                   << " " << GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_Y)
    //                   << " " << GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_Z)
    //                   << ", wp: " << GetGeometry()[node].GetSolutionStepValue(WATER_PRESSURE)
    //                   << ", ap: " << GetGeometry()[node].GetSolutionStepValue(AIR_PRESSURE)
    //                   << std::endl;
    //     }
    // }

//    /////////////////////////////////////////////////////////////////////////
//    //// Integration in space to compute the average of pressure shape function
//    /////////////////////////////////////////////////////////////////////////
//    Vector N_PRESS_averaged;
    // extern Variable<double> GEOMETRICAL_DOMAIN_SIZE;
    // double TotalDomainInitialSize = this->GetValue(GEOMETRICAL_DOMAIN_SIZE);
//    if(mIsStabilized)
//    {
//        N_PRESS_averaged.resize(number_of_nodes_press, false);
//        N_PRESS_averaged = ZeroVector( number_of_nodes_press );
//        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
//        {
//            Weight = integration_points[PointNumber].Weight();
//            DetJ   = MathUtils<double>::Det(J0[PointNumber]);

//            // Shape Functions on current spatial quadrature point
//            noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

//            noalias( N_PRESS_averaged ) += N_PRESS * Weight * DetJ;
//        }
//        N_PRESS_averaged /= TotalDomainInitialSize;
//    }

//    /////////////////////////////////////////////////////////////////////////
//    //// Integration in space to compute the average of capillaryPressure_Dt
//    /////////////////////////////////////////////////////////////////////////
//    double averageCapillaryPressure_dt = 0.0;
//    if(mIsStabilized)
//    {
//        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
//        {
//            Weight = integration_points[PointNumber].Weight();
//            DetJ   = MathUtils<double>::Det(J0[PointNumber]);

//            // Shape Functions on current spatial quadrature point
//            noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );
//            capillaryPressure_Dt = -GetDerivativeDWaterPressureDt(N_PRESS);

//            averageCapillaryPressure_dt += capillaryPressure_Dt * Weight * DetJ;
//        }
//        averageCapillaryPressure_dt /= TotalDomainInitialSize;
//    }

    if ( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
    {
        density_soil = GetValue( DENSITY );
    }
    else
    {
        density_soil = GetProperties()[DENSITY];
    }

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    {
        MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ, DetJ );

        noalias( DN_DX_PRESS ) = prod( DN_De_Pressure[PointNumber], InvJ );

        noalias( DN_DX_DISP ) = prod( DN_De_Displacement[PointNumber], InvJ );

        Weight = integration_points[PointNumber].Weight();

        //modify integration weight in case of 2D
        if ( dim == 2 ) Weight *= GetProperties()[THICKNESS];

        // Shape Functions on current spatial quadrature point
        noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

        noalias( N_DISP ) = row( Ncontainer_Displacement, PointNumber );

        //Initializing B_Operator at the current integration point
        CalculateBoperator( B, DN_DX_DISP );

        //Calculate the current strain vector using the B-Operator
        CalculateStrain( B, CurrentDisp, StrainVector );

        waterPressure = GetWaterPressure( N_PRESS );

        airPressure = GetAirPressure( N_PRESS );

        RetentionParameters.SetFluidPressure( waterPressure );

        RetentionParameters.SetAirPressure( airPressure );

        waterPressure_Dt = GetDerivativeDWaterPressureDt( N_PRESS );

        airPressure_Dt = GetDerivativeDAirPressureDt( N_PRESS );

        Dpc_Dt = airPressure_Dt - waterPressure_Dt;

        saturation = GetSaturation( PointNumber, RetentionParameters );

        DS_Dpc = GetDerivativeDSaturationDpc( PointNumber, RetentionParameters );

        D2S_Dpc2 = GetSecondDerivativeD2SaturationDpc2( PointNumber, RetentionParameters );

        porosity = GetPorosity( DN_DX_DISP );

        Dn_DdivU = GetDerivativeDPorosityDDivU( DN_DX_DISP );

        divU_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

        this->GetWaterValuesOnIntegrationPoints(PointNumber, permeability_water, density_water);

        noalias( grad_water ) = GetGradientWater( DN_DX_PRESS );

        noalias( flow_water ) = GetFlowWater( PointNumber, RetentionParameters, grad_water, permeability_water, density_water );

        this->GetAirValuesOnIntegrationPoints(PointNumber, permeability_air);

        mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, DENSITY_AIR, density_air );

        mRetentionLawVector[PointNumber]->CalculateDerivative( RetentionParameters, DENSITY_AIR, bulk_air );

        noalias( grad_air ) = GetGradientAir( DN_DX_PRESS );

        noalias( flow_air ) = GetFlowAir( PointNumber, RetentionParameters, grad_air, permeability_air, density_air );

        density = GetAveragedDensity( saturation, porosity, density_soil, density_water, density_air );

        Drho_DdivU = Dn_DdivU * ( -density_soil + ( 1 - saturation ) * density_air + saturation * density_water );

        Drho_Dpw = porosity * ( density_air - density_water ) * DS_Dpc;

        Drho_Dpa = porosity * (( -density_air + density_water ) * DS_Dpc + ( 1.0 - saturation ) * bulk_air );

        noalias( DflowWater_Dpw ) = GetDerivativeDWaterFlowDpw( PointNumber, RetentionParameters, flow_water, DS_Dpc );

        DflowWater_Dgradpw = GetDerivativeDWaterFlowDGradpw( PointNumber, RetentionParameters, permeability_water, density_water );

        noalias( DflowWater_Dpa ) = GetDerivativeDWaterFlowDpa( PointNumber, RetentionParameters, flow_water, DS_Dpc );

        noalias( DflowAir_Dpw ) = GetDerivativeDAirFlowDpw( PointNumber, RetentionParameters, flow_air, DS_Dpc );

        DflowAir_Dgradpa = GetDerivativeDAirFlowDGradpa( PointNumber, RetentionParameters, permeability_air, density_air );

        noalias( DflowAir_Dpa ) = GetDerivativeDAirFlowDpa( PointNumber, RetentionParameters, flow_air, DS_Dpc, permeability_air, density_air, bulk_air );

        // KRATOS_WATCH(Dpc_Dt)
        // KRATOS_WATCH(flow_air)
        // KRATOS_WATCH(flow_water)
        // KRATOS_WATCH(density)
        // KRATOS_WATCH(density_water)
        // KRATOS_WATCH(density_air)
        // KRATOS_WATCH(Drho_DdivU)
        // KRATOS_WATCH(Drho_Dpw)
        // KRATOS_WATCH(Drho_Dpa)
        // KRATOS_WATCH(DflowWater_Dpw)
        // KRATOS_WATCH(DflowWater_Dpa)
        // KRATOS_WATCH(DflowWater_Dgradpw)
        // KRATOS_WATCH(DflowAir_Dpw)
        // KRATOS_WATCH(DflowAir_Dgradpa)
        // KRATOS_WATCH(DflowAir_Dpa)

        CalculateStressAndTangentialStiffness( PointNumber, StrainVector, StressVector, TanC_U,
                waterPressure, airPressure, saturation, DS_Dpc, TanC_W, TanC_A, rCurrentProcessInfo );

        if ( CalculateStiffnessMatrixFlag == true )
        {
            noalias( Help_K_UU ) = ZeroMatrix( MatSizeU, MatSizeU );
            CalculateStiffnessMatrixUU( Help_K_UU, TanC_U, B, DN_DX_DISP, N_DISP, Drho_DdivU );
            Help_K_UU *= (Weight * DetJ);

            noalias( Help_K_UW ) = ZeroMatrix( MatSizeU, MatSizeW );
            CalculateStiffnessMatrixUW( Help_K_UW, TanC_W, DN_DX_DISP, N_DISP, N_PRESS, Drho_Dpw );
            Help_K_UW *= (Weight * DetJ);

            noalias( Help_K_UA ) = ZeroMatrix( MatSizeU, MatSizeA );
            CalculateStiffnessMatrixUA( Help_K_UA, TanC_A, DN_DX_DISP, N_DISP, N_PRESS, Drho_Dpa );
            Help_K_UA *= (Weight * DetJ);

            noalias( Help_K_WU ) = ZeroMatrix( MatSizeW, MatSizeU );
            CalculateStiffnessMatrixWU( Help_K_WU, DN_DX_DISP, DN_DX_PRESS, N_PRESS, DS_Dpc, Dpc_Dt, Dn_DdivU );
            Help_K_WU *= Weight * DetJ; // * GetProperties()[SCALE]

            noalias( Help_K_WW ) = ZeroMatrix( MatSizeW, MatSizeW );
            CalculateStiffnessMatrixWW( Help_K_WW, DN_DX_DISP, DN_DX_PRESS, N_PRESS, porosity, DS_Dpc, Dpc_Dt, D2S_Dpc2, divU_Dt, DflowWater_Dgradpw, DflowWater_Dpw );
            Help_K_WW *= Weight * DetJ; // * GetProperties()[SCALE]

            noalias( Help_K_WA ) = ZeroMatrix( MatSizeW, MatSizeA );
            CalculateStiffnessMatrixWA( Help_K_WA, DN_DX_DISP, DN_DX_PRESS, N_PRESS, porosity, DS_Dpc, Dpc_Dt, D2S_Dpc2, divU_Dt, DflowWater_Dpa );
            Help_K_WA *= Weight * DetJ; // * GetProperties()[SCALE]

            noalias( Help_K_AU ) = ZeroMatrix( MatSizeA, MatSizeU );
            CalculateStiffnessMatrixAU( Help_K_AU, DN_DX_DISP, DN_DX_PRESS, N_PRESS, saturation, DS_Dpc, Dpc_Dt, Dn_DdivU, density_air, bulk_air, airPressure_Dt );
            Help_K_AU *= Weight * DetJ; // * GetProperties()[SCALE]

            noalias( Help_K_AW ) = ZeroMatrix( MatSizeA, MatSizeW );
            CalculateStiffnessMatrixAW( Help_K_AW, DN_DX_DISP, DN_DX_PRESS, N_PRESS, porosity, DS_Dpc, Dpc_Dt, D2S_Dpc2, divU_Dt, density_air, bulk_air, airPressure_Dt, DflowAir_Dpw, grad_air );
            Help_K_AW *= Weight * DetJ; // * GetProperties()[SCALE]

            noalias( Help_K_AA ) = ZeroMatrix( MatSizeA, MatSizeA );
            CalculateStiffnessMatrixAA( Help_K_AA, DN_DX_DISP, DN_DX_PRESS, N_PRESS, saturation, porosity, DS_Dpc, Dpc_Dt, D2S_Dpc2, divU_Dt, density_air, bulk_air, airPressure_Dt, flow_air, DflowAir_Dgradpa, DflowAir_Dpa, grad_air );
            Help_K_AA *= Weight * DetJ; // * GetProperties()[SCALE]

            // assemble to elemental stiffness matrix

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air,
                dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air,
                Help_K_UU );

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air,
                1, number_of_nodes_press, 1, 0,
                Help_K_UW );

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                    dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air,
                    1, number_of_nodes_air, 1, 1,
                    Help_K_UA );

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                    1, number_of_nodes_press, 1, 0,
                    dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air,
                    Help_K_WU );

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                    1, number_of_nodes_press, 1, 0,
                    1, number_of_nodes_press, 1, 0,
                    Help_K_WW );

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                    1, number_of_nodes_press, 1, 0,
                    1, number_of_nodes_air, 1, 1,
                    Help_K_WA );

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                    1, number_of_nodes_air, 1, 1,
                    dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air,
                    Help_K_AU );

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                    1, number_of_nodes_air, 1, 1,
                    1, number_of_nodes_press, 1, 0,
                    Help_K_AW );

            AssembleStiffnessFromSubMatrices( rLeftHandSideMatrix,
                    1, number_of_nodes_air, 1, 1,
                    1, number_of_nodes_air, 1, 1,
                    Help_K_AA );
        }

        if ( CalculateResidualVectorFlag == true )
        {
            //Calculation of spatial Loadvector
            noalias( Help_R_U ) = ZeroVector( MatSizeU );
            CalculateBodyForcesToRHSVectorU( Help_R_U, N_DISP, density );
            CalculateInternalForcesToRHSU( Help_R_U, B, StressVector );
            Help_R_U *= Weight * DetJ;

            noalias( Help_R_W ) = ZeroVector( MatSizeW );
            CalculateInternalForcesToRHSW( Help_R_W, DN_DX_DISP, DN_DX_PRESS, N_PRESS, saturation, porosity, DS_Dpc, Dpc_Dt, divU_Dt, flow_water );
            Help_R_W *= Weight * DetJ;// * GetProperties()[SCALE];

            noalias( Help_R_A ) = ZeroVector( MatSizeA );
            CalculateInternalForcesToRHSA( Help_R_A, DN_DX_DISP, DN_DX_PRESS, N_PRESS, saturation, porosity, density_air, bulk_air, DS_Dpc, Dpc_Dt, divU_Dt, flow_air, grad_air, airPressure_Dt );
            Help_R_A *= Weight * DetJ;// * GetProperties()[SCALE];

            // assemble to elemental internal force vector

            // if (Id() == 1 && PointNumber == 0)
            // {
            // //     KRATOS_WATCH(Weight)
            // //     KRATOS_WATCH(DetJ)
            //     KRATOS_WATCH(airPressure)
            //     KRATOS_WATCH(airPressure_Dt)
            //     KRATOS_WATCH(waterPressure)
            //     KRATOS_WATCH(waterPressure_Dt)
            //     KRATOS_WATCH(porosity)
            //     KRATOS_WATCH(capillaryPressure)
            //     KRATOS_WATCH(saturation)
            //     KRATOS_WATCH(porosity)
            //     KRATOS_WATCH(density_air)
            //     KRATOS_WATCH(permeability_air)
            //     KRATOS_WATCH(bulk_air)
            //     KRATOS_WATCH(DS_Dpc)
            //     KRATOS_WATCH(Dpc_Dt)
            //     KRATOS_WATCH(divU_Dt)
            //     KRATOS_WATCH(flow_air)
            //     KRATOS_WATCH(grad_air)
            //     KRATOS_WATCH(airPressure_Dt)
            //     // KRATOS_WATCH(StressVector)
            //     // KRATOS_WATCH(TanC_U)
            //     KRATOS_WATCH(norm_2(Help_R_U))
            //     KRATOS_WATCH(norm_2(Help_R_W))
            //     KRATOS_WATCH(norm_2(Help_R_A))
            //     KRATOS_WATCH(norm_frobenius(Help_K_UU))
            //     KRATOS_WATCH(norm_frobenius(Help_K_UW))
            //     KRATOS_WATCH(norm_frobenius(Help_K_UA))
            //     KRATOS_WATCH(norm_frobenius(Help_K_WU))
            //     KRATOS_WATCH(norm_frobenius(Help_K_WW))
            //     KRATOS_WATCH(norm_frobenius(Help_K_WA))
            //     KRATOS_WATCH(norm_frobenius(Help_K_AU))
            //     KRATOS_WATCH(norm_frobenius(Help_K_AW))
            //     KRATOS_WATCH(norm_frobenius(Help_K_AA))
            // }

            AssembleRHSFromSubVectors( rRightHandSideVector, dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air, Help_R_U );

            AssembleRHSFromSubVectors( rRightHandSideVector, 1, number_of_nodes_press, 1, 0, Help_R_W );

            AssembleRHSFromSubVectors( rRightHandSideVector, 1, number_of_nodes_air, 1, 1, Help_R_A );
        }
    }
    ///////////////////////////////////////////////////////////////////////
    // END Integration in space sum_(beta=0)^(number of quadrature points)
    ///////////////////////////////////////////////////////////////////////

    #ifdef DEBUG_AIR
    // if(Id() == 1)
    // {
    //     KRATOS_WATCH(Id())
    // //     KRATOS_WATCH(GetProperties().Id())
    // //     // KRATOS_WATCH(GetValue(ACTIVATION_LEVEL))
    // //     // KRATOS_WATCH(GetValue(IS_INACTIVE))
    // //     // KRATOS_WATCH(Is(ACTIVE))
    // //     KRATOS_WATCH(Help_R_U)
    // //     KRATOS_WATCH(Help_R_W)
    //     KRATOS_WATCH(Help_R_A)
    //     // KRATOS_WATCH(Help_K_UU)
    //     KRATOS_WATCH(Help_K_UW)
    //     KRATOS_WATCH(Help_K_UA)
    //     KRATOS_WATCH(Help_K_WU)
    //     KRATOS_WATCH(Help_K_WW)
    //     KRATOS_WATCH(Help_K_WA)
    //     KRATOS_WATCH(Help_K_AU)
    //     KRATOS_WATCH(Help_K_AW)
    //     KRATOS_WATCH(Help_K_AA)

    // //     // KRATOS_WATCH(rRightHandSideVector)
    // //     // KRATOS_WATCH(rLeftHandSideMatrix)

    // //     KRATOS_THROW_ERROR(std::runtime_error, "Stop here", "")

    //     EquationIdVectorType eq_ids;
    //     this->EquationIdVector( eq_ids, rCurrentProcessInfo );
    //     std::cout << "equation id:" << std::endl;
    //     for (unsigned int i = 0; i < eq_ids.size(); ++i)
    //         std::cout << " " << eq_ids[i];
    //     std::cout << std::endl;

    //     KRATOS_WATCH(norm_2(Help_R_U))
    //     KRATOS_WATCH(norm_2(Help_R_W))
    //     KRATOS_WATCH(norm_frobenius(Help_K_UU))
    //     KRATOS_WATCH(norm_frobenius(Help_K_UW))
    //     KRATOS_WATCH(norm_frobenius(Help_K_WU))
    //     KRATOS_WATCH(norm_frobenius(Help_K_WW))

    //     KRATOS_WATCH(norm_2(rRightHandSideVector))
    //     KRATOS_WATCH(norm_frobenius(rLeftHandSideMatrix))

    //     KRATOS_WATCH("--------------------")
    // }
    #endif

    #ifdef ENABLE_BEZIER_GEOMETRY
    // finalize the geometry
    GetGeometry().Clean();
    mpSubGeometry->Clean();
    #endif

    KRATOS_CATCH( "" )
}

void UPaPwSmallStrainElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    // TODO
    KRATOS_THROW_ERROR(std::logic_error, "TRANSIENT ANALYSIS is not implemented at UPaPwSmallStrainElement::", __FUNCTION__)
}

void UPaPwSmallStrainElement::CalculateDampingMatrix(MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int number_of_nodes_disp = GetGeometry().size();
    unsigned int number_of_nodes_press = mpSubGeometry->size();
    unsigned int number_of_nodes_air = mpSubGeometry->size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSizeU = number_of_nodes_disp * dim;
    unsigned int MatSizeW = number_of_nodes_press;
    unsigned int MatSizeA = number_of_nodes_air;
    unsigned int MatSize = MatSizeU + MatSizeW + MatSizeA;

    if ( rDampMatrix.size1() != MatSize || rDampMatrix.size2() != MatSize )
        rDampMatrix.resize( MatSize, MatSize, false );

    noalias( rDampMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS

    #ifdef ENABLE_BEZIER_GEOMETRY
    // initialize the geometry
    GetGeometry().Initialize(mThisIntegrationMethod);
    #endif

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    #ifdef ENABLE_BEZIER_GEOMETRY
    mpSubGeometry->Initialize(integration_points); // this is to make sure that the Bezier geometry of the pressure will have the same number of integration points as the geometry for displacement
    #endif

    const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =
        GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    const Matrix& Ncontainer_Pressure = mpSubGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

    double waterPressure;

    double airPressure;

    double saturation;

    double porosity;

    double DS_Dpc;

    double density_air;

    double bulk_air;

    Matrix Help_D_WU( MatSizeW, MatSizeU );

    Matrix Help_D_WW( MatSizeW, MatSizeW );

    Matrix Help_D_WA( MatSizeW, MatSizeA );

    Matrix Help_D_AU( MatSizeA, MatSizeU );

    Matrix Help_D_AW( MatSizeA, MatSizeW );

    Matrix Help_D_AA( MatSizeA, MatSizeA );

    Matrix DN_DX_DISP( number_of_nodes_disp, dim );

    Vector N_PRESS( number_of_nodes_press );

    Matrix InvJ(dim, dim);

    double Weight;

    double DetJ = 0.0;

    //initializing the Jacobian in the reference configuration
    GeometryType::JacobiansType J0;
    Matrix DeltaPosition(GetGeometry().size(), 3);

    for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
    {
        noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
    }

    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

    // /////////////////////////////////////////////////////////////////////////
    // //// Integration in space to compute the average of pressure shape function
    // /////////////////////////////////////////////////////////////////////////
    // Vector N_PRESS_averaged;
    // extern Variable<double> GEOMETRICAL_DOMAIN_SIZE;
    // double TotalDomainInitialSize = this->GetValue(GEOMETRICAL_DOMAIN_SIZE);
    // if(mIsStabilized)
    // {
    //     N_PRESS_averaged.resize( number_of_nodes_press, false );
    //     N_PRESS_averaged = ZeroVector( number_of_nodes_press );
    //     for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
    //     {
    //         Weight = integration_points[PointNumber].Weight();
    //         DetJ   = MathUtils<double>::Det(J0[PointNumber]);

    //         // Shape Functions on current spatial quadrature point
    //         noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

    //         noalias( N_PRESS_averaged ) += N_PRESS * Weight * DetJ;
    //     }
    //     N_PRESS_averaged /= TotalDomainInitialSize;
    // }

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        // Jacobian on current quadrature point
        MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ, DetJ );

        noalias( DN_DX_DISP ) = prod( DN_De_Displacement[PointNumber], InvJ );

        Weight = integration_points[PointNumber].Weight();

        // Shape Functions on current spatial quadrature point
        noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

        waterPressure = GetWaterPressure( N_PRESS );

        airPressure = GetAirPressure( N_PRESS );

        RetentionParameters.SetFluidPressure( waterPressure );

        RetentionParameters.SetAirPressure( airPressure );

        saturation = GetSaturation( PointNumber, RetentionParameters );

        porosity = GetPorosity( DN_DX_DISP );

        DS_Dpc = GetDerivativeDSaturationDpc( PointNumber, RetentionParameters );

        mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, DENSITY_AIR, density_air );

        mRetentionLawVector[PointNumber]->CalculateDerivative( RetentionParameters, DENSITY_AIR, bulk_air );

        // Calculation of spatial Stiffnes and Mass Matrix
        noalias( Help_D_WU ) = ZeroMatrix( MatSizeW, MatSizeU );
        CalculateDampingMatrixWU( Help_D_WU, DN_DX_DISP, N_PRESS, saturation );
        Help_D_WU *= Weight * DetJ; // * GetProperties()[SCALE]

        noalias( Help_D_WW ) = ZeroMatrix( MatSizeW, MatSizeW );
        CalculateDampingMatrixWW( Help_D_WW, DN_DX_DISP, N_PRESS, porosity, DS_Dpc );
        Help_D_WW *= Weight * DetJ; // * GetProperties()[SCALE]

        noalias( Help_D_WA ) = ZeroMatrix( MatSizeW, MatSizeA );
        CalculateDampingMatrixWA( Help_D_WA, DN_DX_DISP, N_PRESS, porosity, DS_Dpc );
        Help_D_WA *= Weight * DetJ; // * GetProperties()[SCALE]

        noalias( Help_D_AU ) = ZeroMatrix( MatSizeA, MatSizeU );
        CalculateDampingMatrixAU( Help_D_AU, DN_DX_DISP, N_PRESS, saturation );
        Help_D_AU *= Weight * DetJ; // * GetProperties()[SCALE]

        noalias( Help_D_AW ) = ZeroMatrix( MatSizeA, MatSizeW );
        CalculateDampingMatrixAW( Help_D_AW, DN_DX_DISP, N_PRESS, porosity, DS_Dpc );
        Help_D_AW *= Weight * DetJ; // * GetProperties()[SCALE]

        noalias( Help_D_AA ) = ZeroMatrix( MatSizeA, MatSizeA );
        CalculateDampingMatrixAA( Help_D_AA, DN_DX_DISP, N_PRESS, saturation, porosity, DS_Dpc, density_air, bulk_air );
        Help_D_AA *= Weight * DetJ; // * GetProperties()[SCALE]

        // if(mIsStabilized)
        //     CalculateDampingMatrixWWs( Help_D_WW, DN_DX_DISP, N_PRESS,
        //                               N_PRESS_averaged, Weight, DetJ );

        // assemble to elemental damping matrix

        AssembleStiffnessFromSubMatrices( rDampMatrix,
                1, number_of_nodes_press, 1, 0,
                dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air,
                Help_D_WU );

        AssembleStiffnessFromSubMatrices( rDampMatrix,
                1, number_of_nodes_press, 1, 0,
                1, number_of_nodes_press, 1, 0,
                Help_D_WW );

        AssembleStiffnessFromSubMatrices( rDampMatrix,
                1, number_of_nodes_press, 1, 0,
                1, number_of_nodes_air, 1, 1,
                Help_D_WA );

        AssembleStiffnessFromSubMatrices( rDampMatrix,
                1, number_of_nodes_air, 1, 1,
                dim, number_of_nodes_disp, 0, number_of_nodes_press + number_of_nodes_air,
                Help_D_AU );

        AssembleStiffnessFromSubMatrices( rDampMatrix,
                1, number_of_nodes_air, 1, 1,
                1, number_of_nodes_press, 1, 0,
                Help_D_AW );

        AssembleStiffnessFromSubMatrices( rDampMatrix,
                1, number_of_nodes_air, 1, 1,
                1, number_of_nodes_air, 1, 1,
                Help_D_AA );

        ///

        // if (Id() == 1 && PointNumber == 0)
        // {
        //     KRATOS_WATCH(airPressure)
        //     KRATOS_WATCH(waterPressure)
        //     KRATOS_WATCH(capillaryPressure)
        //     KRATOS_WATCH(saturation)
        //     KRATOS_WATCH(porosity)
        //     KRATOS_WATCH(density_air)
        //     KRATOS_WATCH(bulk_air)
        //     KRATOS_WATCH(norm_frobenius(Help_D_WU))
        //     KRATOS_WATCH(norm_frobenius(Help_D_WW))
        //     KRATOS_WATCH(norm_frobenius(Help_D_WA))
        //     KRATOS_WATCH(norm_frobenius(Help_D_AU))
        //     KRATOS_WATCH(norm_frobenius(Help_D_AW))
        //     KRATOS_WATCH(norm_frobenius(Help_D_AA))
        // }
    }
    ///////////////////////////////////////////////////////////////////////
    // END Integration in space sum_(beta=0)^(number of quadrature points)
    ///////////////////////////////////////////////////////////////////////

    // if(Id() == 1)
    // {
    //     std::cout << "Damping terms:" << std::endl;
    // //     KRATOS_WATCH(Id())
    //     KRATOS_WATCH(norm_frobenius(Help_D_WU))
    //     KRATOS_WATCH(norm_frobenius(Help_D_WW))
    // //     KRATOS_WATCH(rDampMatrix)

    //     KRATOS_WATCH(norm_frobenius(rDampMatrix))

    //     KRATOS_THROW_ERROR(std::runtime_error, "stop here", "")
    // }

    #ifdef ENABLE_BEZIER_GEOMETRY
    // finalize the geometry
    GetGeometry().Clean();
    mpSubGeometry->Clean();
    #endif

    KRATOS_CATCH( "" )
}

double UPaPwSmallStrainElement::GetAveragedDensity( const double& saturation,
    const double& porosity, const double& density_soil,
    const double& density_water, const double& density_air ) const
{
    double result;

    result = ( 1.0 - porosity ) * density_soil + porosity * ( saturation * density_water + ( 1.0 - saturation ) * density_air );

    return result;
}

Vector UPaPwSmallStrainElement::GetGradientWater( const Matrix& DN_DX_PRESS ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    Vector result( dim );
    noalias( result ) = ZeroVector( dim );

    for ( unsigned int i = 0 ; i < mpSubGeometry->size() ; i++ )
    {
        const double& presW = (*mpSubGeometry)[i].GetSolutionStepValue( WATER_PRESSURE );

        for ( unsigned int k = 0; k < dim; k++ )
        {
            result( k ) += presW * DN_DX_PRESS( i, k );
        }
    }

    return result;
}

Vector UPaPwSmallStrainElement::GetGradientAir( const Matrix& DN_DX_PRESS ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    Vector result( dim );
    noalias( result ) = ZeroVector( dim );

    for ( unsigned int i = 0 ; i < mpSubGeometry->size() ; i++ )
    {
        const double& presA_alpha = (*mpSubGeometry)[i].GetSolutionStepValue( AIR_PRESSURE );

        for ( unsigned int k = 0; k < dim; k++ )
        {
            result( k ) += presA_alpha * DN_DX_PRESS( i, k );
        }
    }

    return result;
}

double UPaPwSmallStrainElement::GetWaterPressure( const Vector& N_PRESS ) const
{
    double waterPressure = 0.0;

    for ( unsigned int i = 0 ; i < mpSubGeometry->size() ; i++ )
    {
        const double& presW_alpha = (*mpSubGeometry)[i].GetSolutionStepValue( WATER_PRESSURE );

        waterPressure += presW_alpha * N_PRESS( i );
    }

    return waterPressure;
}

double UPaPwSmallStrainElement::GetAirPressure( const Vector& N_PRESS ) const
{
    double airPressure = 0.0;

    for ( unsigned int i = 0 ; i < mpSubGeometry->size() ; i++ )
    {
        const double& presA_alpha = (*mpSubGeometry)[i].GetSolutionStepValue( AIR_PRESSURE );

        airPressure += presA_alpha * N_PRESS( i );
    }

    return airPressure;
}

double UPaPwSmallStrainElement::GetDerivativeDWaterPressureDt( const Vector& N_PRESS ) const
{
    double waterPressure_Dt = 0.0;

    for ( unsigned int i = 0 ; i < mpSubGeometry->size() ; i++ )
    {
        const double& presW_alpha_Dt = (*mpSubGeometry)[i].GetSolutionStepValue( WATER_PRESSURE_DT );

        waterPressure_Dt += presW_alpha_Dt * N_PRESS( i );
    }

    return waterPressure_Dt;
}

double UPaPwSmallStrainElement::GetDerivativeDAirPressureDt( const Vector& N_PRESS ) const
{
    double airPressure_Dt = 0.0;

    for ( unsigned int i = 0 ; i < mpSubGeometry->size() ; i++ )
    {
        const double& presA_alpha_Dt = (*mpSubGeometry)[i].GetSolutionStepValue( AIR_PRESSURE_DT );

        airPressure_Dt += presA_alpha_Dt * N_PRESS( i );
    }

    return airPressure_Dt;
}

double UPaPwSmallStrainElement::GetPorosity( const Matrix& DN_DX_DISP ) const
{
    double initialPorosity;
    bool fixPorosity;

    if ( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
    {
        initialPorosity = GetValue(POROSITY);
        fixPorosity = GetValue(FIX_POROSITY);
    }
    else
    {
        initialPorosity = GetProperties()[POROSITY];
        fixPorosity = GetProperties()[FIX_POROSITY];
    }

    double porosity;

    if (fixPorosity)
    {
        porosity = initialPorosity;
    }
    else
    {
        double divu = GetDivU(DN_DX_DISP);
        porosity = 1.0 - (1.0 - initialPorosity)*exp(-divu);
    }

    return porosity;
}

double UPaPwSmallStrainElement::GetDerivativeDPorosityDDivU( const Matrix& DN_DX_DISP ) const
{
    double initialPorosity;
    bool fixPorosity;

    if ( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
    {
        initialPorosity = GetValue(POROSITY);
        fixPorosity = GetValue(FIX_POROSITY);
    }
    else
    {
        initialPorosity = GetProperties()[POROSITY];
        fixPorosity = GetProperties()[FIX_POROSITY];
    }

    double porosity_divu;

    if (fixPorosity)
    {
        porosity_divu = 0.0;
    }
    else
    {
        double divu = GetDivU(DN_DX_DISP);
        porosity_divu = (1.0 - initialPorosity)*exp(-divu);
    }

    return porosity_divu;
}

//************************************************************************************
//************************************************************************************

double UPaPwSmallStrainElement::GetDivU( const Matrix& DN_DX_DISP ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double div = 0.0;

    Vector u_alpha( 3 );

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; i++ )
    {
        noalias( u_alpha ) =
        (
            GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT )
            - row(mInitialDisp, i)
        );

        for ( unsigned int k = 0; k < dim; k++ )
        {
            div += ( u_alpha( k ) )
                   * DN_DX_DISP( i, k );
        }
    }

    return div;
}

//************************************************************************************
//************************************************************************************

double UPaPwSmallStrainElement::GetDerivativeDDivUDt( const Matrix& DN_DX_DISP ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double div = 0.0;

    Vector u_alpha_Dt( dim );

    for ( unsigned int i = 0 ; i < GetGeometry().size() ; i++ )
    {
        noalias( u_alpha_Dt ) =
            GetGeometry()[i].GetSolutionStepValue( VELOCITY );

        for ( unsigned int k = 0; k < dim; k++ )
        {
            div += u_alpha_Dt( k ) * DN_DX_DISP( i, k );
        }
    }

    return div;
}

//************************************************************************************
//************************************************************************************
/// Soil-Water Characteristic curve
//************************************************************************************
//************************************************************************************
double UPaPwSmallStrainElement::GetSaturation( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters ) const
{
    double capillaryPressure;
    RetentionParameters.GetCapillaryPressure(capillaryPressure);
    if (capillaryPressure < 0.0) // no matric suction
        return 1.0;

    // // KRATOS_WATCH(capillaryPressure)
    double saturation;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, SATURATION, saturation );
    // // KRATOS_WATCH(saturation)
    return saturation;
}

double UPaPwSmallStrainElement::GetDerivativeDSaturationDpc( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters ) const
{
    double capillaryPressure;
    RetentionParameters.GetCapillaryPressure(capillaryPressure);
    if (capillaryPressure < 0.0) // no matric suction
        return 0.0;

    double DS_Dpc;
    mRetentionLawVector[PointNumber]->CalculateDerivative( RetentionParameters, SATURATION, DS_Dpc );
    return DS_Dpc;
}

double UPaPwSmallStrainElement::GetSecondDerivativeD2SaturationDpc2( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters ) const
{
    double capillaryPressure;
    RetentionParameters.GetCapillaryPressure(capillaryPressure);
    if (capillaryPressure < 0.0) // no matric suction
        return 0.0;

    double D2S_Dpc2;
    mRetentionLawVector[PointNumber]->CalculateSecondDerivative( RetentionParameters, SATURATION, D2S_Dpc2 );
    return D2S_Dpc2;
}

void UPaPwSmallStrainElement::GetWaterValuesOnIntegrationPoints( const std::size_t& ipoint,
    double& permeability_water, double& density_water ) const
{
    if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
    {
        permeability_water = GetValue(PERMEABILITY_WATER);
        density_water = GetValue(DENSITY_WATER);
    }
    else
    {
        permeability_water = GetProperties()[PERMEABILITY_WATER];
        density_water = GetProperties()[DENSITY_WATER];
    }
}

Vector UPaPwSmallStrainElement::GetFlowWater( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters,
        const Vector& grad_water, const double& permeability_water, const double& density_water ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double relPerm;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, PERMEABILITY_WATER, relPerm );

    const array_1d<double, 3>& gravity = GetProperties()[GRAVITY];

    // we get the G_CONSTANT instead taking the norm of the gravity because the G_CONSTANT is non-zero
    // while the gravity can be zero for some applications
    const double& g_constant = GetProperties()[G_CONSTANT];

    Vector result( dim );

    #ifdef CHECK_NAN
    if (density_water == 0.0)
        KRATOS_THROW_ERROR(std::logic_error, "DENSITY_WATER is zero at element", Id())
    #endif

    for ( unsigned int i = 0; i < dim; i++ )
    {
        result( i ) = -relPerm * permeability_water / ( density_water * g_constant ) * ( grad_water( i ) - density_water * gravity( i ) );
    }

    return result;
}

Vector UPaPwSmallStrainElement::GetDerivativeDWaterFlowDpw( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters,
        const Vector& flow_water, const double& DS_Dpc ) const
{
    double relPerm;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, PERMEABILITY_WATER, relPerm );

    double relPerm_pw;
    mRetentionLawVector[PointNumber]->CalculateDerivative( RetentionParameters, PERMEABILITY_WATER, relPerm_pw );
    relPerm_pw *= DS_Dpc * (-1.0);

    Vector result = relPerm_pw/relPerm * flow_water;

    return result;
}

Vector UPaPwSmallStrainElement::GetDerivativeDWaterFlowDpa( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters,
        const Vector& flow_water, const double& DS_Dpc ) const
{
    double relPerm;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, PERMEABILITY_WATER, relPerm );

    double relPerm_pa;
    mRetentionLawVector[PointNumber]->CalculateDerivative( RetentionParameters, PERMEABILITY_WATER, relPerm_pa );
    relPerm_pa *= DS_Dpc;

    Vector result = relPerm_pa/relPerm * flow_water;

    return result;
}

double UPaPwSmallStrainElement::GetDerivativeDWaterFlowDGradpw( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters,
        const double& permeability_water, const double& density_water ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double relPerm;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, PERMEABILITY_WATER, relPerm );

    const double& g_constant = GetProperties()[G_CONSTANT];

    double result;

    result = -relPerm * permeability_water / ( density_water * g_constant );

    return result;
}

void UPaPwSmallStrainElement::GetAirValuesOnIntegrationPoints( const std::size_t& ipoint,
    double& permeability_air ) const
{
    if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
    {
        permeability_air = GetValue(PERMEABILITY_AIR);
    }
    else
    {
        permeability_air = GetProperties()[PERMEABILITY_AIR];
    }
}

Vector UPaPwSmallStrainElement::GetFlowAir( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters,
        const Vector& grad_air, const double& permeability_air, const double& density_air ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double relPerm;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, PERMEABILITY_AIR, relPerm );

    const array_1d<double, 3>& gravity = GetProperties()[GRAVITY];

    const double& g_constant = GetProperties()[G_CONSTANT];

    Vector result( dim );

    noalias( result ) = ZeroVector( dim );

    for ( unsigned int i = 0; i < dim; i++ )
    {
        result( i ) = -relPerm * permeability_air / ( density_air * g_constant ) * ( grad_air( i ) - density_air * gravity( i ) );
    }

    return result;
}

Vector UPaPwSmallStrainElement::GetDerivativeDAirFlowDpa( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters,
        const Vector& flow_air, const double& DS_Dpc, const double& permeability_air, const double& density_air, const double& bulk_air ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    double relPerm;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, PERMEABILITY_AIR, relPerm );

    double relPerm_pa;
    mRetentionLawVector[PointNumber]->CalculateDerivative( RetentionParameters, PERMEABILITY_AIR, relPerm_pa );
    relPerm_pa *= (-DS_Dpc);

    const array_1d<double, 3>& gravity = GetProperties()[GRAVITY];

    const double& g_constant = GetProperties()[G_CONSTANT];

    Vector result( dim );

    for ( unsigned int i = 0; i < dim; i++ )
    {
        result( i ) = relPerm * permeability_air / ( density_air * g_constant ) * ( bulk_air * gravity[i] )
                      + relPerm_pa/relPerm * flow_air[i] - (bulk_air/density_air) * flow_air[i];
    }

    return result;
}

Vector UPaPwSmallStrainElement::GetDerivativeDAirFlowDpw( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters,
        const Vector& flow_air, const double& DS_Dpc ) const
{
    double relPerm;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, PERMEABILITY_AIR, relPerm );

    double relPerm_pw;
    mRetentionLawVector[PointNumber]->CalculateDerivative( RetentionParameters, PERMEABILITY_AIR, relPerm_pw );
    relPerm_pw *= DS_Dpc;

    Vector result = relPerm_pw/relPerm * flow_air;

    return result;
}

double UPaPwSmallStrainElement::GetDerivativeDAirFlowDGradpa( const std::size_t& PointNumber, RetentionLaw::Parameters& RetentionParameters,
        const double& permeability_air, const double& density_air ) const
{
    double relPerm;
    mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, PERMEABILITY_AIR, relPerm );

    const double& g_constant = GetProperties()[G_CONSTANT];

    double result;

    result = -relPerm * permeability_air / ( density_air * g_constant );

    return result;
}

void UPaPwSmallStrainElement::CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector ) const
{
    KRATOS_TRY

    unsigned int Dim = GetGeometry().WorkingSpaceDimension();
    unsigned int strain_size = Dim * (Dim+1) / 2;

    noalias( StrainVector ) = ZeroVector( strain_size );

    for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
    {
        for ( unsigned int item = 0; item < strain_size; item++ )
            for ( unsigned int dim = 0; dim < Dim; dim++ )
                StrainVector[item] += B( item, Dim * node + dim ) * ( Displacements( node, dim ) - mInitialDisp( node, dim ) );
    }

    KRATOS_CATCH( "" )
}

void UPaPwSmallStrainElement::CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX ) const
{
    KRATOS_TRY

    unsigned int Dim = GetGeometry().WorkingSpaceDimension();
    unsigned int strain_size = Dim * (Dim+1) / 2;
    unsigned int number_of_nodes_disp = GetGeometry().size();

    noalias( B_Operator ) = ZeroMatrix( strain_size, number_of_nodes_disp * Dim );

    if ( Dim == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes_disp; ++i )
        {
            B_Operator( 0, i*2 ) = DN_DX( i, 0 );
            B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*2 ) = DN_DX( i, 1 );
            B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
        }
    }
    else if ( Dim == 3 )
    {
        for ( unsigned int i = 0; i < number_of_nodes_disp; ++i )
        {
            B_Operator( 0, i*3 )     = DN_DX( i, 0 );
            B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
            B_Operator( 3, i*3 )     = DN_DX( i, 1 );
            B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
            B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
            B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
            B_Operator( 5, i*3 )     = DN_DX( i, 2 );
            B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
        }
    }

    KRATOS_CATCH( "" )
}

void UPaPwSmallStrainElement::CalculateEffectiveStress( Vector& StressVector, Matrix& tanC_W, Matrix& tanC_A,
        const double& waterPressure, const double& airPressure, const double& saturation, const double& DS_Dpc ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    noalias( tanC_W ) = ZeroMatrix( dim, dim );

    noalias( tanC_A ) = ZeroMatrix( dim, dim );

    for ( unsigned int i = 0; i < dim; ++i )
    {
        StressVector( i ) -= (( 1.0 - saturation ) * airPressure + saturation * waterPressure );
    }

    for ( unsigned int i = 0; i < dim; ++i )
    {
        tanC_W( i, i ) = ( DS_Dpc * ( waterPressure - airPressure ) - saturation );

        tanC_A( i, i ) = ( DS_Dpc * ( airPressure - waterPressure ) - ( 1.0 - saturation ) );
    }
}

void UPaPwSmallStrainElement::CalculateStressAndTangentialStiffness( const int& PointNumber,
        const Vector& StrainVector, Vector& StressVector, Matrix& tanC_U,
        const double& waterPressure, const double& airPressure,
        const double& saturation, const double& DS_Dpc,
        Matrix& tanC_W, Matrix& tanC_A,
        const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int strain_size = dim * (dim+1) / 2;

    if ( tanC_W.size1() != dim || tanC_W.size2() != dim )
        tanC_W.resize( dim, dim, false );

    noalias( tanC_W ) = ZeroMatrix( dim, dim );

    if ( tanC_A.size1() != dim || tanC_A.size2() != dim )
        tanC_A.resize( dim, dim, false );

    noalias( tanC_A ) = ZeroMatrix( dim, dim );

    if ( tanC_U.size1() != strain_size || tanC_U.size2() != strain_size )
        tanC_U.resize( strain_size, strain_size, false );

    noalias( tanC_U ) = ZeroMatrix( strain_size, strain_size );

    if ( StressVector.size() != strain_size )
        StressVector.resize( strain_size, false );

    //constitutive law
    ConstitutiveLaw::Parameters const_params;
    Vector StrainVectorCopy = StrainVector;
    const_params.SetStrainVector(StrainVectorCopy);
    const_params.SetStressVector(StressVector);
    const_params.SetConstitutiveMatrix(tanC_U);
    const_params.SetProcessInfo(rCurrentProcessInfo);
    const_params.SetMaterialProperties(GetProperties());
    const_params.SetElementGeometry(GetGeometry());

    //Set suction in const. law
    mConstitutiveLawVector[PointNumber]->SetValue( SUCTION, airPressure - waterPressure, rCurrentProcessInfo );

    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(const_params);

    CalculateEffectiveStress( StressVector, tanC_W, tanC_A, waterPressure, airPressure, saturation, DS_Dpc );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//CALCULATE EXTERNAL FORCEVECTORS DISPLACEMENT****************************************
//************************************************************************************
void UPaPwSmallStrainElement::CalculateBodyForcesToRHSVectorU( Vector& R, const Vector& N_DISP, const double& density ) const
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    const array_1d<double, 3>& gravity = GetProperties()[GRAVITY];

    #ifdef REPORT_ZERO_DENSITY
    if (this->Is(ACTIVE))
    {
        if (fabs(density) < 1.0e-6)
            std::cout << "Element " << Id() << " has zero DENSITY" << std::endl;
    }
    #endif

    for ( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
    {
        for ( unsigned int i = 0; i < dim; i++ )
        {
            R( prim*dim + i ) +=
                N_DISP( prim ) * density * gravity( i );
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//CALCULATE INTERNAL FORCEVECTORS DISPLACEMENT****************************************
//************************************************************************************

void UPaPwSmallStrainElement::CalculateInternalForcesToRHSU( Vector& R, const Matrix& B_Operator, const Vector& StressVector ) const
{
    KRATOS_TRY

    noalias( R ) -= prod( trans( B_Operator ), StressVector );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//CALCULATE STIFFNESS MATRICES DISPLACEMENT*******************************************
//************************************************************************************

void UPaPwSmallStrainElement::CalculateStiffnessMatrixUU( Matrix& Help_K_UU,
        const Matrix& tanC, const Matrix& B_Operator, const Matrix& DN_DX_DISP, const Vector& N_DISP,
        const double& Drho_DdivU ) const
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int number_of_nodes_disp = GetGeometry().size();

    const array_1d<double, 3>& gravity = GetProperties()[GRAVITY];

    for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
    {
        for ( unsigned int i = 0; i < dim; i++ )
        {
            for ( unsigned int sec = 0; sec < number_of_nodes_disp; sec++ )
            {
                for ( unsigned int j = 0; j < dim; j++ )
                {
                    Help_K_UU( prim*dim + i, sec*dim + j ) -= N_DISP( prim ) * Drho_DdivU * gravity( i ) * DN_DX_DISP( sec, j );
                }
            }
        }
    }

    noalias( Help_K_UU ) += prod( trans( B_Operator ), Matrix( prod( tanC, B_Operator ) ) );

    KRATOS_CATCH( "" )
}

void UPaPwSmallStrainElement::CalculateStiffnessMatrixUW( Matrix& Help_K_UW,
        const Matrix& tanC_W, const Matrix& DN_DX_DISP, const Vector& N_DISP, const Vector& N_PRESS,
        const double& Drho_Dpw ) const
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    unsigned int displacement_size = GetGeometry().size();

    const array_1d<double, 3>& gravity = GetProperties()[GRAVITY];

    for ( unsigned int prim = 0; prim < displacement_size; prim++ )
    {
        for ( unsigned int i = 0; i < dim; i++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_K_UW( prim*dim + i, sec ) -=
                    N_DISP( prim ) * Drho_Dpw * gravity( i ) * N_PRESS( sec ); //

                for ( unsigned int gamma = 0; gamma < dim; gamma++ )
                {
                    Help_K_UW( prim*dim + i, sec ) +=
                        ( DN_DX_DISP( prim, gamma ) * tanC_W( i, gamma ) * N_PRESS( sec ) ); //
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void UPaPwSmallStrainElement::CalculateStiffnessMatrixUA( Matrix& Help_K_UA,
        const Matrix& tanC_A, const Matrix& DN_DX_DISP, const Vector& N_DISP, const Vector& N_PRESS,
        const double& Drho_Dpa ) const
{
    KRATOS_TRY

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    unsigned int displacement_size = GetGeometry().size();

    const array_1d<double, 3>& gravity = GetProperties()[GRAVITY];

    for ( unsigned int prim = 0; prim < displacement_size; prim++ )
    {
        for ( unsigned int i = 0; i < dim; i++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_K_UA( prim*dim + i, sec ) -=
                    N_DISP( prim ) * Drho_Dpa * gravity( i ) * N_PRESS( sec ); //

                for ( unsigned int gamma = 0; gamma < dim; gamma++ )
                {
                    Help_K_UA( prim*dim + i, sec ) +=
                        ( DN_DX_DISP( prim, gamma ) * tanC_A( i, gamma ) * N_PRESS( sec ) ); //
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void UPaPwSmallStrainElement::CalculateInternalForcesToRHSW( Vector& Help_R_W,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
        const double& saturation, const double& porosity,
        const double& DS_Dpc, const double& Dpc_Dt, const double& divU_Dt,
        const Vector& flow_water ) const
{
    unsigned int pressure_size = mpSubGeometry->size();

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        Help_R_W( prim ) +=
            N_PRESS( prim ) * porosity * DS_Dpc * Dpc_Dt; //

        Help_R_W( prim ) +=
            N_PRESS( prim ) * saturation * divU_Dt; //

        Help_R_W( prim ) -=
            inner_prod( row( DN_DX_PRESS, prim ), flow_water ); //
    }
}

void UPaPwSmallStrainElement::CalculateStiffnessMatrixWU( Matrix& Help_K_WU,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
        const double& DS_Dpc, const double& Dpc_Dt, const double& Dn_DdivU ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    unsigned int displacement_size = GetGeometry().size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < displacement_size; sec++ )
        {
            for ( unsigned int j = 0; j < dim; j++ )
            {
                Help_K_WU( prim, sec*dim + j ) -=
                    N_PRESS( prim ) * Dn_DdivU * DS_Dpc * Dpc_Dt * DN_DX_DISP( sec, j ); //
            }
        }
    }
}

void UPaPwSmallStrainElement::CalculateStiffnessMatrixWW( Matrix& Help_K_WW,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
        const double& porosity, const double& DS_Dpc, const double& Dpc_Dt, const double& D2S_Dpc2,
        const double& divU_Dt, const double& DflowWater_Dgradpw, const Vector& DflowWater_Dpw ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_K_WW( prim, sec ) +=
                N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt * N_PRESS( sec ); //

            Help_K_WW( prim, sec ) +=
                N_PRESS( prim ) * DS_Dpc * divU_Dt * N_PRESS( sec ); //

            for ( unsigned int gamma = 0; gamma < dim; gamma++ )
            {
                Help_K_WW( prim, sec ) +=
                    DN_DX_PRESS( prim, gamma ) * DflowWater_Dpw( gamma )
                    * N_PRESS( sec ); //

                Help_K_WW( prim, sec ) +=
                    DN_DX_PRESS( prim, gamma ) * DflowWater_Dgradpw
                    * DN_DX_PRESS( sec, gamma ); //
            }
        }
    }
}

void UPaPwSmallStrainElement::CalculateStiffnessMatrixWA( Matrix& Help_K_WA,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
        const double& porosity, const double& DS_Dpc, const double& Dpc_Dt, const double& D2S_Dpc2,
        const double& divU_Dt, const Vector& DflowWater_Dpa ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_K_WA( prim, sec ) -=
                N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt * N_PRESS( sec ); //

            Help_K_WA( prim, sec ) -=
                N_PRESS( prim ) * DS_Dpc * divU_Dt * N_PRESS( sec ); //

            for ( unsigned int gamma = 0; gamma < dim; gamma++ )
            {
                Help_K_WA( prim, sec ) +=
                    DN_DX_PRESS( prim, gamma ) * DflowWater_Dpa( gamma ) * N_PRESS( sec ); //
            }
        }
    }
}

void UPaPwSmallStrainElement::CalculateDampingMatrixWU( Matrix& Help_D_WU,
        const Matrix& DN_DX_DISP, const Vector& N_PRESS,
        const double& saturation ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    unsigned int displacement_size = GetGeometry().size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < displacement_size; sec++ )
        {
            for ( unsigned int j = 0; j < dim; j++ )
            {
                Help_D_WU( prim, sec*dim + j ) -=
                    N_PRESS( prim ) * saturation * DN_DX_DISP( sec, j ); //
            }
        }
    }
}

void UPaPwSmallStrainElement::CalculateDampingMatrixWW( Matrix& Help_D_WW,
        const Matrix& DN_DX_DISP, const Vector& N_PRESS, const double& porosity,
        const double& DS_Dpc ) const
{
    unsigned int pressure_size = mpSubGeometry->size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_D_WW( prim, sec ) +=
                N_PRESS( prim ) * porosity * DS_Dpc * N_PRESS( sec ); //
        }
    }
}

void UPaPwSmallStrainElement::CalculateDampingMatrixWA( Matrix& Help_D_WA,
        const Matrix& DN_DX_DISP, const Vector& N_PRESS, const double& porosity,
        const double& DS_Dpc ) const
{
    unsigned int pressure_size = mpSubGeometry->size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_D_WA( prim, sec ) -=
                N_PRESS( prim ) * porosity * DS_Dpc * N_PRESS( sec ); //
        }
    }
}

void UPaPwSmallStrainElement::CalculateInternalForcesToRHSA( Vector& Help_R_A,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
        const double& saturation, const double& porosity,
        const double& density_air, const double& bulk_air,
        const double& DS_Dpc, const double& Dpc_Dt, const double& divU_Dt,
        const Vector& flow_air, const Vector& grad_air,
        const double& airPressure_Dt ) const
{
    unsigned int pressure_size = mpSubGeometry->size();

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        Help_R_A( prim ) -=
            N_PRESS( prim ) * porosity * DS_Dpc * Dpc_Dt; //

        Help_R_A( prim ) +=
            N_PRESS( prim ) * ( 1.0 - saturation ) * divU_Dt; //

        Help_R_A( prim ) +=
            N_PRESS( prim ) * porosity * ( 1.0 - saturation ) / density_air
            * bulk_air * airPressure_Dt; //

        Help_R_A( prim ) +=
            N_PRESS( prim ) * 1.0 / density_air
            * bulk_air * inner_prod( grad_air, flow_air ); //

        Help_R_A( prim ) -=
                inner_prod ( row( DN_DX_PRESS, prim ), flow_air); //
    }
}

void UPaPwSmallStrainElement::CalculateStiffnessMatrixAU( Matrix& Help_K_AU,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
        const double& saturation, const double& DS_Dpc, const double& Dpc_Dt,
        const double& Dn_DdivU, const double& density_air, const double& bulk_air,
        const double& airPressure_Dt ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    unsigned int displacement_size = GetGeometry().size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < displacement_size; sec++ )
        {
            for ( unsigned int j = 0; j < dim; j++ )
            {
                Help_K_AU( prim, sec*dim + j ) +=
                    N_PRESS( prim )
                    * Dn_DdivU * DS_Dpc * Dpc_Dt
                    * DN_DX_DISP( sec, j ); //

                Help_K_AU( prim, sec*dim + j ) -=
                    N_PRESS( prim )
                    * Dn_DdivU * ( 1.0 - saturation ) / density_air * bulk_air
                    * airPressure_Dt * DN_DX_DISP( sec, j ); //
            }
        }
    }
}

void UPaPwSmallStrainElement::CalculateStiffnessMatrixAW( Matrix& Help_K_AW,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
        const double& porosity, const double& DS_Dpc, const double& Dpc_Dt, const double& D2S_Dpc2,
        const double& divU_Dt, const double& density_air, const double& bulk_air,
        const double& airPressure_Dt, const Vector& DflowAir_Dpw, const Vector& grad_air ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_K_AW( prim, sec ) -=
                N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt
                * N_PRESS( sec ); //

            Help_K_AW( prim, sec ) -=
                N_PRESS( prim ) * DS_Dpc * divU_Dt
                * N_PRESS( sec ); //

            Help_K_AW( prim, sec ) -=
                N_PRESS( prim ) * DS_Dpc * porosity / density_air * bulk_air * airPressure_Dt
                * N_PRESS( sec ); //

            Help_K_AW( prim, sec ) -=
                    N_PRESS( prim ) * 1 / density_air * bulk_air
                    * inner_prod( grad_air, DflowAir_Dpw )
                    * N_PRESS( sec ); //

            Help_K_AW( prim, sec ) +=
                    inner_prod( row( DN_DX_PRESS, prim ), DflowAir_Dpw )
                    * N_PRESS( sec ); //
        }
    }
}

void UPaPwSmallStrainElement::CalculateStiffnessMatrixAA( Matrix& Help_K_AA,
        const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, const Vector& N_PRESS,
        const double& saturation,
        const double& porosity, const double& DS_Dpc, const double& Dpc_Dt, const double& D2S_Dpc2,
        const double& divU_Dt, const double& density_air, const double& bulk_air,
        const double& airPressure_Dt, const Vector& flow_air,
        const double& DflowAir_Dgradpa, const Vector& DflowAir_Dpa,
        const Vector& grad_air ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_K_AA( prim, sec ) +=
                N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt
                * N_PRESS( sec ); //

            Help_K_AA( prim, sec ) +=
                N_PRESS( prim ) * DS_Dpc * divU_Dt
                * N_PRESS( sec ); //

            Help_K_AA( prim, sec ) +=
                N_PRESS( prim ) * DS_Dpc * porosity / density_air
                * bulk_air * airPressure_Dt
                * N_PRESS( sec ); //

            Help_K_AA( prim, sec ) +=
                N_PRESS( prim )
                * porosity * ( 1.0 - saturation ) / pow( density_air, 2 )
                * pow( bulk_air, 2 ) * airPressure_Dt
                * N_PRESS( sec ); //

            for ( unsigned int gamma = 0; gamma < dim; gamma++ )
            {
                ///
                Help_K_AA( prim, sec ) +=
                    N_PRESS( prim ) *
                    1.0 / pow( density_air, 2 ) * ( pow( bulk_air, 2 ) )
                    * grad_air( gamma ) * flow_air( gamma )
                    * N_PRESS( sec ); //

                Help_K_AA( prim, sec ) -=
                    N_PRESS( prim ) *
                    1.0 / density_air * bulk_air * flow_air( gamma )
                    * DN_DX_PRESS( sec, gamma ); //

                Help_K_AA( prim, sec ) -=
                    N_PRESS( prim ) *
                    1.0 / density_air * bulk_air * grad_air( gamma ) * DflowAir_Dpa( gamma )
                    * N_PRESS( sec ); //

                Help_K_AA( prim, sec ) -=
                    N_PRESS( prim ) *
                    1.0 / density_air * bulk_air * grad_air( gamma ) * DflowAir_Dgradpa
                    * DN_DX_PRESS( sec, gamma ); //

                ///
                Help_K_AA( prim, sec ) +=
                    DN_DX_PRESS( prim, gamma ) * DflowAir_Dpa( gamma )
                    * N_PRESS( sec ); //

                Help_K_AA( prim, sec ) +=
                    DN_DX_PRESS( prim, gamma ) * DflowAir_Dgradpa
                    * DN_DX_PRESS( sec, gamma ); //
            }
        }
    }
}

void UPaPwSmallStrainElement::CalculateDampingMatrixAU( Matrix& Help_D_AU,
        const Matrix& DN_DX_DISP, const Vector& N_PRESS, const double& saturation ) const
{
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int pressure_size = mpSubGeometry->size();

    unsigned int displacement_size = GetGeometry().size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < displacement_size; sec++ )
        {
            for ( unsigned int j = 0; j < dim; j++ )
            {
                Help_D_AU( prim, sec*dim + j ) -=
                    N_PRESS( prim ) * ( 1 - saturation ) * DN_DX_DISP( sec, j ); //
            }
        }
    }
}

void UPaPwSmallStrainElement::CalculateDampingMatrixAW( Matrix& Help_D_AW,
        const Matrix& DN_DX_DISP, const Vector& N_PRESS,
        const double& porosity, const double& DS_Dpc ) const
{
    unsigned int pressure_size = mpSubGeometry->size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_D_AW( prim, sec ) -=
                N_PRESS( prim ) * porosity * DS_Dpc * N_PRESS( sec ); //
        }
    }
}

void UPaPwSmallStrainElement::CalculateDampingMatrixAA( Matrix& Help_D_AA,
        const Matrix& DN_DX_DISP, const Vector& N_PRESS,
        const double& saturation, const double& porosity, const double& DS_Dpc,
        const double& density_air, const double& bulk_air ) const
{
    unsigned int pressure_size = mpSubGeometry->size();

    for ( unsigned int prim = 0; prim < pressure_size; prim++ )
    {
        for ( unsigned int sec = 0; sec < pressure_size; sec++ )
        {
            Help_D_AA( prim, sec ) +=
                N_PRESS( prim ) * porosity * DS_Dpc * N_PRESS( sec ); //

            Help_D_AA( prim, sec ) -=
                N_PRESS( prim ) * porosity * ( 1 - saturation ) / density_air
                * bulk_air * N_PRESS( sec ); //
        }
    }
}

void UPaPwSmallStrainElement::SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable, const std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
    {
        std::stringstream ss;
        ss << "Error at UPaPwSmallStrainElement element " << Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl;
        ss << "rValues.size(): " << rValues.size() << std::endl;
        ss << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
    }
}

/**
 * Set a Vector Variable from outside
 * @param rVariable Global name of the variable to be calculated
 * @param rValues Vector of the values on the quadrature points
 * @param rCurrentProcessInfo
 */
void UPaPwSmallStrainElement::SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable, const std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
    {
        std::stringstream ss;
        ss << "Error at UPaPwSmallStrainElement element " << Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl;
        ss << "rValues.size(): " << rValues.size() << std::endl;
        ss << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }
}

/**
 * Set a Double Variable from outside
 * @param rVariable Global name of the variable to be calculated
 * @param rValue value on the quadrature points
 * @param rCurrentProcessInfo
 */
void UPaPwSmallStrainElement::SetValuesOnIntegrationPoints( const Variable<double>& rVariable,
        const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == K0 )
    {
        SetValue( K0, rValues[0] );
    }
    else
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            std::stringstream ss;
            ss << "Error at UPaPwSmallStrainElement element " << Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl;
            ss << "rValues.size(): " << rValues.size() << std::endl;
            ss << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
        }
    }
}
/**
 * Set an Int Variable from outside
 * @param rVariable Global name of the variable to be calculated
 * @param rValue value on the quadrature points
 * @param rCurrentProcessInfo
 */
void UPaPwSmallStrainElement::SetValuesOnIntegrationPoints( const Variable<int>& rVariable,
        const std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
    {
        std::stringstream ss;
        ss << "Error at UPaPwSmallStrainElement element " << Id() << ", The size of rValues and mConstitutiveLawVector is incompatible" << std::endl;
        ss << "rValues.size(): " << rValues.size() << std::endl;
        ss << "mConstitutiveLawVector.size(): " << mConstitutiveLawVector.size() << std::endl;
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
    {
        mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
    }
}

void UPaPwSmallStrainElement::SetValuesOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
        const std::vector< ConstitutiveLaw::Pointer >& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == CONSTITUTIVE_LAW )
    {
        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        if( mConstitutiveLawVector.size() != rValues.size() )
        {
            mConstitutiveLawVector.resize( rValues.size() );
        }

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            mConstitutiveLawVector[i] = rValues[i];
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Clean();
        #endif
    }
}

void UPaPwSmallStrainElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if (rVariable == DENSITY_AIR)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        if ( mConstitutiveLawVector.size() == 0 )
            return;

        unsigned int number_of_nodes_press = mpSubGeometry->size();

        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Initialize(mThisIntegrationMethod);
        mpSubGeometry->Initialize(mThisIntegrationMethod);
        #endif

        const Matrix& Ncontainer_Pressure = mpSubGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

        Vector N_PRESS( number_of_nodes_press );

        double capillaryPressure;

        double waterPressure;

        double airPressure;

        double saturation;

        // create general parametes of retention law
        RetentionLaw::Parameters RetentionParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space sum_(beta=0)^(number of quadrature points)
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
        {
            // Shape Functions on current spatial quadrature point
            if ( N_PRESS.size() != number_of_nodes_press )
                N_PRESS.resize( number_of_nodes_press );

            noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

            waterPressure = GetWaterPressure( N_PRESS );

            airPressure = GetAirPressure( N_PRESS );

            RetentionParameters.SetFluidPressure( waterPressure );

            RetentionParameters.SetAirPressure( airPressure );

            mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, DENSITY_AIR, rValues[PointNumber] );
        }

        return;
    }

    BaseType::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

void UPaPwSmallStrainElement::CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable,
        std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    //To Plot Fluid Flows
    if ( rVariable == WATER_FLOW || rVariable == AIR_FLOW )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int number_of_nodes_press = mpSubGeometry->size();

        Vector N_PRESS( number_of_nodes_press );

        Matrix DN_DX_PRESS( number_of_nodes_press, dim );

        Matrix DN_DX_DISP( dim, dim );

        double capillaryPressure;

        double waterPressure;

        double airPressure;

        double saturation;

        double permeability_water;

        double density_water;

        double permeability_air;

        double density_air;

        Vector grad_flow(dim);

        Vector fluid_flow( dim );

        Matrix InvJ(dim, dim);

        double DetJ;

        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Initialize(mThisIntegrationMethod);
        mpSubGeometry->Initialize(mThisIntegrationMethod);
        #endif

        const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =
                GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De_Pressure =
                mpSubGeometry->ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer_Pressure = mpSubGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

        //initializing the Jacobian in the reference configuration
        GeometryType::JacobiansType J0;

        Matrix DeltaPosition(GetGeometry().size(), dim);

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            noalias( row( DeltaPosition, node ) ) = GetGeometry()[node].Coordinates() - GetGeometry()[node].GetInitialPosition();
        }

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod, DeltaPosition );

        // create general parametes of retention law
        RetentionLaw::Parameters RetentionParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
        {
            // Shape Functions on current spatial quadrature point
            if ( N_PRESS.size() != number_of_nodes_press )
                N_PRESS.resize( number_of_nodes_press );

            noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ, DetJ );

            waterPressure = GetWaterPressure( N_PRESS );

            airPressure = GetAirPressure( N_PRESS );

            RetentionParameters.SetFluidPressure( waterPressure );

            RetentionParameters.SetAirPressure( airPressure );

            noalias( DN_DX_PRESS ) = prod( DN_De_Pressure[PointNumber], InvJ );

            noalias( DN_DX_DISP ) = prod( DN_De_Displacement[PointNumber], InvJ );

            if( rVariable == WATER_FLOW )
            {
                this->GetWaterValuesOnIntegrationPoints(PointNumber, permeability_water, density_water);

                noalias( grad_flow ) = GetGradientWater( DN_DX_PRESS );

                noalias( fluid_flow ) = GetFlowWater( PointNumber, RetentionParameters, grad_flow, permeability_water, density_water );
            }
            else if( rVariable == AIR_FLOW )
            {
                this->GetAirValuesOnIntegrationPoints(PointNumber, permeability_air);

                mRetentionLawVector[PointNumber]->CalculateValue( RetentionParameters, DENSITY_AIR, density_air );

                noalias( grad_flow ) = GetGradientAir( DN_DX_PRESS );

                noalias( fluid_flow ) = GetFlowAir( PointNumber, RetentionParameters, grad_flow, permeability_air, density_air );
            }

            rValues[PointNumber][0] = fluid_flow( 0 );

            rValues[PointNumber][1] = fluid_flow( 1 );

            if( dim == 3)
                rValues[PointNumber][2] = fluid_flow( 2 );
            else
                rValues[PointNumber][2] = 0.0;
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        GetGeometry().Clean();
        mpSubGeometry->Clean();
        #endif

        return;
    }

    BaseType::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

void UPaPwSmallStrainElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
    {
        mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
    }

    BaseType::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

void UPaPwSmallStrainElement::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != mConstitutiveLawVector.size() )
        rValues.resize( mConstitutiveLawVector.size() );

    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); ++PointNumber )
    {
        mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
    }

    BaseType::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

int UPaPwSmallStrainElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    BaseType::Check(rCurrentProcessInfo);

    if (this->GetProperties().Has(G_CONSTANT) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "G_CONSTANT is not provided for property ", this->GetProperties().Id());
    }

    KRATOS_CATCH( "" )

    return BaseType::Check(rCurrentProcessInfo);
}

void UPaPwSmallStrainElement::PrintData(std::ostream& rOStream) const
{
    BaseType::PrintData(rOStream);

    double density_soil, density_water, permeability_water, permeability_air;
    if ( GetValue(USE_DISTRIBUTED_PROPERTIES) )
    {
        density_soil = GetValue(DENSITY);
        density_water = GetValue(DENSITY_WATER);
        permeability_water = GetValue(PERMEABILITY_WATER);
        permeability_air = GetValue(PERMEABILITY_AIR);
    }
    else
    {
        density_soil = GetProperties()[DENSITY];
        density_water = GetProperties()[DENSITY_WATER];
        permeability_water = GetProperties()[PERMEABILITY_WATER];
        permeability_air = GetProperties()[PERMEABILITY_AIR];
    }
    rOStream << "DENSITY: " << density_soil << std::endl;
    rOStream << "DENSITY_WATER: " << density_water << std::endl;
    rOStream << "PERMEABILITY_WATER: " << permeability_water << std::endl;
    rOStream << "PERMEABILITY_AIR: " << permeability_air << std::endl;

    const array_1d<double, 3>& gravity = GetProperties()[GRAVITY];
    const double& g_constant = GetProperties()[G_CONSTANT];
    rOStream << "G_CONSTANT: " << g_constant << std::endl;
    rOStream << "GRAVITY: " << gravity << std::endl;
}

/***********************************************/

void UPaPwSmallStrainElement::AssembleRHSFromSubVectors( VectorType& rRightHandSideVector,
        const unsigned int& dim,
        const unsigned int& num,
        const unsigned int& stride,
        const unsigned int& begin,
        const Vector& R_P ) const
{
    for ( unsigned int prim = 0; prim < num; ++prim )
    {
        for ( unsigned int i = 0; i < dim; ++i )
        {
            rRightHandSideVector( begin + prim*(dim+stride) + i ) += R_P( prim*dim + i );
        }
    }
}

//************************************************************************************
//************************************************************************************

void UPaPwSmallStrainElement::AssembleStiffnessFromSubMatrices( MatrixType& rLeftHandSideMatrix,
            const unsigned int& dim1,
            const unsigned int& num1,
            const unsigned int& stride1,
            const unsigned int& begin1,
            const unsigned int& dim2,
            const unsigned int& num2,
            const unsigned int& stride2,
            const unsigned int& begin2,
            const Matrix& K_PQ) const
{
    for ( unsigned int prim = 0; prim < num1; ++prim )
    {
        for ( unsigned int i = 0; i < dim1; ++i )
        {
            for ( unsigned int sec = 0; sec < num2; ++sec )
            {
                for ( unsigned int j = 0; j < dim2; ++j )
                {
                    rLeftHandSideMatrix( begin1 + prim*(dim1+stride1) + i, begin2 + sec*(dim2+stride2) + j )
                            += K_PQ( prim*dim1 + i, sec*dim2 + j );
                }
            }
        }
    }
}

} // Namespace Kratos

#ifdef CHECK_NAN
#undef CHECK_NAN
#endif

#ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
#undef ENABLE_DEBUG_CONSTITUTIVE_LAW
#endif

#ifdef DEBUG_AIR
#undef DEBUG_AIR
#endif
