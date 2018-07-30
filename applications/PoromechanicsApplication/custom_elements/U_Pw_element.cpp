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


// Application includes
#include "custom_elements/U_Pw_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    KRATOS_THROW_ERROR( std::logic_error, "calling the default Create method for a particular element ... illegal operation!!", "" )

    return Element::Pointer( new UPwElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_THROW_ERROR( std::logic_error, "calling the default Create method for a particular element ... illegal operation!!", "" )

    return Element::Pointer( new UPwElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
int UPwElement<TDim,TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();

    // verify nodal variables and dofs
    if ( DISPLACEMENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero at element", this->Id() )
    if ( VELOCITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero at element", this->Id() )
    if ( ACCELERATION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero at element", this->Id() )
    if ( WATER_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "WATER_PRESSURE has Key zero at element", this->Id() )
    if ( DT_WATER_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DT_WATER_PRESSURE has Key zero at element", this->Id() )
    if ( VOLUME_ACCELERATION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VOLUME_ACCELERATION has Key zero at element", this->Id() )

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        if ( Geom[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", Geom[i].Id() )
        if ( Geom[i].SolutionStepsDataHas( VELOCITY ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable VELOCITY on node ", Geom[i].Id() )
        if ( Geom[i].SolutionStepsDataHas( ACCELERATION ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable ACCELERATION on node ", Geom[i].Id() )
        if ( Geom[i].SolutionStepsDataHas( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable WATER_PRESSURE on node ", Geom[i].Id() )
        if ( Geom[i].SolutionStepsDataHas( DT_WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DT_WATER_PRESSURE on node ", Geom[i].Id() )
        if( Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false )
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VOLUME_ACCELERATION variable on node ", Geom[i].Id() );

        if ( Geom[i].HasDofFor( DISPLACEMENT_X ) == false || Geom[i].HasDofFor( DISPLACEMENT_Y ) == false || Geom[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", Geom[i].Id() )
        if ( Geom[i].HasDofFor( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable WATER_PRESSURE on node ", Geom[i].Id() )
    }

    // Verify ProcessInfo variables
    if ( VELOCITY_COEFFICIENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY_COEFFICIENT has Key zero at element", this->Id() )
    if ( DT_PRESSURE_COEFFICIENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DT_PRESSURE_COEFFICIENT has Key zero at element", this->Id() )
    if ( RAYLEIGH_ALPHA.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"RAYLEIGH_ALPHA has Key zero at element", this->Id() )
    if ( RAYLEIGH_BETA.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"RAYLEIGH_BETA has Key zero at element", this->Id() )

    // Verify properties
    if ( DENSITY_SOLID.Key() == 0 || Prop.Has( DENSITY_SOLID ) == false || Prop[DENSITY_SOLID] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_SOLID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( DENSITY_WATER.Key() == 0 || Prop.Has( DENSITY_WATER ) == false || Prop[DENSITY_WATER] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_WATER has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( BULK_MODULUS_SOLID.Key() == 0 || Prop.Has( BULK_MODULUS_SOLID ) == false || Prop[BULK_MODULUS_SOLID] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"BULK_MODULUS_SOLID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( BULK_MODULUS_FLUID.Key() == 0 || Prop.Has( BULK_MODULUS_FLUID ) == false || Prop[BULK_MODULUS_FLUID] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"BULK_MODULUS_FLUID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( YOUNG_MODULUS.Key() == 0 || Prop.Has( YOUNG_MODULUS ) == false || Prop[YOUNG_MODULUS] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( DYNAMIC_VISCOSITY.Key() == 0 || Prop.Has( DYNAMIC_VISCOSITY ) == false || Prop[DYNAMIC_VISCOSITY] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DYNAMIC_VISCOSITY has Key zero, is not defined or has an invalid value at element", this->Id() )
    const double& Porosity = Prop[POROSITY];
    if ( POROSITY.Key() == 0 || Prop.Has( POROSITY ) == false || Porosity < 0.0 || Porosity > 1.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"POROSITY has Key zero, is not defined or has an invalid value at element", this->Id() )
    const double& PoissonRatio = Prop[POISSON_RATIO];
    if ( POISSON_RATIO.Key() == 0 || Prop.Has( POISSON_RATIO ) == false || PoissonRatio < 0.0 || PoissonRatio >= 0.5 )
        KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO has Key zero, is not defined or has an invalid value at element", this->Id() )
    // Verify thickness in a 2D case
    if ( TDim == 2 )
        if ( THICKNESS.Key() == 0 || Prop.Has( THICKNESS ) == false || Prop[THICKNESS] <= 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero, is not defined or has an invalid value at element", this->Id() )

    // If this is a 2D problem, check that nodes are in XY plane
    if ( TDim == 2 )
    {
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            if(Geom[i].Z() != 0.0)
                KRATOS_THROW_ERROR( std::logic_error," Node with non-zero Z coordinate found. Id: ", Geom[i].Id() )
        }
    }

    return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::Initialize()
{
    KRATOS_TRY

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != NumGPoints )
        mConstitutiveLawVector.resize( NumGPoints );

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        mConstitutiveLawVector[i] = Prop[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->InitializeMaterial( Prop, Geom,row( Geom.ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = TNumNodes * (TDim + 1);
    unsigned int index = 0;

    if (rElementalDofList.size() != element_size)
      rElementalDofList.resize( element_size );

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        if(TDim>2)
            rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index++] = rGeom[i].pGetDof(WATER_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != element_size )
        rLeftHandSideMatrix.resize( element_size, element_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_THROW_ERROR(std::logic_error,"UPwElement::CalculateLeftHandSide not implemented","");

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwElement<2,3>::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 3 * (2 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 3; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwElement<2,4>::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 4 * (2 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 4; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<  >
void UPwElement<3,4>::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 4 * (3 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 4; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<  >
void UPwElement<3,6>::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 6 * (3 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 6; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<  >
void UPwElement<3,8>::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 8 * (3 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 8; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resizing mass matrix
    if ( rMassMatrix.size1() != element_size )
        rMassMatrix.resize( element_size, element_size, false );
    noalias( rMassMatrix ) = ZeroMatrix( element_size, element_size );

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    //Defining shape functions and the determinant of the jacobian at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer,mThisIntegrationMethod);

    //Defining necessary variables
    double IntegrationCoefficient;
    const double& Porosity = Prop[POROSITY];
    const double Density = Porosity*Prop[DENSITY_WATER] + (1.0-Porosity)*Prop[DENSITY_SOLID];
    BoundedMatrix<double,TDim+1, TNumNodes*(TDim+1)> Nut = ZeroMatrix(TDim+1, TNumNodes*(TDim+1));

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
    {
        PoroElementUtilities::CalculateNuElementMatrix(Nut,NContainer,GPoint);

        //calculating weighting coefficient for integration
        this->CalculateIntegrationCoefficient( IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Adding contribution to Mass matrix
        noalias(rMassMatrix) += Density*prod(trans(Nut),Nut)*IntegrationCoefficient;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Rayleigh Method (Damping Matrix = alpha*M + beta*K)

    const unsigned int element_size = TNumNodes * (TDim + 1);

    // Compute Mass Matrix
    MatrixType MassMatrix(element_size,element_size);

    this->CalculateMassMatrix(MassMatrix,rCurrentProcessInfo);

    // Compute Stiffness matrix
    MatrixType StiffnessMatrix(element_size,element_size);

    this->CalculateStiffnessMatrix(StiffnessMatrix,rCurrentProcessInfo);

    // Compute Damping Matrix
    if ( rDampingMatrix.size1() != element_size )
        rDampingMatrix.resize( element_size, element_size, false );
    noalias( rDampingMatrix ) = ZeroMatrix( element_size, element_size );

    noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_ALPHA] * MassMatrix;
    noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_BETA] * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::GetValuesVector( Vector& rValues, int Step )
{
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int element_size = TNumNodes * (TDim + 1);
    unsigned int index = 0;

    if ( rValues.size() != element_size )
        rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );
        if ( TDim > 2 )
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
        rValues[index++] = 0.0;
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int element_size = TNumNodes * (TDim + 1);
    unsigned int index = 0;

    if ( rValues.size() != element_size )
        rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
        if ( TDim > 2 )
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
        rValues[index++] = 0.0;
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int element_size = TNumNodes * (TDim + 1);
    unsigned int index = 0;

    if ( rValues.size() != element_size )
        rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
        if ( TDim > 2 )
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
        rValues[index++] = 0.0;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::SetValueOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,
                                                                const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::GetValueOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,
                                                                const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == VON_MISES_STRESS )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else //if(rVariable == DAMAGE_VARIABLE || rVariable == STATE_VARIABLE)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); i++ )
        {
            rValues[i] = 0.0;
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
        }
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,std::vector<array_1d<double,3>>& rValues,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == FLUID_FLUX_VECTOR)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); i++ )
        {
            noalias(rValues[i]) = ZeroVector(3);
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
        }
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,std::vector<Matrix>& rValues,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    if( rVariable == CAUCHY_STRESS_TENSOR || rVariable == TOTAL_STRESS_TENSOR || rVariable == GREEN_LAGRANGE_STRAIN_TENSOR || rVariable == PERMEABILITY_MATRIX )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); i++ )
        {
            rValues[i].resize(TDim,TDim,false);
            noalias(rValues[i]) = ZeroMatrix(TDim,TDim);
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
        }
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,std::vector<ConstitutiveLaw::Pointer>& rValues,
                                                                const ProcessInfo& rCurrentProcessInfo )
{
    if(rVariable == CONSTITUTIVE_LAW)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for(unsigned int i=0; i < mConstitutiveLawVector.size(); i++)
            rValues[i] = mConstitutiveLawVector[i];
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateStiffnessMatrix method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateAll method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateRHS method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwElement<2,3>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ * this->GetProperties()[THICKNESS];
}

//----------------------------------------------------------------------------------------

template< >
void UPwElement<2,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ * this->GetProperties()[THICKNESS];
}

//----------------------------------------------------------------------------------------

template< >
void UPwElement<3,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------

template< >
void UPwElement<3,6>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------

template< >
void UPwElement<3,8>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwElement<2,3>;
template class UPwElement<2,4>;
template class UPwElement<3,4>;
template class UPwElement<3,6>;
template class UPwElement<3,8>;

} // Namespace Kratos
