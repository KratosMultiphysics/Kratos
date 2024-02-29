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
#include "custom_elements/one-phase_flow/U_Pl_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPlElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    KRATOS_THROW_ERROR( std::logic_error, "calling the default Create method for a particular element ... illegal operation!!", "" )

    return Element::Pointer( new UPlElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPlElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_THROW_ERROR( std::logic_error, "calling the default Create method for a particular element ... illegal operation!!", "" )

    return Element::Pointer( new UPlElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
int UPlElement<TDim,TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo ) const
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
    if ( LIQUID_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "LIQUID_PRESSURE has Key zero at element", this->Id() )
    if ( DT_LIQUID_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DT_LIQUID_PRESSURE has Key zero at element", this->Id() )
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
        if ( Geom[i].SolutionStepsDataHas( LIQUID_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable LIQUID_PRESSURE on node ", Geom[i].Id() )
        if ( Geom[i].SolutionStepsDataHas( DT_LIQUID_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DT_LIQUID_PRESSURE on node ", Geom[i].Id() )
        if( Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false )
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VOLUME_ACCELERATION variable on node ", Geom[i].Id() );

        if ( Geom[i].HasDofFor( DISPLACEMENT_X ) == false || Geom[i].HasDofFor( DISPLACEMENT_Y ) == false || Geom[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", Geom[i].Id() )
        if ( Geom[i].HasDofFor( LIQUID_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable LIQUID_PRESSURE on node ", Geom[i].Id() )
    }

    // Verify ProcessInfo variables
    if ( VELOCITY_COEFFICIENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY_COEFFICIENT has Key zero at element", this->Id() )
    if ( DT_LIQUID_PRESSURE_COEFFICIENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DT_LIQUID_PRESSURE_COEFFICIENT has Key zero at element", this->Id() )
    if ( RAYLEIGH_ALPHA.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"RAYLEIGH_ALPHA has Key zero at element", this->Id() )
    if ( RAYLEIGH_BETA.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"RAYLEIGH_BETA has Key zero at element", this->Id() )

    // Verify properties
    if ( DENSITY_SOLID.Key() == 0 || Prop.Has( DENSITY_SOLID ) == false || Prop[DENSITY_SOLID] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_SOLID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( DENSITY_LIQUID.Key() == 0 || Prop.Has( DENSITY_LIQUID ) == false || Prop[DENSITY_LIQUID] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_LIQUID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( BULK_MODULUS_SOLID.Key() == 0 || Prop.Has( BULK_MODULUS_SOLID ) == false || Prop[BULK_MODULUS_SOLID] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"BULK_MODULUS_SOLID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( BULK_MODULUS_LIQUID.Key() == 0 || Prop.Has( BULK_MODULUS_LIQUID ) == false || Prop[BULK_MODULUS_LIQUID] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"BULK_MODULUS_LIQUID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( YOUNG_MODULUS.Key() == 0 || Prop.Has( YOUNG_MODULUS ) == false || Prop[YOUNG_MODULUS] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( DYNAMIC_VISCOSITY_LIQUID.Key() == 0 || Prop.Has( DYNAMIC_VISCOSITY_LIQUID ) == false || Prop[DYNAMIC_VISCOSITY_LIQUID] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DYNAMIC_VISCOSITY_LIQUID has Key zero, is not defined or has an invalid value at element", this->Id() )
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
void UPlElement<TDim,TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != NumGPoints )
        mConstitutiveLawVector.resize( NumGPoints );

    //Imposed Z strain vector initialisation
    if ( mImposedZStrainVector.size() != NumGPoints )
        mImposedZStrainVector.resize( NumGPoints );

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        mConstitutiveLawVector[i] = Prop[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->InitializeMaterial( Prop, Geom,row( Geom.ShapeFunctionsValues( mThisIntegrationMethod ), i ) );

        mImposedZStrainVector[i] = 0.0;
    }

    // Initializing the intrinsic permeability matrix from the properties
    PoroElementUtilities::CalculatePermeabilityMatrix(mIntrinsicPermeability,Prop,TDim);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = TNumNodes * (TDim + 1);
    unsigned int index = 0;

    if (rElementalDofList.size() != element_size)
      rElementalDofList.resize( element_size );

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        if constexpr (TDim>2)
            rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index++] = rGeom[i].pGetDof(LIQUID_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
GeometryData::IntegrationMethod UPlElement<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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
void UPlElement<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_THROW_ERROR(std::logic_error,"UPlElement::CalculateLeftHandSide not implemented","");

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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
void UPlElement<2,3>::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 3 * (2 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 3; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlElement<2,4>::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 4 * (2 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 4; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<  >
void UPlElement<3,4>::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 4 * (3 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 4; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<  >
void UPlElement<3,6>::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 6 * (3 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 6; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<  >
void UPlElement<3,8>::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = 8 * (3 + 1);
    unsigned int index = 0;

    if (rResult.size() != element_size)
      rResult.resize( element_size, false );

    for (unsigned int i = 0; i < 8; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
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
    const double Density = Porosity*Prop[DENSITY_LIQUID] + (1.0-Porosity)*Prop[DENSITY_SOLID];
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
void UPlElement<TDim,TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
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
void UPlElement<TDim,TNumNodes>::GetValuesVector( Vector& rValues, int Step ) const
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
void UPlElement<TDim,TNumNodes>::GetFirstDerivativesVector( Vector& rValues, int Step ) const
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
void UPlElement<TDim,TNumNodes>::GetSecondDerivativesVector( Vector& rValues, int Step ) const
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
void UPlElement<TDim,TNumNodes>::SetValuesOnIntegrationPoints( const Variable<double>& rVariable, const std::vector<double>& rValues,
                                                                const ProcessInfo& rCurrentProcessInfo )
{
    if (rVariable == IMPOSED_Z_STRAIN_VALUE) {
        for ( unsigned int i = 0; i < mImposedZStrainVector.size(); ++i ) {
            mImposedZStrainVector[i] = rValues[i];
        }

    } else {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable,const std::vector<Matrix>& rValues,
                                                            const ProcessInfo& rCurrentProcessInfo )
{
    if (rVariable == PERMEABILITY_MATRIX) {
        // Permeability is set only on the element, not on every GP
        noalias(mIntrinsicPermeability) = rValues[0];

    } else {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,
                                                                const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number, false );

    for ( unsigned int i = 0;  i < integration_points_number; i++ )
    {
        rValues[i] = 0.0;
        rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,std::vector<array_1d<double,3>>& rValues,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number );
    
    for ( unsigned int i = 0;  i < integration_points_number; i++ )
    {
        noalias(rValues[i]) = ZeroVector(3);
        rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,std::vector<Matrix>& rValues,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number );

    for ( unsigned int i = 0;  i < integration_points_number; i++ )
    {
        rValues[i].resize(TDim,TDim,false);
        noalias(rValues[i]) = ZeroMatrix(TDim,TDim);
        rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,std::vector<ConstitutiveLaw::Pointer>& rValues,
                                                                const ProcessInfo& rCurrentProcessInfo )
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

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateStiffnessMatrix method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateAll method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateRHS method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlElement<2,3>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ * this->GetProperties()[THICKNESS];
}

//----------------------------------------------------------------------------------------

template< >
void UPlElement<2,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ * this->GetProperties()[THICKNESS];
}

//----------------------------------------------------------------------------------------

template< >
void UPlElement<3,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------

template< >
void UPlElement<3,6>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------

template< >
void UPlElement<3,8>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateFluxResidual( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateFluxResidual method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateMixBodyForce( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateMixBodyForce method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateNegInternalForce( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateNegInternalForce method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateExplicitContributions (VectorType& rFluxResidual, VectorType& rBodyForce, VectorType& rNegInternalForces, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateExplicitContributions method for a particular element ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateLumpedMassMatrix( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType element_size = TNumNodes * (TDim + 1);

    //Resizing mass matrix
    if ( rLeftHandSideMatrix.size1() != element_size )
        rLeftHandSideMatrix.resize( element_size, element_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );

    const double& porosity = r_prop[POROSITY];
    const double density = porosity*r_prop[DENSITY_LIQUID] + (1.0-porosity)*r_prop[DENSITY_SOLID];

    const double thickness = (TDim == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    // LUMPED MASS MATRIX
    const double total_mass = r_geom.DomainSize() * density * thickness;

    Vector lumping_factors;
    lumping_factors = r_geom.LumpingFactors( lumping_factors );

    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        const double temp = lumping_factors[i] * total_mass;
        for ( IndexType j = 0; j < TDim; ++j ) {
            IndexType index = i * (TDim + 1) + j;
            rLeftHandSideMatrix(index,index) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateDampingMatrixWithLumpedMass( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Rayleigh Method (Damping Matrix = alpha*M + beta*K)

    const unsigned int element_size = TNumNodes * (TDim + 1);

    // Compute Mass Matrix
    MatrixType MassMatrix(element_size,element_size);

    this->CalculateLumpedMassMatrix(MassMatrix,rCurrentProcessInfo);

    // Compute Stiffness matrix
    MatrixType StiffnessMatrix(element_size,element_size);

    this->CalculateStiffnessMatrix(StiffnessMatrix,rCurrentProcessInfo);

    // Compute Damping Matrix
    if ( rLeftHandSideMatrix.size1() != element_size )
        rLeftHandSideMatrix.resize( element_size, element_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );

    noalias(rLeftHandSideMatrix) += rCurrentProcessInfo[RAYLEIGH_ALPHA] * MassMatrix;
    noalias(rLeftHandSideMatrix) += rCurrentProcessInfo[RAYLEIGH_BETA] * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateInertialForce( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

    // Compute Mass Matrix
    MatrixType MassMatrix(element_size,element_size);
    this->CalculateLumpedMassMatrix(MassMatrix,rCurrentProcessInfo);

    VectorType AccelerationVector(element_size);
    this->GetSecondDerivativesVector(AccelerationVector,0);

    noalias(rRightHandSideVector) = prod(MassMatrix,AccelerationVector);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::CalculateDampingForce( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

    // Compute Damping Matrix
    MatrixType DampingMatrix(element_size,element_size);
    this->CalculateDampingMatrixWithLumpedMass(DampingMatrix,rCurrentProcessInfo);

    VectorType VelocityVector(element_size);
    this->GetFirstDerivativesVector(VelocityVector,0);

    noalias(rRightHandSideVector) = prod(DampingMatrix,VelocityVector);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == NODAL_MASS ) {

        GeometryType& rGeom = GetGeometry();

        const unsigned int element_size = TNumNodes * (TDim + 1);

        Matrix lumped_mass_matrix = ZeroMatrix(element_size,element_size);
        this->CalculateLumpedMassMatrix(lumped_mass_matrix,rCurrentProcessInfo);

        for(SizeType i=0; i< TNumNodes; ++i) {

            SizeType index = (TDim + 1) * i;

            #pragma omp atomic
            rGeom[i].GetValue(NODAL_MASS) += lumped_mass_matrix(index,index);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlElement<TDim,TNumNodes>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double,3> >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();

    const unsigned int element_size = TNumNodes * (TDim + 1);

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL ) {
        // CD

        Vector flux_residual = ZeroVector(element_size);
        // this->CalculateFluxResidual(flux_residual,rCurrentProcessInfo);
        Vector body_force = ZeroVector(element_size);
        // this->CalculateMixBodyForce(body_force,rCurrentProcessInfo);
        Vector neg_internal_force = ZeroVector(element_size);
        // this->CalculateNegInternalForce(neg_internal_force,rCurrentProcessInfo);
        this->CalculateExplicitContributions(flux_residual,body_force,neg_internal_force,rCurrentProcessInfo);

        for(SizeType i=0; i< TNumNodes; ++i) {

            SizeType index = (TDim + 1) * i;

            array_1d<double, 3 >& r_external_force = rGeom[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
            array_1d<double, 3 >& r_internal_force = rGeom[i].FastGetSolutionStepValue(INTERNAL_FORCE);

            for(SizeType j=0; j<TDim; ++j) {
                #pragma omp atomic
                r_external_force[j] += body_force[index + j];

                #pragma omp atomic
                r_internal_force[j] += -neg_internal_force[index + j];
            }
        }
    } else if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == DAMPING_FORCE ) {
        // VV

        Vector flux_residual = ZeroVector(element_size);
        Vector body_force = ZeroVector(element_size);
        Vector neg_internal_force = ZeroVector(element_size);
        this->CalculateExplicitContributions(flux_residual,body_force,neg_internal_force,rCurrentProcessInfo);
        Vector damping_force = ZeroVector(element_size);
        this->CalculateDampingForce(damping_force,rCurrentProcessInfo);

        for(SizeType i=0; i< TNumNodes; ++i) {

            SizeType index = (TDim + 1) * i;

            array_1d<double, 3 >& r_external_force = rGeom[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
            array_1d<double, 3 >& r_internal_force = rGeom[i].FastGetSolutionStepValue(INTERNAL_FORCE);
            array_1d<double, 3 >& r_damping_force = rGeom[i].FastGetSolutionStepValue(DAMPING_FORCE);

            for(SizeType j=0; j<TDim; ++j) {
                #pragma omp atomic
                r_external_force[j] += body_force[index + j];

                #pragma omp atomic
                r_internal_force[j] += -neg_internal_force[index + j];

                #pragma omp atomic
                r_damping_force[j] += damping_force[index + j];
            }
        }
    } else if(rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == REACTION ) {
        // Residual/Reactions

        Vector flux_residual = ZeroVector(element_size);
        Vector body_force = ZeroVector(element_size);
        Vector neg_internal_force = ZeroVector(element_size);
        this->CalculateExplicitContributions(flux_residual,body_force,neg_internal_force,rCurrentProcessInfo);
        // Vector damping_force = ZeroVector(element_size);
        // this->CalculateDampingForce(damping_force,rCurrentProcessInfo);
        // Vector inertial_force = ZeroVector(element_size);
        // this->CalculateInertialForce(inertial_force,rCurrentProcessInfo);

        for(SizeType i=0; i< TNumNodes; ++i) {

            SizeType index = (TDim + 1) * i;

            array_1d<double, 3 >& r_force_residual = rGeom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            double& r_flux_residual = rGeom[i].FastGetSolutionStepValue(FLUX_RESIDUAL);

            for(SizeType j=0; j<TDim; ++j) {
                #pragma omp atomic
                r_force_residual[j] += body_force[index + j] + neg_internal_force[index + j]; // - damping_force[index + j] - inertial_force[index + j];
            }

            #pragma omp atomic
            r_flux_residual += flux_residual[index + TDim];
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlElement<2,3>;
template class UPlElement<2,4>;
template class UPlElement<3,4>;
template class UPlElement<3,6>;
template class UPlElement<3,8>;

} // Namespace Kratos
