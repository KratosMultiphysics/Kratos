// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// Application includes
#include "custom_elements/Pw_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer PwElement<TDim,TNumNodes>::Create(IndexType NewId,
                                                               NodesArrayType const& ThisNodes,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new PwElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer PwElement<TDim,TNumNodes>::Create(IndexType NewId,
                                                               GeometryType::Pointer pGeom,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new PwElement( NewId, pGeom, pProperties ) );
}
//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-PwElement::Initialize()") << this->Id() << std::endl;

    const PropertiesType &Prop = this->GetProperties();
    const GeometryType &Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    // pointer to constitutive laws
    if ( mConstitutiveLawVector.size() != NumGPoints )
        mConstitutiveLawVector.resize( NumGPoints );

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        mConstitutiveLawVector[i] = nullptr;
    }

    if ( mRetentionLawVector.size() != NumGPoints )
        mRetentionLawVector.resize( NumGPoints );
    for ( unsigned int i = 0; i < mRetentionLawVector.size(); i++ )
    {
        //RetentionLawFactory::Pointer pRetentionFactory;
        mRetentionLawVector[i] = RetentionLawFactory::Clone(Prop);
        mRetentionLawVector[i]->
            InitializeMaterial( Prop,
                                Geom,
                                row( Geom.ShapeFunctionsValues( this->GetIntegrationMethod() ), i ) );
    }

    mIsInitialised = true;

    KRATOS_CATCH( "" )

    // KRATOS_INFO("1-PwElement::Initialize()") << std::endl;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int PwElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-PwElement::Check()") << this->Id() << std::endl;

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();

    if (Geom.DomainSize() < 1.0e-15)
        KRATOS_THROW_ERROR( std::logic_error, "DomainSize < 1.0e-15 for the element ", this->Id() )

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        if ( Geom[i].SolutionStepsDataHas( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable WATER_PRESSURE on node ", Geom[i].Id() )
        if ( Geom[i].SolutionStepsDataHas( DT_WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DT_WATER_PRESSURE on node ", Geom[i].Id() )
        if ( Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false )
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VOLUME_ACCELERATION variable on node ", Geom[i].Id() );
        if ( Geom[i].HasDofFor( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable WATER_PRESSURE on node ", Geom[i].Id() )
    }

    // Verify ProcessInfo variables

    // Verify properties
    if ( Prop.Has( DENSITY_SOLID ) == false || Prop[DENSITY_SOLID] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_SOLID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( Prop.Has( DENSITY_WATER ) == false || Prop[DENSITY_WATER] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_WATER has Key zero, is not defined or has an invalid value at element", this->Id() )


    if ( Prop.Has( BULK_MODULUS_SOLID ) == false || Prop[BULK_MODULUS_SOLID] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"BULK_MODULUS_SOLID has Key zero, is not defined or has an invalid value at element", this->Id() )

    if ( Prop.Has( POROSITY ) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"POROSITY has Key zero, is not defined or has an invalid value at element", this->Id() )


    if ( TDim == 2 )
    {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            if (Geom[i].Z() != 0.0)
                KRATOS_THROW_ERROR( std::logic_error," Node with non-zero Z coordinate found. Id: ", Geom[i].Id() )
        }
    }

    // Verify specific properties
    if ( Prop.Has( BULK_MODULUS_FLUID ) == false || Prop[BULK_MODULUS_FLUID] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "BULK_MODULUS_FLUID has Key zero, is not defined or has an invalid value at element",
                            this->Id() )

    if ( Prop.Has( DYNAMIC_VISCOSITY ) == false || Prop[DYNAMIC_VISCOSITY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "DYNAMIC_VISCOSITY has Key zero, is not defined or has an invalid value at element",
                            this->Id() )

    if ( Prop.Has( PERMEABILITY_XX ) == false || Prop[PERMEABILITY_XX] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "PERMEABILITY_XX has Key zero, is not defined or has an invalid value at element",
                            this->Id() )

    if ( Prop.Has( PERMEABILITY_YY ) == false || Prop[PERMEABILITY_YY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "PERMEABILITY_YY has Key zero, is not defined or has an invalid value at element",
                            this->Id() )

    if ( Prop.Has( PERMEABILITY_XY ) == false || Prop[PERMEABILITY_XY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "PERMEABILITY_XY has Key zero, is not defined or has an invalid value at element",
                            this->Id() )
    if (TDim > 2)
    {
        if ( Prop.Has( PERMEABILITY_ZZ ) == false || Prop[PERMEABILITY_ZZ] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_ZZ has Key zero, is not defined or has an invalid value at element", this->Id() )

        if ( Prop.Has( PERMEABILITY_YZ ) == false || Prop[PERMEABILITY_YZ] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_YZ has Key zero, is not defined or has an invalid value at element", this->Id() )

        if ( Prop.Has( PERMEABILITY_ZX ) == false || Prop[PERMEABILITY_ZX] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_ZX has Key zero, is not defined or has an invalid value at element", this->Id() )
    }


    // Verify that the constitutive law has the correct dimension

    // Check constitutive law
    if ( mRetentionLawVector.size() > 0 ) {
        return mRetentionLawVector[0]->Check( Prop, rCurrentProcessInfo );
    }

    // KRATOS_INFO("1-PwElement::Check()") << std::endl;

    return 0;

    KRATOS_CATCH( "" );
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-PwElement::InitializeSolutionStep()") << this->Id() << std::endl;

    if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);

    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(Geom, this->GetProperties(), rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        // Initialize retention law
        mRetentionLawVector[GPoint]->InitializeSolutionStep(RetentionParameters);
    }

    // KRATOS_INFO("1-PwElement::InitializeSolutionStep()") << std::endl;
    KRATOS_CATCH("");
}


//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // nothing

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // nothing

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-PwElement::FinalizeSolutionStep()") << std::endl;

    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(Geom, this->GetProperties(), rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
    {
        // retention law
        mRetentionLawVector[GPoint]->FinalizeSolutionStep(RetentionParameters);
    }

    // KRATOS_INFO("1-PwElement::FinalizeSolutionStep()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<double>& rVariable,
                                  std::vector<double>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-PwElement::CalculateOnIntegrationPoints()") << std::endl;

    if (rVariable == DEGREE_OF_SATURATION ||
        rVariable == EFFECTIVE_SATURATION ||
        rVariable == BISHOP_COEFICIENT ||
        rVariable == DERIVATIVE_OF_SATURATION ||
        rVariable == RELATIVE_PERMEABILITY )
    {
        UPwSmallStrainElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
    else
    {
        if ( rOutput.size() != mRetentionLawVector.size() )
            rOutput.resize(mRetentionLawVector.size());

        for ( unsigned int i = 0;  i < mRetentionLawVector.size(); i++ )
        {
            rOutput[i] = mRetentionLawVector[i]->GetValue( rVariable, rOutput[i] );
        }
    }
    
    // KRATOS_INFO("1-PwElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable,
                                  std::vector<array_1d<double,3>>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-PwElement::CalculateOnIntegrationPoints<double,3>()") << std::endl;

    if (rVariable == FLUID_FLUX_VECTOR)
    {
        UPwSmallStrainElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
    else
    {
        if ( rOutput.size() != mRetentionLawVector.size() )
            rOutput.resize(mRetentionLawVector.size());

        for ( unsigned int i = 0;  i < mRetentionLawVector.size(); i++ )
        {
            noalias(rOutput[i]) = ZeroVector(3);
            rOutput[i] = mRetentionLawVector[i]->GetValue( rVariable, rOutput[i] );
        }
    }

    // KRATOS_INFO("1-PwElement::CalculateOnIntegrationPoints<double,3>()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                 std::vector<Matrix>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-PwElement::CalculateOnIntegrationPoints()") << rVariable << std::endl;


    if (rVariable == PERMEABILITY_MATRIX)
    {
        UPwSmallStrainElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
    else
    {
        if ( rOutput.size() != mConstitutiveLawVector.size() )
            rOutput.resize(mConstitutiveLawVector.size());

        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); i++ )
        {
            rOutput[i].resize(TDim,TDim,false);
            noalias(rOutput[i]) = ZeroMatrix(TDim,TDim);
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
        }
    }

    // KRATOS_INFO("1-PwElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& rCurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    // KRATOS_INFO("0-PwElement::CalculateAll()") << std::endl;

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( this->GetIntegrationMethod() );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      rCurrentProcessInfo );

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(Geom, this->GetProperties(), rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        //Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>( Variables.BodyAcceleration,
                                                                Variables.NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = CalculateBulkModulus(Variables.ConstitutiveMatrix);
        this->InitializeBiotCoefficients(Variables, BulkModulus);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient,
                                              Variables.detJ0,
                                              IntegrationPoints[GPoint].Weight());

        //Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if (CalculateResidualVectorFlag)  this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    // KRATOS_INFO("1-PwElement::CalculateAll()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    InitializeElementVariables( ElementVariables& rVariables,
                                const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-PwElement::InitializeElementVariables()") << std::endl;

    //Properties variables
    this->InitializeProperties( rVariables );

    //ProcessInfo variables
    rVariables.VelocityCoefficient = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    //Nodal Variables
    this->InitializeNodalPorePressureVariables( rVariables );

    //Variables computed at each GP
    noalias(rVariables.Nu) = ZeroMatrix(TDim, TNumNodes*TDim);
    rVariables.Np.resize(TNumNodes,false);
    rVariables.GradNpT.resize(TNumNodes,TDim,false);

    //General Variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    // shape functions
    (rVariables.NContainer).resize(NumGPoints, TNumNodes, false);
    rVariables.NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );

    // gradient of shape functions and determinant of Jacobian
    (rVariables.detJContainer).resize(NumGPoints,false);

    Geom.ShapeFunctionsIntegrationPointsGradients( rVariables.DN_DXContainer,
                                                   rVariables.detJContainer,
                                                   this->GetIntegrationMethod() );

    // Retention law
    rVariables.FluidPressure = 0.0;
    rVariables.DegreeOfSaturation = 1.0;
    rVariables.DerivativeOfSaturation = 0.0;
    rVariables.RelativePermeability = 1.0;
    rVariables.BishopCoefficient = 1.0;

    // KRATOS_INFO("1-PwElement::InitializeElementVariables()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-PwElement::CalculateAndAddLHS()") << std::endl;
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);

    // KRATOS_INFO("1-PwElement::CalculateAndAddLHS()") << std::endl;

    KRATOS_CATCH("");
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-PwElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    noalias(rVariables.PMatrix) = - PORE_PRESSURE_SIGN_FACTOR 
                                  * rVariables.DtPressureCoefficient
                                  * rVariables.BiotModulusInverse
                                  * outer_prod(rVariables.Np, rVariables.Np)
                                  * rVariables.IntegrationCoefficient;

    //Distribute compressibility block matrix into the elemental matrix
    GeoElementUtilities::
        AssemblePBlockMatrix< TDim, TNumNodes >(rLeftHandSideMatrix,
                                                rVariables.PMatrix);

    // KRATOS_INFO("1-PwElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-PwElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    noalias(rVariables.PDimMatrix) = - PORE_PRESSURE_SIGN_FACTOR 
                                     * prod(rVariables.GradNpT, rVariables.PermeabilityMatrix);

    noalias(rVariables.PMatrix) =  rVariables.DynamicViscosityInverse
                                 * rVariables.RelativePermeability
                                 * prod(rVariables.PDimMatrix, trans(rVariables.GradNpT))
                                 * rVariables.IntegrationCoefficient;

    //Distribute permeability block matrix into the elemental matrix
    GeoElementUtilities::AssemblePBlockMatrix< TDim, TNumNodes >(rLeftHandSideMatrix, rVariables.PMatrix);

    // KRATOS_INFO("1-PwElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-PwElement::CalculateAndAddRHS()") << std::endl;


    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    // VG: TODO: Check
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    // KRATOS_INFO("1-PwElement::CalculateAndAddRHS()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityFlow( VectorType& rRightHandSideVector,
                                     ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-PwElement::CalculateAndAddPermeabilityFlow()") << std::endl;

    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT, rVariables.PermeabilityMatrix);

    noalias(rVariables.PMatrix) = - PORE_PRESSURE_SIGN_FACTOR
                                  * rVariables.DynamicViscosityInverse
                                  * rVariables.RelativePermeability
                                  * prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))
                                  * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = - prod(rVariables.PMatrix, rVariables.PressureVector);

    //Distribute permeability block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<TDim, TNumNodes>(rRightHandSideVector, rVariables.PVector);

    // KRATOS_INFO("1-PwElement::CalculateAndAddPermeabilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                 ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-PwElement::CalculateAndAddFluidBodyFlow()") << std::endl;

    noalias(rVariables.PDimMatrix) =  prod(rVariables.GradNpT,rVariables.PermeabilityMatrix)
                                    * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) =  rVariables.DynamicViscosityInverse
                                 * rVariables.FluidDensity
                                 * rVariables.RelativePermeability
                                 * prod(rVariables.PDimMatrix,rVariables.BodyAcceleration);

    //Distribute fluid body flow block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.PVector);

    // KRATOS_INFO("1-PwElement::CalculateAndAddFluidBodyFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwElement<TDim,TNumNodes>::
    CalculateKinematics( ElementVariables& rVariables,
                         const unsigned int &PointNumber )

{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateKinematics") << std::endl;

    //Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Np) = row(rVariables.NContainer, PointNumber);
    noalias(rVariables.GradNpT) = rVariables.DN_DXContainer[PointNumber];

    rVariables.detJ0 = rVariables.detJContainer[PointNumber];

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateKinematics") << std::endl;

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------

template class PwElement<2,3>;
template class PwElement<2,4>;
template class PwElement<3,4>;
template class PwElement<3,8>;

template class PwElement<2,6>;
template class PwElement<2,8>;
template class PwElement<2,9>;
template class PwElement<3,10>;
template class PwElement<3,20>;
template class PwElement<3,27>;

} // Namespace Kratos
