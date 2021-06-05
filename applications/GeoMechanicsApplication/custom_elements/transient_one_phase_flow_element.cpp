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
#include "custom_elements/transient_one_phase_flow_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer TransientOnePhaseFlowElement<TDim,TNumNodes>::Create(IndexType NewId,
                                                               NodesArrayType const& ThisNodes,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new TransientOnePhaseFlowElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer TransientOnePhaseFlowElement<TDim,TNumNodes>::Create(IndexType NewId,
                                                               GeometryType::Pointer pGeom,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new TransientOnePhaseFlowElement( NewId, pGeom, pProperties ) );
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    GetDofList( DofsVectorType& rElementalDofList,
                const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int N_DOF = this->GetNumberOfDOF();
    unsigned int index = 0;

    if (rElementalDofList.size() != N_DOF)
      rElementalDofList.resize( N_DOF );

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[index++] = rGeom[i].pGetDof(WATER_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    EquationIdVector(EquationIdVectorType& rResult,
                     const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int N_DOF = this->GetNumberOfDOF();
    unsigned int index = 0;

    if (rResult.size() != N_DOF)
      rResult.resize( N_DOF, false );

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateMassMatrix( MatrixType& rMassMatrix,
                         const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    //Resizing mass matrix
    if ( rMassMatrix.size1() != N_DOF )
        rMassMatrix.resize( N_DOF, N_DOF, false );
    noalias( rMassMatrix ) = ZeroMatrix( N_DOF, N_DOF );

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateDampingMatrix(MatrixType& rDampingMatrix,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Rayleigh Method (Damping Matrix = alpha*M + beta*K)

    const unsigned int N_DOF = this->GetNumberOfDOF();

    // Compute Damping Matrix
    if ( rDampingMatrix.size1() != N_DOF )
        rDampingMatrix.resize( N_DOF, N_DOF, false );
    noalias( rDampingMatrix ) = ZeroMatrix( N_DOF, N_DOF );

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    GetValuesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if ( rValues.size() != N_DOF )
        rValues.resize( N_DOF, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    GetFirstDerivativesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if ( rValues.size() != N_DOF )
        rValues.resize( N_DOF, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    GetSecondDerivativesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    if ( rValues.size() != N_DOF )
        rValues.resize( N_DOF, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[i] = 0.0;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::Initialize()") << this->Id() << std::endl;

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

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::Initialize()") << std::endl;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int TransientOnePhaseFlowElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::Check()") << this->Id() << std::endl;

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
        KRATOS_THROW_ERROR( std::invalid_argument,
            "DENSITY_SOLID does not exist in the material properties or has an invalid value at element", this->Id() )
    if ( Prop.Has( DENSITY_WATER ) == false || Prop[DENSITY_WATER] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
            "DENSITY_WATER does not exist in the material properties or has an invalid value at element", this->Id() )


    if ( Prop.Has( BULK_MODULUS_SOLID ) == false || Prop[BULK_MODULUS_SOLID] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
            "BULK_MODULUS_SOLID does not exist in the material properties or has an invalid value at element", this->Id() )

    if ( Prop.Has( POROSITY ) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
            "POROSITY does not exist in the material properties or has an invalid value at element", this->Id() )


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
                            "BULK_MODULUS_FLUID does not exist in the material properties or has an invalid value at element",
                            this->Id() )

    if ( Prop.Has( DYNAMIC_VISCOSITY ) == false || Prop[DYNAMIC_VISCOSITY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "DYNAMIC_VISCOSITY does not exist in the material properties or has an invalid value at element",
                            this->Id() )

    if ( Prop.Has( PERMEABILITY_XX ) == false || Prop[PERMEABILITY_XX] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "PERMEABILITY_XX does not exist in the material properties or has an invalid value at element",
                            this->Id() )

    if ( Prop.Has( PERMEABILITY_YY ) == false || Prop[PERMEABILITY_YY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "PERMEABILITY_YY does not exist in the material properties or has an invalid value at element",
                            this->Id() )

    if ( Prop.Has( PERMEABILITY_XY ) == false || Prop[PERMEABILITY_XY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "PERMEABILITY_XY does not exist in the material properties or has an invalid value at element",
                            this->Id() )

    if (!Prop.Has( BIOT_COEFFICIENT ))
        KRATOS_THROW_ERROR( std::invalid_argument,
                            "BIOT_COEFFICIENT does not exist in the material properties in element",
                            this->Id() )

    if (TDim > 2)
    {
        if ( Prop.Has( PERMEABILITY_ZZ ) == false || Prop[PERMEABILITY_ZZ] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_ZZ does not exist in the material properties or has an invalid value at element", this->Id() )

        if ( Prop.Has( PERMEABILITY_YZ ) == false || Prop[PERMEABILITY_YZ] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_YZ does not exist in the material properties or has an invalid value at element", this->Id() )

        if ( Prop.Has( PERMEABILITY_ZX ) == false || Prop[PERMEABILITY_ZX] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_ZX does not exist in the material properties or has an invalid value at element", this->Id() )
    }


    // Verify that the constitutive law has the correct dimension

    // Check constitutive law
    if ( mRetentionLawVector.size() > 0 ) {
        return mRetentionLawVector[0]->Check( Prop, rCurrentProcessInfo );
    }

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::Check()") << std::endl;

    return 0;

    KRATOS_CATCH( "" );
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::InitializeSolutionStep()") << this->Id() << std::endl;

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

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::InitializeSolutionStep()") << std::endl;
    KRATOS_CATCH("");
}


//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // nothing

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // nothing

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::FinalizeSolutionStep()") << std::endl;

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

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::FinalizeSolutionStep()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<double>& rVariable,
                                  std::vector<double>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateOnIntegrationPoints()") << std::endl;

    if (rVariable == DEGREE_OF_SATURATION     ||
        rVariable == EFFECTIVE_SATURATION     ||
        rVariable == BISHOP_COEFICIENT        ||
        rVariable == DERIVATIVE_OF_SATURATION ||
        rVariable == RELATIVE_PERMEABILITY    ||
        rVariable == HYDRAULIC_HEAD)
    {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
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
    
    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable,
                                  std::vector<array_1d<double,3>>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateOnIntegrationPoints<double,3>()") << std::endl;

    if (rVariable == FLUID_FLUX_VECTOR)
    {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
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

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateOnIntegrationPoints<double,3>()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                 std::vector<Matrix>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateOnIntegrationPoints()") << rVariable << std::endl;


    if (rVariable == PERMEABILITY_MATRIX)
    {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
    else
    {
        if ( rOutput.size() != mRetentionLawVector.size() )
            rOutput.resize(mRetentionLawVector.size());

        for ( unsigned int i = 0;  i < mRetentionLawVector.size(); i++ )
        {
            rOutput[i].resize(TDim,TDim,false);
            noalias(rOutput[i]) = ZeroMatrix(TDim,TDim);
            rOutput[i] = mRetentionLawVector[i]->GetValue( rVariable, rOutput[i] );
        }
    }

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& rCurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateAll()") << std::endl;

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

    const bool hasBiotCoefficient = Prop.Has(BIOT_COEFFICIENT);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        //Compute Nu and BodyAcceleration
        GeoElementUtilities::
            CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>( Variables.BodyAcceleration,
                                                                Variables.NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient,
                                              Variables.detJ0,
                                              IntegrationPoints[GPoint].Weight());

        //Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if (CalculateResidualVectorFlag)  this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateAll()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    InitializeElementVariables( ElementVariables& rVariables,
                                const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::InitializeElementVariables()") << std::endl;

    //Properties variables
    this->InitializeProperties( rVariables );

    //ProcessInfo variables
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    //Nodal Variables
    this->InitializeNodalPorePressureVariables( rVariables );
    this->InitializeNodalVolumeAccelerationVariables( rVariables );

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

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::InitializeElementVariables()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateAndAddLHS()") << std::endl;
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateAndAddLHS()") << std::endl;

    KRATOS_CATCH("");
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    this->CalculateCompressibilityMatrix(rVariables.PMatrix, rVariables);

    //Distribute compressibility block matrix into the elemental matrix
    GeoElementUtilities::
        AssemblePBlockMatrix< 0, TNumNodes >(rLeftHandSideMatrix,
                                             rVariables.PMatrix);

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix,
                                      ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    this->CalculatePermeabilityMatrix(rVariables.PDimMatrix,
                                      rVariables.PMatrix,
                                      rVariables);

    //Distribute permeability block matrix into the elemental matrix
    GeoElementUtilities::AssemblePBlockMatrix< 0, TNumNodes >(rLeftHandSideMatrix, rVariables.PMatrix);

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateAndAddRHS()") << std::endl;

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    // VG: TODO: Check
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateAndAddRHS()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityFlow( VectorType& rRightHandSideVector,
                                     ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateAndAddPermeabilityFlow()") << std::endl;

    this->CalculatePermeabilityFlow(rVariables.PDimMatrix,
                                    rVariables.PMatrix,
                                    rVariables.PVector,
                                    rVariables);

    //Distribute permeability block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.PVector);

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateAndAddPermeabilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                 ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateAndAddFluidBodyFlow()") << std::endl;

    this->CalculateFluidBodyFlow(rVariables.PDimMatrix,
                                 rVariables.PVector,
                                 rVariables);

    //Distribute fluid body flow block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector,rVariables.PVector);

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateAndAddFluidBodyFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
    CalculateAndAddCompressibilityFlow( VectorType& rRightHandSideVector,
                                        ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientOnePhaseFlowElement::CalculateAndAddCompressibilityFlow()") << std::endl;

    this->CalculateCompressibilityFlow(rVariables.PMatrix,
                                       rVariables.PVector,
                                       rVariables);

    //Distribute compressibility block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.PVector);

    // KRATOS_INFO("1-TransientOnePhaseFlowElement::CalculateAndAddCompressibilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientOnePhaseFlowElement<TDim,TNumNodes>::
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

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
unsigned int TransientOnePhaseFlowElement<TDim,TNumNodes>::GetNumberOfDOF() const
{
    return TNumNodes;
}

//----------------------------------------------------------------------------------------------------

template class TransientOnePhaseFlowElement<2,3>;
template class TransientOnePhaseFlowElement<2,4>;
template class TransientOnePhaseFlowElement<3,4>;
template class TransientOnePhaseFlowElement<3,8>;

template class TransientOnePhaseFlowElement<2,6>;
template class TransientOnePhaseFlowElement<2,8>;
template class TransientOnePhaseFlowElement<2,9>;
template class TransientOnePhaseFlowElement<3,10>;
template class TransientOnePhaseFlowElement<3,20>;
template class TransientOnePhaseFlowElement<3,27>;

} // Namespace Kratos
