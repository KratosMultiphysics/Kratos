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
#include "custom_elements/transient_Pw_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer TransientPwElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new TransientPwElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer TransientPwElement<TDim,TNumNodes>::
    Create(IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new TransientPwElement( NewId, pGeom, pProperties ) );
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
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
void TransientPwElement<TDim,TNumNodes>::
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
void TransientPwElement<TDim,TNumNodes>::
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
void TransientPwElement<TDim,TNumNodes>::
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
void TransientPwElement<TDim,TNumNodes>::
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
void TransientPwElement<TDim,TNumNodes>::
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
void TransientPwElement<TDim,TNumNodes>::
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
void TransientPwElement<TDim,TNumNodes>::
    Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientPwElement::Initialize()") << this->Id() << std::endl;

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

    // KRATOS_INFO("1-TransientPwElement::Initialize()") << std::endl;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int TransientPwElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientPwElement::Check()") << this->Id() << std::endl;

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();

    if (Geom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        if ( Geom[i].SolutionStepsDataHas( WATER_PRESSURE ) == false )
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].SolutionStepsDataHas( DT_WATER_PRESSURE ) == false )
            KRATOS_ERROR << "missing variable DT_WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false )
            KRATOS_ERROR << "missing variable VOLUME_ACCELERATION on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].HasDofFor( WATER_PRESSURE ) == false )
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;
    }

    // Verify ProcessInfo variables

    // Verify properties
    if ( Prop.Has( DENSITY_WATER ) == false || Prop[DENSITY_WATER] < 0.0 )
        KRATOS_ERROR << "DENSITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( BULK_MODULUS_SOLID ) == false || Prop[BULK_MODULUS_SOLID] < 0.0 )
        KRATOS_ERROR << "BULK_MODULUS_SOLID does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( POROSITY ) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0 )
        KRATOS_ERROR << "POROSITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;


    if ( TDim == 2 )
    {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            if (Geom[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
        }
    }

    // Verify specific properties
    if ( Prop.Has( BULK_MODULUS_FLUID ) == false || Prop[BULK_MODULUS_FLUID] < 0.0 )
        KRATOS_ERROR << "BULK_MODULUS_FLUID does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( DYNAMIC_VISCOSITY ) == false || Prop[DYNAMIC_VISCOSITY] < 0.0 )
        KRATOS_ERROR << "DYNAMIC_VISCOSITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( PERMEABILITY_XX ) == false || Prop[PERMEABILITY_XX] < 0.0 )
        KRATOS_ERROR << "PERMEABILITY_XX does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( PERMEABILITY_YY ) == false || Prop[PERMEABILITY_YY] < 0.0 )
        KRATOS_ERROR << "PERMEABILITY_YY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( PERMEABILITY_XY ) == false || Prop[PERMEABILITY_XY] < 0.0 )
        KRATOS_ERROR << "PERMEABILITY_XY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

    if (!Prop.Has( BIOT_COEFFICIENT ))
        KRATOS_ERROR << "BIOT_COEFFICIENT does not exist in the material properties in element" << this->Id() << std::endl;

    if (TDim > 2)
    {
        if ( Prop.Has( PERMEABILITY_ZZ ) == false || Prop[PERMEABILITY_ZZ] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_ZZ does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

        if ( Prop.Has( PERMEABILITY_YZ ) == false || Prop[PERMEABILITY_YZ] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_YZ does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;

        if ( Prop.Has( PERMEABILITY_ZX ) == false || Prop[PERMEABILITY_ZX] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_ZX does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
    }


    // Verify that the constitutive law has the correct dimension

    // Check constitutive law
    if ( mRetentionLawVector.size() > 0 ) {
        return mRetentionLawVector[0]->Check( Prop, rCurrentProcessInfo );
    }

    // KRATOS_INFO("1-TransientPwElement::Check()") << std::endl;

    return 0;

    KRATOS_CATCH( "" );
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientPwElement::InitializeSolutionStep()") << this->Id() << std::endl;

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

    // reset hydraulic discharge
    this->ResetHydraulicDischarge();

    // KRATOS_INFO("1-TransientPwElement::InitializeSolutionStep()") << std::endl;
    KRATOS_CATCH("");
}


//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // nothing

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // nothing

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientPwElement::FinalizeSolutionStep()") << std::endl;

    this->CalculateHydraulicDischarge(rCurrentProcessInfo);

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

    // KRATOS_INFO("1-TransientPwElement::FinalizeSolutionStep()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<double>& rVariable,
                                  std::vector<double>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientPwElement::CalculateOnIntegrationPoints()") << std::endl;

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
    
    // KRATOS_INFO("1-TransientPwElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable,
                                  std::vector<array_1d<double,3>>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientPwElement::CalculateOnIntegrationPoints<double,3>()") << std::endl;

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

    // KRATOS_INFO("1-TransientPwElement::CalculateOnIntegrationPoints<double,3>()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                 std::vector<Matrix>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientPwElement::CalculateOnIntegrationPoints()") << rVariable << std::endl;


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

    // KRATOS_INFO("1-TransientPwElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& rCurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientPwElement::CalculateAll()") << std::endl;

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
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJ0);

        //Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if (CalculateResidualVectorFlag)  this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    // KRATOS_INFO("1-TransientPwElement::CalculateAll()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    InitializeElementVariables( ElementVariables& rVariables,
                                const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-TransientPwElement::InitializeElementVariables()") << std::endl;

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

    // KRATOS_INFO("1-TransientPwElement::InitializeElementVariables()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientPwElement::CalculateAndAddLHS()") << std::endl;
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);

    // KRATOS_INFO("1-TransientPwElement::CalculateAndAddLHS()") << std::endl;

    KRATOS_CATCH("");
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientPwElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    this->CalculateCompressibilityMatrix(rVariables.PMatrix, rVariables);

    //Distribute compressibility block matrix into the elemental matrix
    GeoElementUtilities::
        AssemblePBlockMatrix< 0, TNumNodes >(rLeftHandSideMatrix,
                                             rVariables.PMatrix);

    // KRATOS_INFO("1-TransientPwElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix,
                                      ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientPwElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    this->CalculatePermeabilityMatrix(rVariables.PDimMatrix,
                                      rVariables.PMatrix,
                                      rVariables);

    //Distribute permeability block matrix into the elemental matrix
    GeoElementUtilities::AssemblePBlockMatrix< 0, TNumNodes >(rLeftHandSideMatrix, rVariables.PMatrix);

    // KRATOS_INFO("1-TransientPwElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientPwElement::CalculateAndAddRHS()") << std::endl;

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    // KRATOS_INFO("1-TransientPwElement::CalculateAndAddRHS()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityFlow( VectorType& rRightHandSideVector,
                                     ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientPwElement::CalculateAndAddPermeabilityFlow()") << std::endl;

    this->CalculatePermeabilityFlow(rVariables.PDimMatrix,
                                    rVariables.PMatrix,
                                    rVariables.PVector,
                                    rVariables);

    //Distribute permeability block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.PVector);

    // KRATOS_INFO("1-TransientPwElement::CalculateAndAddPermeabilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                 ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientPwElement::CalculateAndAddFluidBodyFlow()") << std::endl;

    this->CalculateFluidBodyFlow(rVariables.PDimMatrix,
                                 rVariables.PVector,
                                 rVariables);

    //Distribute fluid body flow block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector,rVariables.PVector);

    // KRATOS_INFO("1-TransientPwElement::CalculateAndAddFluidBodyFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
    CalculateAndAddCompressibilityFlow( VectorType& rRightHandSideVector,
                                        ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-TransientPwElement::CalculateAndAddCompressibilityFlow()") << std::endl;

    this->CalculateCompressibilityFlow(rVariables.PMatrix,
                                       rVariables.PVector,
                                       rVariables);

    //Distribute compressibility block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.PVector);

    // KRATOS_INFO("1-TransientPwElement::CalculateAndAddCompressibilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void TransientPwElement<TDim,TNumNodes>::
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
unsigned int TransientPwElement<TDim,TNumNodes>::GetNumberOfDOF() const
{
    return TNumNodes;
}

//----------------------------------------------------------------------------------------------------

template class TransientPwElement<2,3>;
template class TransientPwElement<2,4>;
template class TransientPwElement<3,4>;
template class TransientPwElement<3,8>;

template class TransientPwElement<2,6>;
template class TransientPwElement<2,8>;
template class TransientPwElement<2,9>;
template class TransientPwElement<3,10>;
template class TransientPwElement<3,20>;
template class TransientPwElement<3,27>;

} // Namespace Kratos
