// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// Application includes
#include "custom_elements/U_Pw_base_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwBaseElement<TDim,TNumNodes>::Create(IndexType NewId, 
                                                    NodesArrayType const& ThisNodes,
                                                    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "calling the default Create method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    return Element::Pointer( new UPwBaseElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwBaseElement<TDim,TNumNodes>::Create(IndexType NewId,
                                                    GeometryType::Pointer pGeom,
                                                    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "calling the default Create method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    return Element::Pointer( new UPwBaseElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int UPwBaseElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();

    // verify nodal variables and dofs

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        if ( Geom[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_ERROR << "missing variable DISPLACEMENT on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].SolutionStepsDataHas( VELOCITY ) == false )
            KRATOS_ERROR << "missing variable VELOCITY on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].SolutionStepsDataHas( ACCELERATION ) == false )
            KRATOS_ERROR << "missing variable ACCELERATION on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].SolutionStepsDataHas( WATER_PRESSURE ) == false )
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].SolutionStepsDataHas( DT_WATER_PRESSURE ) == false )
            KRATOS_ERROR << "missing variable DT_WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false )
            KRATOS_ERROR << "missing variable VOLUME_ACCELERATION on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].HasDofFor( DISPLACEMENT_X ) == false ||
             Geom[i].HasDofFor( DISPLACEMENT_Y ) == false ||
             Geom[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].HasDofFor( WATER_PRESSURE ) == false )
            KRATOS_ERROR << "missing the dof for the variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;
    }

    // Verify ProcessInfo variables

    // Verify properties
    if ( Prop.Has( DENSITY_SOLID ) == false || Prop[DENSITY_SOLID] < 0.0 )
        KRATOS_ERROR << "DENSITY_SOLID has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( DENSITY_WATER ) == false || Prop[DENSITY_WATER] < 0.0 )
        KRATOS_ERROR << "DENSITY_WATER has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( YOUNG_MODULUS ) == false )
    {
        if ( Prop.Has( UDSM_NAME ) == false )
        {
            KRATOS_ERROR << "YOUNG_MODULUS has Key zero or is not defined at element" << this->Id() << std::endl;
        }
    } 
    else
    {
        if ( Prop[YOUNG_MODULUS] <= 0.0 )
            KRATOS_ERROR << "YOUNG_MODULUS has an invalid value at element" << this->Id() << std::endl;
    }

    if ( Prop.Has( POISSON_RATIO ) == false )
    {
        if ( Prop.Has( UDSM_NAME ) == false )
        {
            KRATOS_ERROR << "POISSON_RATIO has Key zero or is not defined at element" << this->Id() << std::endl;
        } 
    }
    else
    {
        const double& PoissonRatio = Prop[POISSON_RATIO];
        if ( PoissonRatio < 0.0 || PoissonRatio >= 0.5 )
            KRATOS_ERROR << "POISSON_RATIO has an invalid value at element" << this->Id() << std::endl;
    }

    if ( Prop.Has( BULK_MODULUS_SOLID ) == false || Prop[BULK_MODULUS_SOLID] < 0.0 )
        KRATOS_ERROR << "BULK_MODULUS_SOLID has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

    if ( Prop.Has( POROSITY ) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0 )
        KRATOS_ERROR << "POROSITY has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

    if ( TDim == 2 ) {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i=0; i<TNumNodes; ++i) {
            if (Geom[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
        }
    }

    return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwBaseElement::Initialize()") << this->Id() << std::endl;

    const PropertiesType &Prop = this->GetProperties();
    const GeometryType &Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    // pointer to constitutive laws
    if ( mConstitutiveLawVector.size() != NumGPoints )
        mConstitutiveLawVector.resize( NumGPoints );

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        mConstitutiveLawVector[i] = Prop[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->
            InitializeMaterial( Prop,
                                Geom,
                                row( Geom.ShapeFunctionsValues( this->GetIntegrationMethod() ), i ) );
    }

    // resize mStressVector:
    if ( mStressVector.size() != NumGPoints )
    {
       unsigned int VoigtSize = VOIGT_SIZE_3D;
       if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRAIN;
       mStressVector.resize(NumGPoints);
       for (unsigned int i=0; i < mStressVector.size(); ++i)
       {
           mStressVector[i].resize(VoigtSize);
           std::fill(mStressVector[i].begin(), mStressVector[i].end(), 0.0);
       }
    }

    // resizing and setting state variables
    if (mStateVariablesFinalized.size() != NumGPoints)
       mStateVariablesFinalized.resize(NumGPoints);
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
    {
        int nStateVariables = 0;
        nStateVariables = mConstitutiveLawVector[i]->GetValue( NUMBER_OF_UMAT_STATE_VARIABLES,
                                                               nStateVariables);
        if (nStateVariables > 0)
        {
            //ProcessInfo rCurrentProcessInfo;
            mConstitutiveLawVector[i]->SetValue( STATE_VARIABLES,
                                                 mStateVariablesFinalized[i],
                                                 rCurrentProcessInfo );
        }
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

    // KRATOS_INFO("1-UPwBaseElement::Initialize()") << std::endl;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    ResetConstitutiveLaw()
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwBaseElement::ResetConstitutiveLaw()") << this->Id() << std::endl;

    // erasing stress vectors
    for (unsigned int i=0; i < mStressVector.size(); ++i)
    {
        mStressVector[i].clear();
    }
    mStressVector.clear();

    for (unsigned int i=0; i < mStateVariablesFinalized.size(); ++i)
    {
        mStateVariablesFinalized[i].clear();
    }
    mStateVariablesFinalized.clear();


    KRATOS_CATCH( "" )

    // KRATOS_INFO("1-UPwBaseElement::ResetConstitutiveLaw()") << this->Id() << std::endl;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
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
        rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        if (TDim>2)
            rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index++] = rGeom[i].pGetDof(WATER_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
GeometryData::IntegrationMethod UPwBaseElement<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
    //return GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                          VectorType& rRightHandSideVector,
                          const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != N_DOF )
        rLeftHandSideMatrix.resize( N_DOF, N_DOF, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( N_DOF, N_DOF );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != N_DOF )
        rRightHandSideVector.resize( N_DOF, false );
    noalias( rRightHandSideVector ) = ZeroVector( N_DOF );

    //calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix,
                 rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                           const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = false;
    VectorType TempVector;

    CalculateAll(rLeftHandSideMatrix,
                 TempVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    CalculateRightHandSide( VectorType& rRightHandSideVector,
                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    //Resetting the RHS
    if ( rRightHandSideVector.size() != N_DOF )
        rRightHandSideVector.resize( N_DOF, false );
    noalias( rRightHandSideVector ) = ZeroVector( N_DOF );

    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType TempMatrix = Matrix();

    CalculateAll(TempMatrix,
                 rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    EquationIdVector(EquationIdVectorType& rResult,
                     const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int N_DOF = this->GetNumberOfDOF();
    unsigned int index = 0;

    if (rResult.size() != N_DOF)
      rResult.resize( N_DOF, false );

    if (TDim == 2)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
        }
    }
    else if (TDim == 3)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
        }
    }
    else
    {
        KRATOS_ERROR << "undefined dimension in EquationIdVector... illegal operation!!" << this->Id() << std::endl;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    CalculateMassMatrix( MatrixType& rMassMatrix,
                         const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateMassMatrix method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    CalculateDampingMatrix(MatrixType& rDampingMatrix,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Rayleigh Method (Damping Matrix = alpha*M + beta*K)

    const unsigned int N_DOF = this->GetNumberOfDOF();

    // Compute Mass Matrix
    MatrixType MassMatrix(N_DOF, N_DOF);

    this->CalculateMassMatrix(MassMatrix,rCurrentProcessInfo);

    // Compute Stiffness matrix
    MatrixType StiffnessMatrix(N_DOF, N_DOF);

    this->CalculateMaterialStiffnessMatrix(StiffnessMatrix, rCurrentProcessInfo);

    // Compute Damping Matrix
    if ( rDampingMatrix.size1() != N_DOF )
        rDampingMatrix.resize( N_DOF, N_DOF, false );
    noalias( rDampingMatrix ) = ZeroMatrix( N_DOF, N_DOF );

    const PropertiesType& Prop = this->GetProperties();

    if (Prop.Has( RAYLEIGH_ALPHA ))
        noalias(rDampingMatrix) += Prop[RAYLEIGH_ALPHA] * MassMatrix;
    else
        noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_ALPHA] * MassMatrix;

    if (Prop.Has( RAYLEIGH_BETA ))
        noalias(rDampingMatrix) += Prop[RAYLEIGH_BETA] * StiffnessMatrix;
    else
        noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_BETA] * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    GetValuesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();
    const unsigned int N_DOF = this->GetNumberOfDOF();

    if ( rValues.size() != N_DOF )
        rValues.resize( N_DOF, false );

    if ( TDim > 2 )
    {
        unsigned int index = 0;
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index++] = 0.0;
        }
    }
    else
    {
        unsigned int index = 0;
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );
            rValues[index++] = 0.0;
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    GetFirstDerivativesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();
    const unsigned int N_DOF = this->GetNumberOfDOF();

    if ( rValues.size() != N_DOF )
        rValues.resize( N_DOF, false );


    if ( TDim > 2 )
    {
        unsigned int index = 0;
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index++] = 0.0;
        }
    }
    else
    {
        unsigned int index = 0;
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
            rValues[index++] = 0.0;
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    GetSecondDerivativesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();
    const unsigned int N_DOF = this->GetNumberOfDOF();

    if ( rValues.size() != N_DOF )
        rValues.resize( N_DOF, false );

    unsigned int index = 0;

    if ( TDim > 2 )
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
            rValues[index++] = 0.0;
        }

    }
    else
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
            rValues[index++] = 0.0;
        }
    }

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                 const std::vector<Vector>& rValues,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                 const std::vector<Matrix>& rValues,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );

    KRATOS_CATCH( "" )
}


//-------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                 const std::vector<double>& rValues,
                                 const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );

    KRATOS_CATCH( "" )
}



//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                 std::vector<ConstitutiveLaw::Pointer>& rValues,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == CONSTITUTIVE_LAW)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for (unsigned int i=0; i < mConstitutiveLawVector.size(); i++)
            rValues[i] = mConstitutiveLawVector[i];
    }

    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    CalculateMaterialStiffnessMatrix( MatrixType& rStiffnessMatrix,
                                      const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateMaterialStiffnessMatrix method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwBaseElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateAll method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwBaseElement<TDim,TNumNodes>::
    CalculateIntegrationCoefficient(const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                                    const IndexType& PointNumber,
                                    const double& detJ)

{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Matrix& UPwBaseElement<TDim,TNumNodes>::
    CalculateDeltaDisplacement(Matrix& DeltaDisplacement) const
{
    KRATOS_TRY

    DeltaDisplacement.resize(TNumNodes , TDim, false);

    for ( IndexType i_node = 0; i_node < TNumNodes; i_node++ ) {
        const array_1d<double, 3>& current_displacement  = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3>& previous_displacement = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( IndexType iDim = 0; iDim < TDim; ++iDim )
            DeltaDisplacement(i_node, iDim) = current_displacement[iDim] - previous_displacement[iDim];
    }

    return DeltaDisplacement;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwBaseElement<TDim,TNumNodes>::
    CalculateDerivativesOnInitialConfiguration(const GeometryType& Geometry,
                                               Matrix& DNu_DX0,
                                               const IndexType& GPoint,
                                               IntegrationMethod ThisIntegrationMethod) const
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

    Matrix J0, InvJ0;
    double detJ0;
    GeometryUtils::JacobianOnInitialConfiguration(Geometry, IntegrationPoints[GPoint], J0);
    const Matrix& DN_De = Geometry.ShapeFunctionsLocalGradients(ThisIntegrationMethod)[GPoint];
    MathUtils<double>::InvertMatrix( J0, InvJ0, detJ0 );
    GeometryUtils::ShapeFunctionsGradients(DN_De, InvJ0, DNu_DX0);

    return detJ0;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
unsigned int UPwBaseElement<TDim,TNumNodes>::GetNumberOfDOF() const
{
    return TNumNodes * (TDim + 1);
}

//-------------------------------------------------------------------------------------------------------------------------------------------

template class UPwBaseElement<2,3>;
template class UPwBaseElement<2,4>;
template class UPwBaseElement<3,4>;
template class UPwBaseElement<3,6>;
template class UPwBaseElement<3,8>;

template class UPwBaseElement<2,6>;
template class UPwBaseElement<2,8>;
template class UPwBaseElement<2,9>;
template class UPwBaseElement<3,10>;
template class UPwBaseElement<3,20>;
template class UPwBaseElement<3,27>;

} // Namespace Kratos
