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
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "custom_elements/geo_structural_base_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer GeoStructuralBaseElement<TDim,TNumNodes>::
    Create( IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << "calling the default Create method for a particular element ... illegal operation!!" << std::endl;

    return Element::Pointer( new GeoStructuralBaseElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer GeoStructuralBaseElement<TDim,TNumNodes>::
    Create( IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << "calling the default Create method for a particular element ... illegal operation!!" << std::endl;

    return Element::Pointer( new GeoStructuralBaseElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int GeoStructuralBaseElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

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

        if ( Geom[i].SolutionStepsDataHas( ROTATION ) == false )
            KRATOS_ERROR << "missing variable ROTATION on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].HasDofFor( DISPLACEMENT_X ) == false ||
             Geom[i].HasDofFor( DISPLACEMENT_Y ) == false ||
             Geom[i].HasDofFor( DISPLACEMENT_Z ) == false   )
            KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << Geom[i].Id() << std::endl;

        if ( Geom[i].HasDofFor( ROTATION_X ) == false ||
             Geom[i].HasDofFor( ROTATION_Y ) == false ||
             Geom[i].HasDofFor( ROTATION_Z ) == false )
            KRATOS_ERROR << "missing one of the dofs for the variable ROTATION on node " << Geom[i].Id() << std::endl;
    }

    // Verify ProcessInfo variables
    // Verify properties
    if ( Prop.Has( DENSITY ) == false ||
         Prop[DENSITY] < 0.0 )
        KRATOS_ERROR << "DENSITY has Key zero, is not defined or has an invalid value at element " << this->Id() << std::endl;

    if ( Prop.Has( YOUNG_MODULUS ) == false )
    {
        if ( Prop.Has( UDSM_NAME ) == false )
        {
            KRATOS_ERROR << "YOUNG_MODULUS has Key zero or is not defined at element " << this->Id() << std::endl;
        }
    }
    else
    {
        if ( Prop[YOUNG_MODULUS] <= 0.0 )
            KRATOS_ERROR << "YOUNG_MODULUS has an invalid value at element " << this->Id() << std::endl;
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

    // If this is a 2D problem, check that nodes are in XY plane
    if ( TDim == 2 )
    {
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            if (Geom[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
        }
    }

    return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const PropertiesType &Prop = this->GetProperties();
    const GeometryType &Geom = this->GetGeometry();
    // const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );
    const unsigned int NumGPoints = GetIntegrationPointsNumber();

    if ( mConstitutiveLawVector.size() != NumGPoints )
        mConstitutiveLawVector.resize( NumGPoints );

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        mConstitutiveLawVector[i] = Prop[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->InitializeMaterial( Prop,
                                                       Geom,
                                                       row( Geom.ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }

    // resize mStressVector:
    if ( mStressVector.size() != NumGPoints )
    {
       mStressVector.resize(NumGPoints);
       for (unsigned int i=0; i < mStressVector.size(); ++i)
       {
           mStressVector[i].resize(VoigtSize);
       }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    GetDofList( DofsVectorType& rElementalDofList,
                const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    if (rElementalDofList.size() != N_DOF_ELEMENT)
      rElementalDofList.resize( N_DOF_ELEMENT );

    unsigned int index = 0;
    if (TDim == 3)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[index++] = rGeom[i].pGetDof(ROTATION_X);
            rElementalDofList[index++] = rGeom[i].pGetDof(ROTATION_Y);
            rElementalDofList[index++] = rGeom[i].pGetDof(ROTATION_Z);
        }
    }
    else if (TDim == 2)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[index++] = rGeom[i].pGetDof(ROTATION_Z);
        }
    }
    else
    {
        KRATOS_ERROR << " Unspecified dimension in GetDofList: " << this->Id() << std::endl;
    }


    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
GeometryData::IntegrationMethod GeoStructuralBaseElement<TDim,TNumNodes>::
    GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                          VectorType& rRightHandSideVector,
                          const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != N_DOF_ELEMENT )
        rLeftHandSideMatrix.resize( N_DOF_ELEMENT, N_DOF_ELEMENT, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( N_DOF_ELEMENT, N_DOF_ELEMENT );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != N_DOF_ELEMENT )
        rRightHandSideVector.resize( N_DOF_ELEMENT, false );
    noalias( rRightHandSideVector ) = ZeroVector( N_DOF_ELEMENT );

    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                           const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_ERROR << "GeoStructuralBaseElement::CalculateLeftHandSide not implemented, element: " << this->Id() << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateRightHandSide( VectorType& rRightHandSideVector,
                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Resetting the RHS
    if ( rRightHandSideVector.size() != N_DOF_ELEMENT )
        rRightHandSideVector.resize( N_DOF_ELEMENT, false );
    noalias( rRightHandSideVector ) = ZeroVector( N_DOF_ELEMENT );

    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    EquationIdVector(EquationIdVectorType& rResult,
                     const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    unsigned int index = 0;

    if (rResult.size() != N_DOF_ELEMENT)
      rResult.resize( N_DOF_ELEMENT, false );

    if (TDim == 2)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index++] = rGeom[i].GetDof(ROTATION_Z).EquationId();
        }
    }
    else if (TDim == 3)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index++] = rGeom[i].GetDof(ROTATION_X).EquationId();
            rResult[index++] = rGeom[i].GetDof(ROTATION_Y).EquationId();
            rResult[index++] = rGeom[i].GetDof(ROTATION_Z).EquationId();
        }
    }
    else
    {
        KRATOS_ERROR << "undefined dimension in EquationIdVector... illegal operation!! element: " << this->Id() << std::endl;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateMassMatrix method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Rayleigh Method (Damping Matrix = alpha*M + beta*K)

    // Compute Mass Matrix
    MatrixType MassMatrix(N_DOF_ELEMENT, N_DOF_ELEMENT);

    this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

    // Compute Stiffness matrix
    MatrixType StiffnessMatrix(N_DOF_ELEMENT, N_DOF_ELEMENT);

    this->CalculateStiffnessMatrix(StiffnessMatrix, rCurrentProcessInfo);

    // Compute Damping Matrix
    if ( rDampingMatrix.size1() != N_DOF_ELEMENT )
        rDampingMatrix.resize( N_DOF_ELEMENT, N_DOF_ELEMENT, false );
    noalias( rDampingMatrix ) = ZeroMatrix( N_DOF_ELEMENT, N_DOF_ELEMENT );

    noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_ALPHA] * MassMatrix;
    noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_BETA] * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    GetValuesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();

    if ( rValues.size() != N_DOF_ELEMENT )
        rValues.resize( N_DOF_ELEMENT, false );

    unsigned int index = 0;
    if ( TDim > 2 )
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ROTATION_X,     Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ROTATION_Y,     Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ROTATION_Z,     Step );
        }
    }
    else
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ROTATION_Z,     Step );
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
GetFirstDerivativesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();

    if ( rValues.size() != N_DOF_ELEMENT )
        rValues.resize( N_DOF_ELEMENT, false );

    unsigned int index = 0;

    if ( TDim > 2 )
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ANGULAR_VELOCITY_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ANGULAR_VELOCITY_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
        }
    }
    else
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    GetSecondDerivativesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();

    if ( rValues.size() != N_DOF_ELEMENT )
        rValues.resize( N_DOF_ELEMENT, false );

    unsigned int index = 0;

    if ( TDim > 2 )
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
            
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
        }
    }
    else
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
            rValues[index++] = Geom[i].FastGetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                const std::vector<double>& rValues,
                                const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix,
                              const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateStiffnessMatrix method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateAll method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateRHS( VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateRHS method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    CalculateCrossDirection( Matrix &CrossDirection ) const
{
    KRATOS_TRY;

    KRATOS_ERROR << "calling the default CalculateCrossDirection method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
    InitializeElementVariables( ElementVariables& rVariables,
                                ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                const GeometryType& Geom,
                                const PropertiesType& Prop,
                                const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    //Properties variables

    //ProcessInfo variables

    //Nodal Variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector, Geom, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector,     Geom, VELOCITY);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration, Geom, VOLUME_ACCELERATION);

    rVariables.DofValuesVector.resize(N_DOF_ELEMENT*TNumNodes);
    GetNodalDofValuesVector(rVariables.DofValuesVector, Geom);

    //General Variables
    rVariables.DofValuesVector.resize(N_DOF_ELEMENT*TNumNodes);

    rVariables.CrossDirection.resize(TDim, TNumNodes);
    CalculateCrossDirection(rVariables.CrossDirection);

    //Variables computed at each GP
    rVariables.B.resize(VoigtSize, N_DOF_ELEMENT*TNumNodes, false);

    noalias(rVariables.Nu) = ZeroMatrix(TDim, TNumNodes*TDim);

    rVariables.TransformationMatrix.resize(TDim, VoigtSize, false);
    rVariables.UVoigtMatrix.resize(N_DOF_ELEMENT*TNumNodes, VoigtSize, false);

    //Constitutive Law parameters
    rVariables.StrainVector.resize(VoigtSize, false);
    rVariables.StressVector.resize(VoigtSize, false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);
    rVariables.Np.resize(TNumNodes, false);
    rVariables.GradNpT.resize(TNumNodes, TDim, false);
    rVariables.F.resize(TDim,TDim,false);
    rVariables.detF = 1.0;
    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetStressVector(rVariables.StressVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Np);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.GradNpT);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);
    rConstitutiveParameters.SetDeterminantF(rVariables.detF);
    //Auxiliary variables

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoStructuralBaseElement<TDim,TNumNodes>::
     GetNodalDofValuesVector(Vector &rNodalVariableVector,
                             const GeometryType &rGeom,
                             IndexType SolutionStepIndex) const
{
    unsigned int index = 0;

    if (TDim == 3)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X, SolutionStepIndex);
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y, SolutionStepIndex);
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Z, SolutionStepIndex);
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(ROTATION_X, SolutionStepIndex);
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(ROTATION_Y, SolutionStepIndex);
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(ROTATION_Z, SolutionStepIndex);
        }
    }
    else if (TDim == 2)
    {
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X, SolutionStepIndex);
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y, SolutionStepIndex);
            rNodalVariableVector[index++] = rGeom[i].FastGetSolutionStepValue(ROTATION_Z, SolutionStepIndex);
        }
    }
    else
    {
        KRATOS_ERROR << " Unspecified dimension in GetDofList: " << this->Id() << std::endl;
    }
}

//-------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
SizeType GeoStructuralBaseElement<TDim,TNumNodes>::
    GetIntegrationPointsNumber() const
{

    KRATOS_ERROR << "calling the default GetIntegrationPointsNumber method for a particular element ... illegal operation!!" << this->Id() << std::endl;

    return 0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class GeoStructuralBaseElement<2,3>;

template class GeoStructuralBaseElement<3,3>;
template class GeoStructuralBaseElement<3,4>;
template class GeoStructuralBaseElement<3,6>;
template class GeoStructuralBaseElement<3,8>;

} // Namespace Kratos
