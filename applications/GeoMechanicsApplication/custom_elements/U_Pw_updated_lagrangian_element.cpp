// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz,
//                   Vahid Galavi
//

// External includes

// Project includes
#include "custom_elements/U_Pw_updated_lagrangian_element.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwUpdatedLagrangianElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           GeometryType::Pointer pGeom,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwUpdatedLagrangianElement( NewId, pGeom, pProperties ) );
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    // Verify generic variables
    int ierr = UPwSmallStrainElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);

    return ierr;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    UPwSmallStrainElement<TDim,TNumNodes>::Initialize(rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType &IntegrationPoints =
        this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType NumGPoints = IntegrationPoints.size();

    if ( mDetF0.size() != NumGPoints)
        mDetF0.resize( NumGPoints );

    if ( mF0.size() != NumGPoints)
        mF0.resize( NumGPoints );

    for (IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint)
    {
        mDetF0[GPoint] = 1.0;
        mF0[GPoint] = IdentityMatrix(TDim);
    }

    mF0Computed = false;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    InitializeSolutionStep(const ProcessInfo &rCurrentProcessInfo)
{
    UPwSmallStrainElement<TDim, TNumNodes>::InitializeSolutionStep(rCurrentProcessInfo);

    mF0Computed = false;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(this->GetGeometry(),
                                                       this->GetProperties(),
                                                       rCurrentProcessInfo);
    //ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);


    ElementVariables Variables;
    UPwSmallStrainElement<TDim,TNumNodes>::InitializeElementVariables( Variables,
                                                                       ConstitutiveParameters,
                                                                       this->GetGeometry(),
                                                                       this->GetProperties(),
                                                                       rCurrentProcessInfo );


    // Reading integration points
    for ( IndexType GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
        // Compute element kinematics B, F, GradNpT ...
        this->CalculateKinematics(Variables, GPoint);

        // Call the constitutive law to update material variables
        //Compute constitutive tensor and stresses
        UpdateElementalVariableStressVector(Variables, GPoint);
        mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
        mStateVariablesFinalized[GPoint] = 
            mConstitutiveLawVector[GPoint]->GetValue( STATE_VARIABLES,
                                                      mStateVariablesFinalized[GPoint] );

        // Update the element internal variables
        this->UpdateHistoricalDatabase(Variables, GPoint);
    }

    mF0Computed = true;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    UpdateHistoricalDatabase(ElementVariables& rVariables,
                             const SizeType GPoint)
{
    mDetF0[GPoint] = rVariables.detF;
    noalias(mF0[GPoint]) = rVariables.F;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& rCurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &IntegrationPoints =
        this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    //Containers of variables at all integration points
    const Matrix& NContainer = this->GetGeometry().ShapeFunctionsValues( this->GetIntegrationMethod() );

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(this->GetGeometry(),
                                                       this->GetProperties(),
                                                       rCurrentProcessInfo);
    if (CalculateStiffnessMatrixFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)  ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);


    ElementVariables Variables;
    UPwSmallStrainElement<TDim,TNumNodes>::InitializeElementVariables( Variables,
                                                                       ConstitutiveParameters,
                                                                       this->GetGeometry(),
                                                                       this->GetProperties(),
                                                                       rCurrentProcessInfo );


    // Computing in all integrations points
    for ( IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint )
    {
        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, NContainer, GPoint);
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>( Variables.BodyAcceleration,
                                                                NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );


        // Compute element kinematics B, F, GradNpT ...
        this->CalculateKinematics(Variables, GPoint);

        // Cauchy strain: This needs to be investigated which strain measure should be used
        // In some references, e.g. Bathe, suggested to use Almansi strain measure
        noalias(Variables.StrainVector) = prod(Variables.B, Variables.DisplacementVector);

        //Compute constitutive tensor and stresses
        UpdateElementalVariableStressVector(Variables, GPoint);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        UpdateStressVector(Variables, GPoint);

        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = CalculateBulkModulus(Variables.ConstitutiveMatrix);
        this->InitializeBiotCoefficients(Variables, BulkModulus);

        // Calculating weights for integration on the reference configuration
        this->CalculateIntegrationCoefficient( Variables.IntegrationCoefficient,
                                               Variables.detJ0,
                                               IntegrationPoints[GPoint].Weight() );

        if ( CalculateStiffnessMatrixFlag == true )
        {
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

            /* Geometric stiffness matrix */
            this->CalculateAndAddGeometricStiffnessMatrix( rLeftHandSideMatrix, Variables );
        }

        if (CalculateResidualVectorFlag)
        {
            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateAndAddGeometricStiffnessMatrix( MatrixType& rLeftHandSideMatrix,
                                             ElementVariables& rVariables )
{
    KRATOS_TRY

    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix ReducedKgMatrix = prod( rVariables.GradNpT,
                                   rVariables.IntegrationCoefficient *
                                   Matrix( prod( StressTensor, trans(rVariables.GradNpT) ) ) ); //to be optimized

    Matrix UMatrix(TNumNodes*TDim, TNumNodes*TDim);
    noalias(UMatrix) = ZeroMatrix(TNumNodes*TDim, TNumNodes*TDim);
    MathUtils<double>::ExpandAndAddReducedMatrix( UMatrix, ReducedKgMatrix, TDim );

    //Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUBlockMatrix(rLeftHandSideMatrix, UMatrix, TNumNodes, TDim);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateKinematics(ElementVariables& rVariables,
                        const SizeType GPoint)
{
    rVariables.detJ0 =
        this->CalculateDerivativesOnReferenceConfiguration(rVariables.J0,
                                                           rVariables.InvJ0,
                                                           rVariables.GradNpT,
                                                           GPoint,
                                                           this->GetIntegrationMethod());

    // Calculating jacobian
    Matrix J, inv_J;
    rVariables.detJ0 =
        this->CalculateDerivativesOnCurrentConfiguration(J,
                                                         inv_J,
                                                         rVariables.GradNpT,
                                                         GPoint,
                                                         this->GetIntegrationMethod());

    KRATOS_ERROR_IF(rVariables.detJ0 < 0.0)
     << "ERROR:: ELEMENT ID: "
     << this->Id()
     << " INVERTED. DETJ0: "
     << rVariables.detJ0
     << std::endl;

    // Deformation gradient
    Matrix DF = prod( J, rVariables.InvJ0 );

    const double detDF = MathUtils<double>::Det(DF);
    rVariables.detF = detDF * this->ReferenceConfigurationDeformationGradientDeterminant(GPoint);
    noalias(rVariables.F) = prod(DF, this->ReferenceConfigurationDeformationGradient(GPoint));

    // Calculating operator B
    this->CalculateBMatrix( rVariables.B, rVariables.GradNpT);
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateDerivativesOnReferenceConfiguration(Matrix& J0,
                                                 Matrix& InvJ0,
                                                 Matrix& GradNpT,
                                                 const IndexType &GPoint,
                                                 IntegrationMethod ThisIntegrationMethod) const
{
    J0.clear();

    double detJ0;

    Matrix deltaDisplacement;
    deltaDisplacement = this->CalculateDeltaDisplacement(deltaDisplacement);

    J0 = this->GetGeometry().Jacobian(J0, GPoint, ThisIntegrationMethod, deltaDisplacement);

    const Matrix& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[GPoint];

    MathUtils<double>::InvertMatrix( J0, InvJ0, detJ0 );

    noalias( GradNpT ) = prod( DN_De, InvJ0);

    return detJ0;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
    double UPwUpdatedLagrangianElement<TDim,TNumNodes>::
        CalculateDerivativesOnCurrentConfiguration( Matrix& rJ,
                                                    Matrix& rInvJ,
                                                    Matrix& rDN_DX,
                                                    const IndexType &PointNumber,
                                                    IntegrationMethod ThisIntegrationMethod ) const
{
    double detJ;
    rJ = this->GetGeometry().Jacobian( rJ, PointNumber, ThisIntegrationMethod );
    const Matrix& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    MathUtils<double>::InvertMatrix( rJ, rInvJ, detJ );
    GeometryUtils::ShapeFunctionsGradients(DN_De, rInvJ, rDN_DX);
    return detJ;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Matrix& UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateDeltaDisplacement(Matrix& DeltaDisplacement) const
{
    KRATOS_TRY

    DeltaDisplacement.resize(TNumNodes , TDim, false);

    for ( IndexType iNode = 0; iNode < TNumNodes; iNode++ ) {
        const array_1d<double, 3>& currentDisplacement  = this->GetGeometry()[iNode].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3>& previousDisplacement = this->GetGeometry()[iNode].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( IndexType iDim = 0; iDim < TDim; ++iDim )
            DeltaDisplacement(iNode, iDim) = currentDisplacement[iDim] - previousDisplacement[iDim];
    }

    return DeltaDisplacement;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    ReferenceConfigurationDeformationGradientDeterminant(const IndexType GPoint) const
{
    if (mF0Computed == false)
        return mDetF0[GPoint];

    return 1.0;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Matrix UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    ReferenceConfigurationDeformationGradient(const IndexType GPoint) const
{
    if (mF0Computed == false)
        return mF0[GPoint];

    return IdentityMatrix(TDim);
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                 std::vector<double>& rValues,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        if (rValues.size() != mConstitutiveLawVector.size())
            rValues.resize(mConstitutiveLawVector.size());

        for ( IndexType GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint )
            rValues[GPoint] = mDetF0[GPoint];
    } else {
        UPwSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,
                                                                            rValues,
                                                                            rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                 std::vector<Matrix>& rValues,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
        if (rValues.size() != mConstitutiveLawVector.size())
            rValues.resize(mConstitutiveLawVector.size());

        for ( IndexType GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint )
            rValues[GPoint] = mF0[GPoint];
    } else {
        UPwSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,
                                                                            rValues,
                                                                            rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                 std::vector<double>& rValues,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        KRATOS_ERROR_IF(rValues.size() != mConstitutiveLawVector.size()) 
            << "Can not set REFERENCE_DEFORMATION_GRADIENT_DETERMINANT, expected size: "
            << mConstitutiveLawVector.size() << " current size: " << rValues.size() << std::endl;

        for ( IndexType GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            mDetF0[GPoint] = rValues[GPoint];
        }
    } else {
        UPwSmallStrainElement<TDim,TNumNodes>::SetValuesOnIntegrationPoints(rVariable,
                                                                            rValues,
                                                                            rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                 std::vector<Matrix>& rValues,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
        KRATOS_ERROR_IF(rValues.size() != mConstitutiveLawVector.size())
            << "Can not set REFERENCE_DEFORMATION_GRADIENT, expected size: "
            << mConstitutiveLawVector.size() << " current size: " << rValues.size() << std::endl;

        for ( IndexType GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint )
            mF0[GPoint] = rValues[GPoint];
    } else {
        UPwSmallStrainElement<TDim,TNumNodes>::SetValuesOnIntegrationPoints(rVariable,
                                                                            rValues,
                                                                            rCurrentProcessInfo);
    }
}


//----------------------------------------------------------------------------------------

template class UPwUpdatedLagrangianElement<2,3>;
template class UPwUpdatedLagrangianElement<2,4>;
template class UPwUpdatedLagrangianElement<3,4>;
template class UPwUpdatedLagrangianElement<3,8>;

template class UPwUpdatedLagrangianElement<2,6>;
template class UPwUpdatedLagrangianElement<2,8>;
template class UPwUpdatedLagrangianElement<2,9>;
template class UPwUpdatedLagrangianElement<3,10>;
template class UPwUpdatedLagrangianElement<3,20>;
template class UPwUpdatedLagrangianElement<3,27>;

} // Namespace Kratos


