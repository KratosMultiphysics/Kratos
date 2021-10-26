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
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& rCurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag )
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType &IntegrationPoints =
        this->GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(this->GetGeometry(),
                                                       this->GetProperties(),
                                                       rCurrentProcessInfo);

    // Stiffness matrix is always needed for Biot coefficient
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)  ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    // Computing in all integrations points
    for ( IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint ) {
        // Compute element kinematics B, F, GradNpT ...
        this->CalculateKinematics(Variables, GPoint);

        //Compute Np, Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>( Variables.BodyAcceleration,
                                                                Variables.NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );

        // Cauchy strain: This needs to be investigated which strain measure should be used
        // In some references, e.g. Bathe, suggested to use Almansi strain measure
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        //Compute constitutive tensor and stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = this->CalculateBulkModulus(Variables.ConstitutiveMatrix);
        this->InitializeBiotCoefficients(Variables, BulkModulus);

        // Calculating weights for integration on the reference configuration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJ);

        if (CalculateStiffnessMatrixFlag) {
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

            /* Geometric stiffness matrix */
            if (Variables.ConsiderGeometricStiffness)
                this->CalculateAndAddGeometricStiffnessMatrix( rLeftHandSideMatrix, Variables, GPoint );
        }

        if (CalculateResidualVectorFlag) {
            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateAndAddGeometricStiffnessMatrix( MatrixType& rLeftHandSideMatrix,
                                             ElementVariables& rVariables,
                                             unsigned int GPoint )
{
    KRATOS_TRY

    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( mStressVector[GPoint] );
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

// //----------------------------------------------------------------------------------------
// template< unsigned int TDim, unsigned int TNumNodes >
// void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
//     CalculateKinematics(ElementVariables &rVariables,
//                         const unsigned int &GPoint)
// {
//     noalias(rVariables.Np) = row(rVariables.NContainer, GPoint);

//     Matrix J0, InvJ0;
//     rVariables.detJInitialConfiguration =
//         this->CalculateDerivativesOnReferenceConfiguration(J0,
//                                                            InvJ0,
//                                                            rVariables.GradNpT,
//                                                            GPoint,
//                                                            mThisIntegrationMethod);

//     // Calculating operator B
//     this->CalculateBMatrix( rVariables.B, rVariables.GradNpT, rVariables.Np);

//     // Calculating jacobian
//     Matrix J, InvJ;
//     rVariables.detJ =
//         this->CalculateDerivativesOnCurrentConfiguration(J,
//                                                          InvJ,
//                                                          rVariables.GradNpT,
//                                                          GPoint,
//                                                          mThisIntegrationMethod);

// #ifdef KRATOS_COMPILED_IN_WINDOWS
//     if (rVariables.detJ < 0.0)
//     {
//         KRATOS_INFO("negative detJ")
//         << "ERROR:: ELEMENT ID: "
//         << this->Id()
//         << " INVERTED. DETJ: "
//         << rVariables.detJ
//         << " nodes:" << this->GetGeometry()
//         << std::endl;
//     }
// #endif

//     KRATOS_ERROR_IF(rVariables.detJ < 0.0)
//      << "ERROR:: ELEMENT ID: "
//      << this->Id()
//      << " INVERTED. DETJ: "
//      << rVariables.detJ
//      << " nodes:" << this->GetGeometry()
//      << std::endl;

//     // Deformation gradient
//     // Matrix DF = prod( J, rVariables.InvJ0 );
//     // const double detDF = MathUtils<double>::Det(DF);
//     // rVariables.detF = detDF * this->ReferenceConfigurationDeformationGradientDeterminant(GPoint);
//     // noalias(rVariables.F) = prod(DF, this->ReferenceConfigurationDeformationGradient(GPoint));

//     noalias(rVariables.F) = prod( J, InvJ0 );
//     rVariables.detF = MathUtils<double>::Det(rVariables.F);

// }


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
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                 std::vector<double>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        ElementVariables Variables;
        this->InitializeElementVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            this->CalculateDeformationGradient(Variables, GPoint);
            rOutput[GPoint] = Variables.detF;
        }

    } else {
        UPwSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,
                                                                            rOutput,
                                                                            rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwUpdatedLagrangianElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                 std::vector<Matrix>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        ElementVariables Variables;
        this->InitializeElementVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            this->CalculateDeformationGradient(Variables, GPoint);
            rOutput[GPoint] = Variables.F;
        }
    } else {
        UPwSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,
                                                                            rOutput,
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


