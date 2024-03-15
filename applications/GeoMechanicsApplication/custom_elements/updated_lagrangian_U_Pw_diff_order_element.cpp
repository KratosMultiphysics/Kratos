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

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_Pw_diff_order_element.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
Element::Pointer UpdatedLagrangianUPwDiffOrderElement::Create(IndexType             NewId,
                                                              NodesArrayType const& ThisNodes,
                                                              PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UpdatedLagrangianUPwDiffOrderElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
Element::Pointer UpdatedLagrangianUPwDiffOrderElement::Create(IndexType             NewId,
                                                              GeometryType::Pointer pGeom,
                                                              PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UpdatedLagrangianUPwDiffOrderElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                        VectorType&        rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo,
                                                        const bool CalculateStiffnessMatrixFlag,
                                                        const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const GeometryType&   rGeom = GetGeometry();
    const PropertiesType& rProp = this->GetProperties();

    // Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Stiffness matrix is always needed to calculate Biot coefficient
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(rProp, rCurrentProcessInfo);

    // Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

    // Computing in all integrations points
    for (IndexType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        // Compute element kinematics B, F, DNu_DX ...
        this->CalculateKinematics(Variables, GPoint);

        // Compute strain
        this->CalculateStrain(Variables, GPoint);

        // set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute constitutive tensor and stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // calculate Bulk modulus from stiffness matrix
        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        // Calculating weights for integration on the reference configuration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

        Variables.IntegrationCoefficientInitialConfiguration = this->CalculateIntegrationCoefficient(
            IntegrationPoints, GPoint, Variables.detJInitialConfiguration);

        if (CalculateStiffnessMatrixFlag) {
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

            /* Geometric stiffness matrix */
            if (Variables.ConsiderGeometricStiffness)
                this->CalculateAndAddGeometricStiffnessMatrix(rLeftHandSideMatrix, Variables, GPoint);
        }

        if (CalculateResidualVectorFlag) {
            // Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateAndAddGeometricStiffnessMatrix(MatrixType& rLeftHandSideMatrix,
                                                                                   ElementVariables& rVariables,
                                                                                   unsigned int GPoint)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();

    Matrix StressTensor = MathUtils<double>::StressVectorToTensor(mStressVector[GPoint]);

    Matrix ReducedKgMatrix =
        prod(rVariables.DNu_DX,
             rVariables.IntegrationCoefficient * Matrix(prod(StressTensor, trans(rVariables.DNu_DX)))); // to be optimized

    Matrix UUMatrix(NumUNodes * Dim, NumUNodes * Dim);
    noalias(UUMatrix) = ZeroMatrix(NumUNodes * Dim, NumUNodes * Dim);
    MathUtils<double>::ExpandAndAddReducedMatrix(UUMatrix, ReducedKgMatrix, Dim);

    // Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, UUMatrix);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                        std::vector<double>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            this->CalculateDeformationGradient(Variables, GPoint);
            rOutput[GPoint] = Variables.detF;
        }
    } else {
        SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                        std::vector<Vector>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int& IntegrationPointsNumber =
        rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

    if (rOutput.size() != IntegrationPointsNumber) rOutput.resize(IntegrationPointsNumber);

    if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        // Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            this->CalculateDeformationGradient(Variables, GPoint);
            Variables.StrainVector = this->CalculateGreenLagrangeStrain(Variables.F);

            if (rOutput[GPoint].size() != Variables.StrainVector.size())
                rOutput[GPoint].resize(Variables.StrainVector.size(), false);

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else {
        SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                        std::vector<Matrix>& rOutput,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType      Dim   = rGeom.WorkingSpaceDimension();

    if (rOutput.size() != mConstitutiveLawVector.size())
        rOutput.resize(mConstitutiveLawVector.size());

    if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            this->CalculateDeformationGradient(Variables, GPoint);

            if (rOutput[GPoint].size2() != Dim) rOutput[GPoint].resize(Dim, Dim, false);

            rOutput[GPoint] = Variables.F;
        }
    } else {
        SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainUPwDiffOrderElement)
}

//----------------------------------------------------------------------------------------
void UpdatedLagrangianUPwDiffOrderElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainUPwDiffOrderElement)
}

} // Namespace Kratos
