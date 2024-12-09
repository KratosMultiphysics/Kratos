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
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_utilities/constitutive_law_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/linear_nodal_extrapolator.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "includes/cfd_variables.h"

namespace
{

Kratos::BoundedMatrix<double, 4, 4> GetExtrapolationMatrixFor3D4NElement()
{
    auto result = Kratos::BoundedMatrix<double, 4, 4>{};

    result(0, 0) = -0.309016988749894905;
    result(0, 1) = -0.3090169887498949046;
    result(0, 2) = -0.309016988749894905;
    result(0, 3) = 1.9270509662496847144;

    result(1, 0) = 1.9270509662496847144;
    result(1, 1) = -0.30901698874989490481;
    result(1, 2) = -0.3090169887498949049;
    result(1, 3) = -0.30901698874989490481;

    result(2, 0) = -0.30901698874989490473;
    result(2, 1) = 1.9270509662496847143;
    result(2, 2) = -0.3090169887498949049;
    result(2, 3) = -0.30901698874989490481;

    result(3, 0) = -0.3090169887498949048;
    result(3, 1) = -0.30901698874989490471;
    result(3, 2) = 1.9270509662496847143;
    result(3, 3) = -0.30901698874989490481;

    return result;
}

Kratos::BoundedMatrix<double, 8, 8> GetExtrapolationMatrixFor3D8NElement()
{
    auto result = Kratos::BoundedMatrix<double, 8, 8>{};

    result(0, 0) = 2.549038105676658;
    result(0, 1) = -0.6830127018922192;
    result(0, 2) = 0.18301270189221927;
    result(0, 3) = -0.6830127018922192;
    result(0, 4) = -0.6830127018922192;
    result(0, 5) = 0.18301270189221927;
    result(0, 6) = -0.04903810567665795;
    result(0, 7) = 0.18301270189221927;

    result(1, 0) = -0.6830127018922192;
    result(1, 1) = 2.549038105676658;
    result(1, 2) = -0.6830127018922192;
    result(1, 3) = 0.18301270189221927;
    result(1, 4) = 0.18301270189221927;
    result(1, 5) = -0.6830127018922192;
    result(1, 6) = 0.18301270189221927;
    result(1, 7) = -0.04903810567665795;

    result(2, 0) = 0.18301270189221927;
    result(2, 1) = -0.6830127018922192;
    result(2, 2) = 2.549038105676658;
    result(2, 3) = -0.6830127018922192;
    result(2, 4) = -0.04903810567665795;
    result(2, 5) = 0.18301270189221927;
    result(2, 6) = -0.6830127018922192;
    result(2, 7) = 0.18301270189221927;

    result(3, 0) = -0.6830127018922192;
    result(3, 1) = 0.18301270189221927;
    result(3, 2) = -0.6830127018922192;
    result(3, 3) = 2.549038105676658;
    result(3, 4) = 0.18301270189221927;
    result(3, 5) = -0.04903810567665795;
    result(3, 6) = 0.18301270189221927;
    result(3, 7) = -0.6830127018922192;

    result(4, 0) = -0.6830127018922192;
    result(4, 1) = 0.18301270189221927;
    result(4, 2) = -0.04903810567665795;
    result(4, 3) = 0.18301270189221927;
    result(4, 4) = 2.549038105676658;
    result(4, 5) = -0.6830127018922192;
    result(4, 6) = 0.18301270189221927;
    result(4, 7) = -0.6830127018922192;

    result(5, 0) = 0.18301270189221927;
    result(5, 1) = -0.6830127018922192;
    result(5, 2) = 0.18301270189221927;
    result(5, 3) = -0.04903810567665795;
    result(5, 4) = -0.6830127018922192;
    result(5, 5) = 2.549038105676658;
    result(5, 6) = -0.6830127018922192;
    result(5, 7) = 0.18301270189221927;

    result(6, 0) = -0.04903810567665795;
    result(6, 1) = 0.18301270189221927;
    result(6, 2) = -0.6830127018922192;
    result(6, 3) = 0.18301270189221927;
    result(6, 4) = 0.18301270189221927;
    result(6, 5) = -0.6830127018922192;
    result(6, 6) = 2.549038105676658;
    result(6, 7) = -0.6830127018922192;

    result(7, 0) = 0.18301270189221927;
    result(7, 1) = -0.04903810567665795;
    result(7, 2) = 0.18301270189221927;
    result(7, 3) = -0.6830127018922192;
    result(7, 4) = -0.6830127018922192;
    result(7, 5) = 0.18301270189221927;
    result(7, 6) = -0.6830127018922192;
    result(7, 7) = 2.549038105676658;

    return result;
}

} // namespace

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                NodesArrayType const& ThisNodes,
                                                                PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainElement(NewId, this->GetGeometry().Create(ThisNodes),
                                                      pProperties, this->GetStressStatePolicy().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                GeometryType::Pointer pGeom,
                                                                PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
        new UPwSmallStrainElement(NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int UPwSmallStrainElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    // Verify generic variables
    int ierr = UPwBaseElement::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& rProp = this->GetProperties();
    const GeometryType&   rGeom = this->GetGeometry();

    KRATOS_ERROR_IF(rGeom.DomainSize() < 1.0e-15)
        << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    // Verify specific properties
    KRATOS_ERROR_IF_NOT(rProp.Has(IGNORE_UNDRAINED))
        << "IGNORE_UNDRAINED does not exist in the parameter list" << this->Id() << std::endl;

    if (!rProp[IGNORE_UNDRAINED]) {
        KRATOS_ERROR_IF(!rProp.Has(BULK_MODULUS_FLUID) || rProp[BULK_MODULUS_FLUID] < 0.0)
            << "BULK_MODULUS_FLUID has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        KRATOS_ERROR_IF(!rProp.Has(DYNAMIC_VISCOSITY) || rProp[DYNAMIC_VISCOSITY] < 0.0)
            << "DYNAMIC_VISCOSITY has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_XX) || rProp[PERMEABILITY_XX] < 0.0)
            << "PERMEABILITY_XX has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_YY) || rProp[PERMEABILITY_YY] < 0.0)
            << "PERMEABILITY_YY has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_XY) || rProp[PERMEABILITY_XY] < 0.0)
            << "PERMEABILITY_XY has Key zero, is not defined or has an invalid value at element"
            << this->Id() << std::endl;

        if constexpr (TDim > 2) {
            KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_ZZ) || rProp[PERMEABILITY_ZZ] < 0.0)
                << "PERMEABILITY_ZZ has Key zero, is not defined or has an invalid value at element"
                << this->Id() << std::endl;

            KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_YZ) || rProp[PERMEABILITY_YZ] < 0.0)
                << "PERMEABILITY_YZ has Key zero, is not defined or has an invalid value at element"
                << this->Id() << std::endl;

            KRATOS_ERROR_IF(!rProp.Has(PERMEABILITY_ZX) || rProp[PERMEABILITY_ZX] < 0.0)
                << "PERMEABILITY_ZX has Key zero, is not defined or has an invalid value at element"
                << this->Id() << std::endl;
        }
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has(CONSTITUTIVE_LAW))
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strainSize = this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    if (TDim == 2) {
        KRATOS_ERROR_IF_NOT(strainSize == VOIGT_SIZE_2D_PLANE_STRAIN)
            << "Wrong constitutive law used. This is a 2D element! expected "
               "strain size is "
            << VOIGT_SIZE_2D_PLANE_STRAIN << " But received: " << strainSize
            << " in element id: " << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strainSize == VOIGT_SIZE_3D)
            << "Wrong constitutive law used. This is a 3D element! expected "
               "strain size is "
            << VOIGT_SIZE_3D << " But received: " << strainSize << " in element id: " << this->Id()
            << std::endl;
    }

    // Check constitutive law
    if (mConstitutiveLawVector.size() > 0) {
        return mConstitutiveLawVector[0]->Check(rProp, rGeom, rCurrentProcessInfo);
    }

    // Check retention law
    if (mRetentionLawVector.size() > 0) {
        return mRetentionLawVector[0]->Check(rProp, rCurrentProcessInfo);
    }

    return ierr;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);

    ConstitutiveLaw::Parameters ConstitutiveParameters(this->GetGeometry(), this->GetProperties(),
                                                       rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); // Note: this is for nonlocal damage

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(deformation_gradients);
    const auto strain_vectors = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());

    const auto number_of_integration_points =
        this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);
        Variables.B            = b_matrices[GPoint];
        Variables.F            = deformation_gradients[GPoint];
        Variables.StrainVector = strain_vectors[GPoint];

        ConstitutiveLawUtilities::SetConstitutiveParameters(
            ConstitutiveParameters, Variables.StrainVector, Variables.ConstitutiveMatrix, Variables.Np,
            Variables.GradNpT, Variables.F, determinants_of_deformation_gradients[GPoint]);

        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->InitializeMaterialResponseCauchy(ConstitutiveParameters);
    }

    // Reset hydraulic discharge
    this->ResetHydraulicDischarge();

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::ResetHydraulicDischarge()
{
    KRATOS_TRY

    // Reset hydraulic discharge
    GeometryType& rGeom = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        ThreadSafeNodeWrite(rGeom[i], HYDRAULIC_DISCHARGE, 0.0);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateHydraulicDischarge(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    std::vector<array_1d<double, 3>> FluidFlux;
    this->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, FluidFlux, rCurrentProcessInfo);

    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(mThisIntegrationMethod);

    ElementVariables Variables;
    // Gradient of shape functions and determinant of Jacobian
    Variables.GradNpTInitialConfiguration.resize(TNumNodes, TDim, false);
    Variables.GradNpT.resize(TNumNodes, TDim, false);
    Variables.detJContainer.resize(NumGPoints, false);
    rGeom.ShapeFunctionsIntegrationPointsGradients(Variables.DN_DXContainer,
                                                   Variables.detJContainer, mThisIntegrationMethod);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        noalias(Variables.GradNpT) = Variables.DN_DXContainer[GPoint];
        Variables.detJ             = Variables.detJContainer[GPoint];

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        for (unsigned int node = 0; node < TNumNodes; ++node) {
            double HydraulicDischarge = 0;
            for (unsigned int iDir = 0; iDir < TDim; ++iDir) {
                HydraulicDischarge += Variables.GradNpT(node, iDir) * FluidFlux[GPoint][iDir];
            }

            HydraulicDischarge *= Variables.IntegrationCoefficient;
            HydraulicDischarge += rGeom[node].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE);
            ThreadSafeNodeWrite(this->GetGeometry()[node], HYDRAULIC_DISCHARGE, HydraulicDischarge);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); // Note: this is for nonlocal damage

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NContainer, Variables.DN_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateHydraulicDischarge(rCurrentProcessInfo);

    ConstitutiveLaw::Parameters ConstitutiveParameters(this->GetGeometry(), this->GetProperties(),
                                                       rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(deformation_gradients);
    const auto strain_vectors = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());

    const auto number_of_integration_points =
        this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    Matrix StressContainer(number_of_integration_points, mStressVector[0].size());

    for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);
        Variables.B            = b_matrices[GPoint];
        Variables.F            = deformation_gradients[GPoint];
        Variables.StrainVector = strain_vectors[GPoint];

        ConstitutiveLawUtilities::SetConstitutiveParameters(
            ConstitutiveParameters, Variables.StrainVector, Variables.ConstitutiveMatrix, Variables.Np,
            Variables.GradNpT, Variables.F, determinants_of_deformation_gradients[GPoint]);

        // Compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
        mStateVariablesFinalized[GPoint] =
            mConstitutiveLawVector[GPoint]->GetValue(STATE_VARIABLES, mStateVariablesFinalized[GPoint]);

        if (rCurrentProcessInfo[NODAL_SMOOTHING])
            this->SaveGPStress(StressContainer, mStressVector[GPoint], GPoint);
    }

    if (rCurrentProcessInfo[NODAL_SMOOTHING]) this->ExtrapolateGPValues(StressContainer);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::SaveGPStress(Matrix&       rStressContainer,
                                                          const Vector& StressVector,
                                                          unsigned int  GPoint)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < StressVector.size(); ++i) {
        rStressContainer(GPoint, i) = StressVector[i];
    }

    /* INFO: (Quadrilateral_2D_4 with GI_GAUSS_2)
     *
     *                      |S0-0 S1-0 S2-0|
     * rStressContainer =   |S0-1 S1-1 S2-1|
     *                      |S0-2 S1-2 S2-2|
     *                      |S0-3 S1-3 S2-3|
     *
     * S1-0 = S[1] at GP 0
     */

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::ExtrapolateGPValues(const Matrix& StressContainer)
{
    KRATOS_TRY

    array_1d<double, TNumNodes> DamageContainer;

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) {
        DamageContainer[iNode] = 0.0;
        DamageContainer[iNode] =
            mConstitutiveLawVector[iNode]->GetValue(DAMAGE_VARIABLE, DamageContainer[iNode]);
    }

    GeometryType&               rGeom = this->GetGeometry();
    const double&               Area  = rGeom.Area(); // In 3D this is volume
    array_1d<Vector, TNumNodes> NodalStressVector;    // List with stresses at each node
    array_1d<Matrix, TNumNodes> NodalStressTensor;

    auto const StressTensorSize = this->GetStressStatePolicy().GetStressTensorSize();
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) {
        NodalStressVector[iNode].resize(this->GetStressStatePolicy().GetVoigtSize());
        NodalStressTensor[iNode].resize(StressTensorSize, StressTensorSize);
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    Matrix AuxNodalStress;
    AuxNodalStress.resize(TNumNodes, this->GetStressStatePolicy().GetVoigtSize());
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix, StressContainer);

    /* INFO:
     *
     *                  |S0-0 S1-0 S2-0|
     * AuxNodalStress = |S0-1 S1-1 S2-1|
     *                  |S0-2 S1-2 S2-2|
     *
     * S1-0 = S[1] at node 0
     */

    array_1d<double, TNumNodes> NodalDamage;
    noalias(NodalDamage) = prod(ExtrapolationMatrix, DamageContainer);

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        noalias(NodalStressVector[i]) = row(AuxNodalStress, i) * Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) += NodalDamage[i] * Area;
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                          const std::vector<Vector>& rValues,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == CAUCHY_STRESS_VECTOR) {
        KRATOS_ERROR_IF(rValues.size() != mStressVector.size())
            << "Unexpected number of values for "
               "UPwSmallStrainElement::SetValuesOnIntegrationPoints"
            << std::endl;
        std::copy(rValues.begin(), rValues.end(), mStressVector.begin());
    } else {
        KRATOS_ERROR_IF(rValues.size() < mConstitutiveLawVector.size())
            << "Insufficient number of values for "
               "UPwSmallStrainElement::SetValuesOnIntegrationPoints"
            << std::endl;
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            mConstitutiveLawVector[GPoint]->SetValue(rVariable, rValues[GPoint], rCurrentProcessInfo);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                          std::vector<double>& rOutput,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const auto number_of_integration_points = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

    auto& r_prop = this->GetProperties();
    rOutput.resize(number_of_integration_points);

    if (rVariable == VON_MISES_STRESS) {
        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStress(mStressVector[GPoint]);
        }
    } else if (rVariable == MEAN_EFFECTIVE_STRESS) {
        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateMeanStress(mStressVector[GPoint]);
        }
    } else if (rVariable == MEAN_STRESS) {
        std::vector<Vector> StressVector;
        CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateMeanStress(StressVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_VON_MISES_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStrain(StrainVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_VOLUMETRIC_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateTrace(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VON_MISES_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateVonMisesStrain(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VOLUMETRIC_STRAIN) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = StressStrainUtilities::CalculateTrace(StrainVector[GPoint]);
        }
    } else if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION ||
               rVariable == BISHOP_COEFFICIENT || rVariable == DERIVATIVE_OF_SATURATION ||
               rVariable == RELATIVE_PERMEABILITY) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        RetentionLaw::Parameters RetentionParameters(this->GetProperties());

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            // Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            RetentionParameters.SetFluidPressure(GeoTransportEquationUtilities::CalculateFluidPressure(
                Variables.Np, Variables.PressureVector));

            if (rVariable == DEGREE_OF_SATURATION)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateSaturation(RetentionParameters);
            else if (rVariable == EFFECTIVE_SATURATION)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateEffectiveSaturation(RetentionParameters);
            else if (rVariable == BISHOP_COEFFICIENT)
                rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(RetentionParameters);
            else if (rVariable == DERIVATIVE_OF_SATURATION)
                rOutput[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(RetentionParameters);
            else if (rVariable == RELATIVE_PERMEABILITY)
                rOutput[GPoint] =
                    mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);
        }
    } else if (rVariable == HYDRAULIC_HEAD) {
        const PropertiesType& rProp = this->GetProperties();

        // Defining the shape functions, the Jacobian and the shape functions local gradients containers
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);

        const auto NodalHydraulicHead =
            GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(rGeom, rProp);

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            double HydraulicHead = 0.0;
            for (unsigned int node = 0; node < TNumNodes; ++node)
                HydraulicHead += NContainer(GPoint, node) * NodalHydraulicHead[node];

            rOutput[GPoint] = HydraulicHead;
        }
    } else if (rVariable == CONFINED_STIFFNESS || rVariable == SHEAR_STIFFNESS) {
        size_t variable_index = 0;
        if (rVariable == CONFINED_STIFFNESS) {
            if (TDim == 2) {
                variable_index = INDEX_2D_PLANE_STRAIN_XX;
            } else if (TDim == 3) {
                variable_index = INDEX_3D_XX;
            } else {
                KRATOS_ERROR << "CONFINED_STIFFNESS can not be retrieved for dim " << TDim
                             << " in element: " << this->Id() << std::endl;
            }
        } else if (rVariable == SHEAR_STIFFNESS) {
            if (TDim == 2) {
                variable_index = INDEX_2D_PLANE_STRAIN_XY;
            } else if (TDim == 3) {
                variable_index = INDEX_3D_XZ;
            } else {
                KRATOS_ERROR << "SHEAR_STIFFNESS can not be retrieved for dim " << TDim
                             << " in element: " << this->Id() << std::endl;
            }
        }

        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, this->GetStressStatePolicy().GetVoigtSize());

        const PropertiesType& rProp = this->GetProperties();

        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        std::vector<Matrix> constitutive_matrices;
        this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                             Variables.NContainer, Variables.DN_DXContainer,
                                             strain_vectors, mStressVector, constitutive_matrices);

        std::transform(constitutive_matrices.begin(), constitutive_matrices.end(), rOutput.begin(),
                       [variable_index](const Matrix& constitutive_matrix) {
            return constitutive_matrix(variable_index, variable_index);
        });
    } else if (r_prop.Has(rVariable)) {
        // Map initial material property to gauss points, as required for the output
        rOutput.clear();
        std::fill_n(std::back_inserter(rOutput), number_of_integration_points, r_prop.GetValue(rVariable));
    } else {
        for (unsigned int i = 0; i < number_of_integration_points; ++i) {
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom               = this->GetGeometry();
    const auto number_of_integration_points = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);
    rOutput.resize(number_of_integration_points);

    if (rVariable == FLUID_FLUX_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        const auto strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, this->GetStressStatePolicy().GetVoigtSize());

        const auto fluid_fluxes =
            CalculateFluidFluxes(GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(
                                     strain_vectors, this->GetProperties()),
                                 rCurrentProcessInfo);
        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], fluid_fluxes[GPoint]);
        }
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            rOutput[i] = ZeroVector(3);
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                          std::vector<Vector>& rOutput,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = this->GetGeometry();
    const auto number_of_integration_points = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);
    const auto& r_properties = this->GetProperties();

    rOutput.resize(number_of_integration_points);

    if (rVariable == CAUCHY_STRESS_VECTOR) {
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size() != mStressVector[GPoint].size())
                rOutput[GPoint].resize(mStressVector[GPoint].size(), false);

            rOutput[GPoint] = mStressVector[GPoint];
        }
    } else if (rVariable == TOTAL_STRESS_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        ConstitutiveLaw::Parameters ConstitutiveParameters(r_geometry, this->GetProperties(), rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, this->GetStressStatePolicy().GetVoigtSize());
        std::vector<Matrix> constitutive_matrices;
        this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                             Variables.NContainer, Variables.DN_DXContainer,
                                             strain_vectors, mStressVector, constitutive_matrices);
        const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
            constitutive_matrices, this->GetProperties());
        const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
            Variables.NContainer, Variables.PressureVector);
        const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] =
                mStressVector[GPoint] + PORE_PRESSURE_SIGN_FACTOR * biot_coefficients[GPoint] *
                                            bishop_coefficients[GPoint] * fluid_pressures[GPoint] *
                                            this->GetStressStatePolicy().GetVoigtVector();
        }
    } else if (rVariable == ENGINEERING_STRAIN_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            noalias(Variables.Np) = row(Variables.NContainer, GPoint);

            Matrix J0;
            Matrix InvJ0;
            this->CalculateDerivativesOnInitialConfiguration(
                Variables.detJInitialConfiguration, J0, InvJ0, Variables.GradNpTInitialConfiguration, GPoint);

            // Calculating operator B
            Variables.B = this->CalculateBMatrix(Variables.GradNpTInitialConfiguration, Variables.Np);

            // Compute infinitesimal strain
            Variables.StrainVector =
                StressStrainUtilities::CalculateCauchyStrain(Variables.B, Variables.DisplacementVector);

            if (rOutput[GPoint].size() != Variables.StrainVector.size())
                rOutput[GPoint].resize(Variables.StrainVector.size(), false);

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        rOutput                          = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, this->GetStressStatePolicy().GetVoigtSize());
    } else if (r_properties.Has(rVariable)) {
        // Map initial material property to Gauss points, as required for the output
        rOutput.clear();
        std::fill_n(std::back_inserter(rOutput), mConstitutiveLawVector.size(), r_properties.GetValue(rVariable));
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                          std::vector<Matrix>& rOutput,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = this->GetGeometry();
    const auto number_of_integration_points = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

    rOutput.resize(number_of_integration_points);

    if (rVariable == CAUCHY_STRESS_TENSOR) {
        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(mStressVector[GPoint]);
        }
    } else if (rVariable == TOTAL_STRESS_TENSOR) {
        std::vector<Vector> StressVector;
        this->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
        }
    } else if (rVariable == ENGINEERING_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else if (rVariable == PERMEABILITY_MATRIX) {
        // If the permeability of the element is a given property
        BoundedMatrix<double, TDim, TDim> PermeabilityMatrix;
        GeoElementUtilities::FillPermeabilityMatrix(PermeabilityMatrix, this->GetProperties());

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            rOutput[GPoint] = PermeabilityMatrix;
        }
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            rOutput[i] = ZeroMatrix(TDim, TDim);
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix,
                                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const IndexType N_DOF = TNumNodes * (TDim + 1);

    if (rStiffnessMatrix.size1() != N_DOF) rStiffnessMatrix.resize(N_DOF, N_DOF, false);
    noalias(rStiffnessMatrix) = ZeroMatrix(N_DOF, N_DOF);

    // Previous definitions
    const PropertiesType&                           rProp = this->GetProperties();
    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(mThisIntegrationMethod);

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);
    const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NContainer, Variables.DN_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);

    const auto stiffness_matrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        b_matrices, constitutive_matrices, integration_coefficients);

    GeoElementUtilities::AssembleUUBlockMatrix(rStiffnessMatrix, stiffness_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddGeometricStiffnessMatrix(
    MatrixType& rLeftHandSideMatrix, const Vector& rStressVector, const Matrix& rGradNpt, double IntegrationCoefficient) const
{
    KRATOS_TRY

    Matrix reduced_Kg_matrix =
        prod(rGradNpt, Matrix(prod(MathUtils<double>::StressVectorToTensor(rStressVector), trans(rGradNpt)))) *
        IntegrationCoefficient;

    Matrix geometric_stiffness_matrix = ZeroMatrix(TNumNodes * TDim, TNumNodes * TDim);
    MathUtils<double>::ExpandAndAddReducedMatrix(geometric_stiffness_matrix, reduced_Kg_matrix, TDim);

    // Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, geometric_stiffness_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geom             = this->GetGeometry();
    const auto          integration_method = this->GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geom.IntegrationPoints(integration_method);
    const auto N_container = r_geom.ShapeFunctionsValues(integration_method);

    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        N_container, this->GetPressureSolutionVector());
    const auto degrees_saturation = this->CalculateDegreesOfSaturation(fluid_pressures);

    const auto solid_densities =
        GeoTransportEquationUtilities::CalculateSoilDensities(degrees_saturation, this->GetProperties());

    const auto det_Js_initial_configuration =
        GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(r_geom, integration_method);

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(integration_points, det_Js_initial_configuration);

    const auto mass_matrix_u = GeoEquationOfMotionUtilities::CalculateMassMatrix(
        r_geom.WorkingSpaceDimension(), r_geom.PointsNumber(), integration_points.size(),
        r_geom.ShapeFunctionsValues(integration_method), solid_densities, integration_coefficients);

    rMassMatrix = ZeroMatrix(this->GetNumberOfDOF(), this->GetNumberOfDOF());
    GeoElementUtilities::AssembleUUBlockMatrix(rMassMatrix, mass_matrix_u);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                          VectorType&        rRightHandSideVector,
                                                          const ProcessInfo& rCurrentProcessInfo,
                                                          bool CalculateStiffnessMatrixFlag,
                                                          bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const PropertiesType&                           rProp = this->GetProperties();
    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);

    // Stiffness matrix is needed to calculate Biot coefficient
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    RetentionLaw::Parameters RetentionParameters(rProp);

    const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);
    const auto det_Js_initial_configuration =
        GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(rGeom, this->GetIntegrationMethod());
    const auto integration_coefficients_on_initial_configuration =
        this->CalculateIntegrationCoefficients(IntegrationPoints, det_Js_initial_configuration);

    const auto deformation_gradients = CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NContainer, Variables.DN_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);
    const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
        constitutive_matrices, this->GetProperties());
    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Variables.NContainer, Variables.PressureVector);
    const auto degrees_of_saturation     = CalculateDegreesOfSaturation(fluid_pressures);
    const auto derivatives_of_saturation = CalculateDerivativesOfSaturation(fluid_pressures);
    const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
        biot_coefficients, degrees_of_saturation, derivatives_of_saturation, rProp);
    auto relative_permeability_values = this->CalculateRelativePermeabilityValues(fluid_pressures);
    const auto permeability_update_factors = GetOptionalPermeabilityUpdateFactors(strain_vectors);
    std::transform(permeability_update_factors.cbegin(), permeability_update_factors.cend(),
                   relative_permeability_values.cbegin(), relative_permeability_values.begin(),
                   std::multiplies<>{});

    const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);

    for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);
        Variables.B                  = b_matrices[GPoint];
        Variables.F                  = deformation_gradients[GPoint];
        Variables.StrainVector       = strain_vectors[GPoint];
        Variables.ConstitutiveMatrix = constitutive_matrices[GPoint];

        // Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        Variables.RelativePermeability = relative_permeability_values[GPoint];
        Variables.BishopCoefficient    = bishop_coefficients[GPoint];

        Variables.BiotCoefficient        = biot_coefficients[GPoint];
        Variables.BiotModulusInverse     = biot_moduli_inverse[GPoint];
        Variables.DegreeOfSaturation     = degrees_of_saturation[GPoint];
        Variables.IntegrationCoefficient = integration_coefficients[GPoint];
        Variables.IntegrationCoefficientInitialConfiguration =
            integration_coefficients_on_initial_configuration[GPoint];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<double> UPwSmallStrainElement<TDim, TNumNodes>::GetOptionalPermeabilityUpdateFactors(
    const std::vector<Vector>& rStrainVectors) const
{
    return GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(rStrainVectors,
                                                                             this->GetProperties());
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<double> UPwSmallStrainElement<TDim, TNumNodes>::CalculateDerivativesOfSaturation(
    const std::vector<double>& rFluidPressures) const
{
    KRATOS_ERROR_IF(rFluidPressures.size() != mRetentionLawVector.size());
    std::vector<double> result;

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};
    std::transform(rFluidPressures.begin(), rFluidPressures.end(), mRetentionLawVector.begin(),
                   std::back_inserter(result), [&retention_law_params](auto fluid_pressure, auto pRetentionLaw) {
        retention_law_params.SetFluidPressure(fluid_pressure);
        return pRetentionLaw->CalculateDerivativeOfSaturation(retention_law_params);
    });

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<double> UPwSmallStrainElement<TDim, TNumNodes>::CalculateDegreesOfSaturation(
    const std::vector<double>& rFluidPressures) const
{
    KRATOS_ERROR_IF(rFluidPressures.size() != mRetentionLawVector.size());
    std::vector<double> result;

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};
    std::transform(rFluidPressures.begin(), rFluidPressures.end(), mRetentionLawVector.begin(),
                   std::back_inserter(result), [&retention_law_params](auto fluid_pressure, auto pRetentionLaw) {
        retention_law_params.SetFluidPressure(fluid_pressure);
        return pRetentionLaw->CalculateSaturation(retention_law_params);
    });

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<array_1d<double, TDim>> UPwSmallStrainElement<TDim, TNumNodes>::CalculateFluidFluxes(
    const std::vector<double>& rPermeabilityUpdateFactors, const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());

    std::vector<array_1d<double, TDim>> FluidFluxes;
    ElementVariables                    Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    const PropertiesType& rProp = this->GetProperties();

    auto relative_permeability_values = this->CalculateRelativePermeabilityValues(
        GeoTransportEquationUtilities::CalculateFluidPressures(Variables.NContainer, Variables.PressureVector));
    std::transform(relative_permeability_values.cbegin(), relative_permeability_values.cend(),
                   rPermeabilityUpdateFactors.cbegin(), relative_permeability_values.begin(),
                   std::multiplies<>{});

    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);

        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        Variables.RelativePermeability = relative_permeability_values[GPoint];

        array_1d<double, TDim> GradPressureTerm = prod(trans(Variables.GradNpT), Variables.PressureVector);
        GradPressureTerm += PORE_PRESSURE_SIGN_FACTOR * rProp[DENSITY_WATER] * Variables.BodyAcceleration;

        FluidFluxes.push_back(PORE_PRESSURE_SIGN_FACTOR * Variables.DynamicViscosityInverse *
                              Variables.RelativePermeability *
                              prod(Variables.PermeabilityMatrix, GradPressureTerm));
    }

    return FluidFluxes;
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Properties variables
    this->InitializeProperties(rVariables);

    // ProcessInfo variables
    rVariables.VelocityCoefficient   = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    // Nodal Variables
    this->InitializeNodalDisplacementVariables(rVariables);
    this->InitializeNodalPorePressureVariables(rVariables);
    this->InitializeNodalVolumeAccelerationVariables(rVariables);

    // Variables computed at each GP
    rVariables.Nu = ZeroMatrix(TDim, TNumNodes * TDim);
    rVariables.Np.resize(TNumNodes, false);
    rVariables.GradNpT.resize(TNumNodes, TDim, false);
    rVariables.F = identity_matrix<double>(TDim);

    rVariables.B = ZeroMatrix(this->GetStressStatePolicy().GetVoigtSize(), TNumNodes * TDim);

    const GeometryType& rGeom      = this->GetGeometry();
    const IndexType     NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    // Shape functions
    rVariables.NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);

    // Gradient of shape functions and determinant of Jacobian
    rVariables.detJContainer.resize(NumGPoints, false);

    rGeom.ShapeFunctionsIntegrationPointsGradients(
        rVariables.DN_DXContainer, rVariables.detJContainer, mThisIntegrationMethod);

    // Constitutive Law parameters
    rVariables.StressVector.resize(this->GetStressStatePolicy().GetVoigtSize(), false);
    rVariables.StrainVector.resize(this->GetStressStatePolicy().GetVoigtSize(), false);
    rVariables.ConstitutiveMatrix.resize(this->GetStressStatePolicy().GetVoigtSize(),
                                         this->GetStressStatePolicy().GetVoigtSize(), false);

    // Auxiliary variables
    rVariables.UVoigtMatrix.resize(TNumNodes * TDim, this->GetStressStatePolicy().GetVoigtSize(), false);

    // Retention law
    rVariables.DegreeOfSaturation   = 1.0;
    rVariables.RelativePermeability = 1.0;
    rVariables.BishopCoefficient    = 1.0;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
Matrix UPwSmallStrainElement<TDim, TNumNodes>::CalculateBMatrix(const Matrix& rDN_DX, const Vector& rN) const
{
    return this->GetStressStatePolicy().CalculateBMatrix(rDN_DX, rN, this->GetGeometry());
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<Matrix> UPwSmallStrainElement<TDim, TNumNodes>::CalculateBMatrices(
    const GeometryType::ShapeFunctionsGradientsType& rDN_DXContainer, const Matrix& rNContainer) const
{
    std::vector<Matrix> result;
    for (unsigned int GPoint = 0; GPoint < rDN_DXContainer.size(); ++GPoint) {
        result.push_back(this->CalculateBMatrix(rDN_DXContainer[GPoint], row(rNContainer, GPoint)));
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);

        const auto permeability_matrix =
            GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
                rVariables.GradNpT, rVariables.DynamicViscosityInverse, rVariables.PermeabilityMatrix,
                rVariables.RelativePermeability, rVariables.IntegrationCoefficient);
        GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, permeability_matrix);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix,
                                                                            const ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto stiffness_matrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrixGPoint(
        rVariables.B, rVariables.ConstitutiveMatrix, rVariables.IntegrationCoefficient);
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, stiffness_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix,
                                                                           const ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix(
        rVariables.B, this->GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
        rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);
    GeoElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix, coupling_matrix);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        const BoundedMatrix<double, TNumNodes, TNumNodes * TDim> transposed_coupling_matrix =
            PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient * rVariables.VelocityCoefficient *
            trans(coupling_matrix);
        GeoElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix, transposed_coupling_matrix);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                                  const ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    // Distribute compressibility block matrix into the elemental matrix
    GeoElementUtilities::AssemblePPBlockMatrix(
        rLeftHandSideMatrix, compressibility_matrix * rVariables.DtPressureCoefficient);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                ElementVariables& rVariables,
                                                                unsigned int      GPoint)
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables, GPoint);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector,
                                                                           const ElementVariables& rVariables,
                                                                           unsigned int GPoint)
{
    KRATOS_TRY
    const array_1d<double, TNumNodes * TDim> stiffness_force =
        -1.0 * prod(trans(rVariables.B), mStressVector[GPoint]) * rVariables.IntegrationCoefficient;

    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, stiffness_force);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector,
                                                                         ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateSoilGamma(rVariables);

    const array_1d<double, TNumNodes * TDim> mix_body_force =
        prod(trans(rVariables.Nu), rVariables.SoilGamma) * rVariables.IntegrationCoefficientInitialConfiguration;

    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, mix_body_force);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateSoilGamma(ElementVariables& rVariables)
{
    KRATOS_TRY

    noalias(rVariables.SoilGamma) = GeoTransportEquationUtilities::CalculateSoilDensity(
                                        rVariables.DegreeOfSaturation, this->GetProperties()) *
                                    rVariables.BodyAcceleration;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector,
                                                                          const ElementVariables& rVariables)
{
    KRATOS_TRY

    const BoundedMatrix<double, TNumNodes * TDim, TNumNodes> coupling_matrix =
        (-1.0) * GeoTransportEquationUtilities::CalculateCouplingMatrix(
                     rVariables.B, this->GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
                     rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);

    const array_1d<double, TNumNodes * TDim> coupling_terms = prod(coupling_matrix, rVariables.PressureVector);
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, coupling_terms);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        const array_1d<double, TNumNodes> transposed_coupling_matrix =
            PORE_PRESSURE_SIGN_FACTOR * SaturationCoefficient *
            prod(trans(coupling_matrix), rVariables.VelocityVector);
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, transposed_coupling_matrix);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TNumNodes> UPwSmallStrainElement<TDim, TNumNodes>::CalculateCompressibilityFlow(const ElementVariables& rVariables) const
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);
    return -prod(compressibility_matrix, rVariables.DtPressureVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector,
                                                                                const ElementVariables& rVariables)
{
    KRATOS_TRY

    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector,
                                              this->CalculateCompressibilityFlow(rVariables));

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TNumNodes> UPwSmallStrainElement<TDim, TNumNodes>::CalculatePermeabilityFlow(const ElementVariables& rVariables) const
{
    KRATOS_TRY

    const auto permeability_matrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
        rVariables.GradNpT, rVariables.DynamicViscosityInverse, rVariables.PermeabilityMatrix,
        rVariables.RelativePermeability, rVariables.IntegrationCoefficient);

    return -prod(permeability_matrix, rVariables.PressureVector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<double> UPwSmallStrainElement<TDim, TNumNodes>::CalculateRelativePermeabilityValues(
    const std::vector<double>& rFluidPressures) const
{
    KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

    auto result = std::vector<double>{};
    std::transform(mRetentionLawVector.begin(), mRetentionLawVector.end(), rFluidPressures.begin(),
                   std::back_inserter(result), [&retention_law_params](auto pRetentionLaw, auto FluidPressure) {
        retention_law_params.SetFluidPressure(FluidPressure);
        return pRetentionLaw->CalculateRelativePermeability(retention_law_params);
    });
    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<double> UPwSmallStrainElement<TDim, TNumNodes>::CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const
{
    KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

    auto result = std::vector<double>{};
    std::transform(mRetentionLawVector.begin(), mRetentionLawVector.end(), rFluidPressures.begin(),
                   std::back_inserter(result), [&retention_law_params](auto pRetentionLaw, auto FluidPressure) {
        retention_law_params.SetFluidPressure(FluidPressure);
        return pRetentionLaw->CalculateBishopCoefficient(retention_law_params);
    });
    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                                                             const ElementVariables& rVariables)
{
    KRATOS_TRY

    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, this->CalculatePermeabilityFlow(rVariables));

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TNumNodes> UPwSmallStrainElement<TDim, TNumNodes>::CalculateFluidBodyFlow(const ElementVariables& rVariables) const
{
    KRATOS_TRY

    const BoundedMatrix<double, TNumNodes, TDim> temp_matrix =
        prod(rVariables.GradNpT, rVariables.PermeabilityMatrix) * rVariables.IntegrationCoefficient;

    return rVariables.DynamicViscosityInverse * this->GetProperties()[DENSITY_WATER] * rVariables.BishopCoefficient *
           rVariables.RelativePermeability * prod(temp_matrix, rVariables.BodyAcceleration);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                          const ElementVariables& rVariables)
{
    KRATOS_TRY

    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, this->CalculateFluidBodyFlow(rVariables));

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
Vector UPwSmallStrainElement<TDim, TNumNodes>::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
{
    return this->GetStressStatePolicy().CalculateGreenLagrangeStrain(rDeformationGradient);
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<Matrix> UPwSmallStrainElement<TDim, TNumNodes>::CalculateDeformationGradients() const
{
    std::vector<Matrix> result;
    for (unsigned int GPoint = 0;
         GPoint < this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()); ++GPoint) {
        result.push_back(CalculateDeformationGradient(GPoint));
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
Matrix UPwSmallStrainElement<TDim, TNumNodes>::CalculateDeformationGradient(unsigned int GPoint) const
{
    KRATOS_TRY

    // Calculation of derivative of shape function with respect to reference
    // configuration derivative of shape function (displacement)
    Matrix J0;
    Matrix InvJ0;
    Matrix DNu_DX0;
    double detJ0;
    this->CalculateDerivativesOnInitialConfiguration(detJ0, J0, InvJ0, DNu_DX0, GPoint);

    // Calculating current Jacobian in order to find deformation gradient
    Matrix J;
    Matrix InvJ;
    double detJ;
    this->CalculateJacobianOnCurrentConfiguration(detJ, J, InvJ, GPoint);

    KRATOS_ERROR_IF(detJ < 0.0)
        << "ERROR:: Element " << this->Id() << " is inverted. DetJ: " << detJ << std::endl
        << "This usually indicates that the deformations are too large for the mesh size." << std::endl;

    return prod(J, InvJ0);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNodalPorePressureVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    // Nodal variables
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rVariables.PressureVector[i]   = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.DtPressureVector[i] = rGeom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNodalDisplacementVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    // Nodal variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector, rGeom, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector, rGeom, VELOCITY);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNodalVolumeAccelerationVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    // Nodal variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration,
                                                                 rGeom, VOLUME_ACCELERATION);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeProperties(ElementVariables& rVariables)
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    rVariables.IgnoreUndrained = rProp[IGNORE_UNDRAINED];
    rVariables.UseHenckyStrain = false;
    if (rProp.Has(USE_HENCKY_STRAIN)) rVariables.UseHenckyStrain = rProp[USE_HENCKY_STRAIN];

    rVariables.ConsiderGeometricStiffness = false;
    if (rProp.Has(CONSIDER_GEOMETRIC_STIFFNESS))
        rVariables.ConsiderGeometricStiffness = rProp[CONSIDER_GEOMETRIC_STIFFNESS];

    rVariables.DynamicViscosityInverse = 1.0 / rProp[DYNAMIC_VISCOSITY];
    GeoElementUtilities::FillPermeabilityMatrix(rVariables.PermeabilityMatrix, rProp);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateKinematics(ElementVariables& rVariables, unsigned int GPoint)
{
    KRATOS_TRY

    // Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Np)      = row(rVariables.NContainer, GPoint);
    noalias(rVariables.GradNpT) = rVariables.DN_DXContainer[GPoint];

    rVariables.detJ = rVariables.detJContainer[GPoint];

    Matrix J0, InvJ0;
    this->CalculateDerivativesOnInitialConfiguration(rVariables.detJInitialConfiguration, J0, InvJ0,
                                                     rVariables.GradNpTInitialConfiguration, GPoint);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateExtrapolationMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rExtrapolationMatrix)
{
    KRATOS_TRY

    if constexpr ((TDim == 2u) && ((TNumNodes == 3u) || (TNumNodes == 4u))) {
        const auto extrapolator = LinearNodalExtrapolator{};
        const auto result       = extrapolator.CalculateElementExtrapolationMatrix(
            this->GetGeometry(), this->GetIntegrationMethod());
        KRATOS_ERROR_IF_NOT(result.size1() == TNumNodes)
            << "Extrapolation matrix has unexpected number of rows: " << result.size1()
            << " (expected " << TNumNodes << ")" << std::endl;
        KRATOS_ERROR_IF_NOT(result.size2() == TNumNodes)
            << "Extrapolation matrix has unexpected number of columns: " << result.size2()
            << " (expected " << TNumNodes << ")" << std::endl;
        rExtrapolationMatrix = result;
    } else if constexpr ((TDim == 3u) && (TNumNodes == 4u)) {
        rExtrapolationMatrix = ::GetExtrapolationMatrixFor3D4NElement();
    } else if constexpr ((TDim == 3u) && (TNumNodes == 8u)) {
        rExtrapolationMatrix = ::GetExtrapolationMatrixFor3D8NElement();
    } else {
        KRATOS_ERROR << "undefined number of nodes in CalculateExtrapolationMatrix "
                        "... TNumNodes:"
                     << TNumNodes << " element: " << this->Id() << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
Vector UPwSmallStrainElement<TDim, TNumNodes>::GetPressureSolutionVector()
{
    Vector result(TNumNodes);
    std::transform(this->GetGeometry().begin(), this->GetGeometry().end(), result.begin(),
                   [](const auto& node) { return node.FastGetSolutionStepValue(WATER_PRESSURE); });
    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAnyOfMaterialResponse(
    const std::vector<Matrix>&                       rDeformationGradients,
    ConstitutiveLaw::Parameters&                     rConstitutiveParameters,
    const Matrix&                                    rNuContainer,
    const GeometryType::ShapeFunctionsGradientsType& rDNu_DXContainer,
    std::vector<Vector>&                             rStrainVectors,
    std::vector<Vector>&                             rStressVectors,
    std::vector<Matrix>&                             rConstitutiveMatrices)
{
    if (rStrainVectors.size() != rDeformationGradients.size()) {
        rStrainVectors.resize(rDeformationGradients.size());
        std::fill(rStrainVectors.begin(), rStrainVectors.end(),
                  ZeroVector(this->GetStressStatePolicy().GetVoigtSize()));
    }
    if (rStressVectors.size() != rDeformationGradients.size()) {
        rStressVectors.resize(rDeformationGradients.size());
        std::fill(rStressVectors.begin(), rStressVectors.end(),
                  ZeroVector(this->GetStressStatePolicy().GetVoigtSize()));
    }
    if (rConstitutiveMatrices.size() != rDeformationGradients.size()) {
        rConstitutiveMatrices.resize(rDeformationGradients.size());
        std::fill(rConstitutiveMatrices.begin(), rConstitutiveMatrices.end(),
                  ZeroMatrix(this->GetStressStatePolicy().GetVoigtSize(),
                             this->GetStressStatePolicy().GetVoigtSize()));
    }

    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(rDeformationGradients);

    for (unsigned int GPoint = 0; GPoint < rDeformationGradients.size(); ++GPoint) {
        ConstitutiveLawUtilities::SetConstitutiveParameters(
            rConstitutiveParameters, rStrainVectors[GPoint], rConstitutiveMatrices[GPoint],
            row(rNuContainer, GPoint), rDNu_DXContainer[GPoint], rDeformationGradients[GPoint],
            determinants_of_deformation_gradients[GPoint]);
        rConstitutiveParameters.SetStressVector(rStressVectors[GPoint]);

        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(rConstitutiveParameters);
    }
}

template class UPwSmallStrainElement<2, 3>;
template class UPwSmallStrainElement<2, 4>;
template class UPwSmallStrainElement<3, 4>;
template class UPwSmallStrainElement<3, 8>;

template class UPwSmallStrainElement<2, 6>;
template class UPwSmallStrainElement<2, 8>;
template class UPwSmallStrainElement<2, 9>;
template class UPwSmallStrainElement<2, 10>;
template class UPwSmallStrainElement<2, 15>;
template class UPwSmallStrainElement<3, 10>;
template class UPwSmallStrainElement<3, 20>;
template class UPwSmallStrainElement<3, 27>;

} // Namespace Kratos
