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
#include "custom_elements/U_Pw_small_strain_element.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/hydraulic_discharge.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/node_utilities.h"
#include "custom_utilities/output_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"

#include <numeric>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
UPwSmallStrainElement<TDim, TNumNodes>::UPwSmallStrainElement(IndexType             NewId,
                                                              const NodesArrayType& ThisNodes,
                                                              std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                                              std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : UPwBaseElement(NewId, ThisNodes, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
{
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwSmallStrainElement<TDim, TNumNodes>::UPwSmallStrainElement(IndexType             NewId,
                                                              GeometryType::Pointer pGeometry,
                                                              std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                                              std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : UPwBaseElement(NewId, pGeometry, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
{
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwSmallStrainElement<TDim, TNumNodes>::UPwSmallStrainElement(IndexType               NewId,
                                                              GeometryType::Pointer   pGeometry,
                                                              PropertiesType::Pointer pProperties,
                                                              std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                                              std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : UPwBaseElement(NewId, pGeometry, pProperties, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                NodesArrayType const& ThisNodes,
                                                                PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainElement(NewId, this->GetGeometry().Create(ThisNodes),
                                                      pProperties, this->GetStressStatePolicy().Clone(),
                                                      this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                GeometryType::Pointer pGeom,
                                                                PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainElement(NewId, pGeom, pProperties,
                                                      this->GetStressStatePolicy().Clone(),
                                                      this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int UPwSmallStrainElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    // Verify generic variables
    if (const auto ierr = UPwBaseElement::Check(rCurrentProcessInfo); ierr != 0) return ierr;

    const auto& r_properties = this->GetProperties();
    const auto& r_geometry   = this->GetGeometry();

    CheckUtilities::CheckDomainSize(r_geometry.DomainSize(), this->Id());

    const CheckProperties check_properties(r_properties, "property", this->Id(),
                                           CheckProperties::Bounds::AllInclusive);
    check_properties.CheckAvailability(IGNORE_UNDRAINED);
    if (!r_properties[IGNORE_UNDRAINED]) {
        check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(BULK_MODULUS_FLUID);
        check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(DYNAMIC_VISCOSITY);
        check_properties.CheckPermeabilityProperties(r_geometry.WorkingSpaceDimension());
    }

    check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW);
    const auto expected_size = this->GetStressStatePolicy().GetVoigtSize();
    ConstitutiveLawUtilities::CheckStrainSize(r_properties, expected_size, this->Id());

    const auto error_code = r_properties[CONSTITUTIVE_LAW]->Check(r_properties, r_geometry, rCurrentProcessInfo);

    return error_code + RetentionLaw::Check(mRetentionLawVector, r_properties, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    this->ResetHydraulicDischarge();
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::ResetHydraulicDischarge()
{
    KRATOS_TRY

    // Reset hydraulic discharge
    for (auto& r_node : this->GetGeometry()) {
        NodeUtilities::ThreadSafeNodeWrite(r_node, HYDRAULIC_DISCHARGE, 0.0);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateHydraulicDischarge(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    std::vector<array_1d<double, 3>> fluid_flux;
    this->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, fluid_flux, rCurrentProcessInfo);

    GeometryType& r_geometry = this->GetGeometry();
    const IndexType number_of_integration_points = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(mThisIntegrationMethod);

    ElementVariables Variables;
    // Gradient of shape functions and determinant of Jacobian
    Variables.GradNpTInitialConfiguration.resize(TNumNodes, TDim, false);
    Variables.GradNpT.resize(TNumNodes, TDim, false);
    Variables.detJContainer.resize(number_of_integration_points, false);
    r_geometry.ShapeFunctionsIntegrationPointsGradients(
        Variables.DN_DXContainer, Variables.detJContainer, mThisIntegrationMethod);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);

    HydraulicDischarge::CalculateHydraulicDischarge(
        fluid_flux, integration_coefficients, Variables.DN_DXContainer, mThisIntegrationMethod, r_geometry);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConstitutiveLaw::Parameters ConstitutiveParameters(this->GetGeometry(), this->GetProperties(),
                                                       rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE);

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

    const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(deformation_gradients);
    const auto strain_vectors = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());

    const auto number_of_integration_points =
        this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

    for (unsigned int integration_point = 0; integration_point < number_of_integration_points; ++integration_point) {
        this->CalculateKinematics(Variables, integration_point);
        Variables.B            = b_matrices[integration_point];
        Variables.F            = deformation_gradients[integration_point];
        Variables.StrainVector = strain_vectors[integration_point];

        ConstitutiveLawUtilities::SetConstitutiveParameters(
            ConstitutiveParameters, Variables.StrainVector, Variables.ConstitutiveMatrix, Variables.Np,
            Variables.GradNpT, Variables.F, determinants_of_deformation_gradients[integration_point]);

        // Compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[integration_point];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[integration_point]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
        mStateVariablesFinalized[integration_point] = mConstitutiveLawVector[integration_point]->GetValue(
            STATE_VARIABLES, mStateVariablesFinalized[integration_point]);
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
        KRATOS_ERROR_IF(rValues.size() != this->GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod))
            << "Unexpected number of values for "
               "UPwSmallStrainElement::SetValuesOnIntegrationPoints"
            << std::endl;
        mStressVector.resize(rValues.size());
        std::copy(rValues.begin(), rValues.end(), mStressVector.begin());
    } else {
        KRATOS_ERROR_IF(rValues.size() < mConstitutiveLawVector.size())
            << "Insufficient number of values for "
               "UPwSmallStrainElement::SetValuesOnIntegrationPoints"
            << std::endl;
        for (unsigned int integration_point = 0; integration_point < mConstitutiveLawVector.size();
             ++integration_point) {
            mConstitutiveLawVector[integration_point]->SetValue(
                rVariable, rValues[integration_point], rCurrentProcessInfo);
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

    const GeometryType& r_geometry = this->GetGeometry();
    const auto          number_of_integration_points =
        r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());

    auto& r_properties = this->GetProperties();
    rOutput.resize(number_of_integration_points);

    if (rVariable == VON_MISES_STRESS) {
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] =
                StressStrainUtilities::CalculateVonMisesStress(mStressVector[integration_point]);
        }
    } else if (rVariable == MEAN_EFFECTIVE_STRESS) {
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] =
                StressStrainUtilities::CalculateMeanStress(mStressVector[integration_point]);
        }
    } else if (rVariable == MEAN_STRESS) {
        std::vector<Vector> stress_vector;
        CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, stress_vector, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] =
                StressStrainUtilities::CalculateMeanStress(stress_vector[integration_point]);
        }
    } else if (rVariable == ENGINEERING_VON_MISES_STRAIN) {
        std::vector<Vector> strain_vector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] =
                StressStrainUtilities::CalculateVonMisesStrain(strain_vector[integration_point]);
        }
    } else if (rVariable == ENGINEERING_VOLUMETRIC_STRAIN) {
        std::vector<Vector> strain_vector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] =
                StressStrainUtilities::CalculateTrace(strain_vector[integration_point]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VON_MISES_STRAIN) {
        std::vector<Vector> strain_vector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] =
                StressStrainUtilities::CalculateVonMisesStrain(strain_vector[integration_point]);
        }
    } else if (rVariable == GREEN_LAGRANGE_VOLUMETRIC_STRAIN) {
        std::vector<Vector> strain_vector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] =
                StressStrainUtilities::CalculateTrace(strain_vector[integration_point]);
        }
    } else if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION ||
               rVariable == BISHOP_COEFFICIENT || rVariable == DERIVATIVE_OF_SATURATION ||
               rVariable == RELATIVE_PERMEABILITY) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        RetentionLaw::Parameters RetentionParameters(r_properties);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            // Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, integration_point);

            RetentionParameters.SetFluidPressure(GeoTransportEquationUtilities::CalculateFluidPressure(
                Variables.Np, Variables.PressureVector));
            rOutput[integration_point] = mRetentionLawVector[integration_point]->CalculateValue(
                RetentionParameters, rVariable, rOutput[integration_point]);
        }
    } else if (rVariable == HYDRAULIC_HEAD) {
        // Defining the shape functions, the Jacobian and the shape functions local gradients containers
        const Matrix& n_container = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        const auto nodal_hydraulic_head =
            GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(r_geometry, r_properties);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            const auto& shape_function_values = row(n_container, integration_point);
            rOutput[integration_point] =
                std::inner_product(shape_function_values.begin(), shape_function_values.end(),
                                   nodal_hydraulic_head.begin(), 0.0);
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

        ConstitutiveLaw::Parameters ConstitutiveParameters(r_geometry, r_properties, rCurrentProcessInfo);
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
    } else if (rVariable == GEO_SHEAR_CAPACITY) {
        OutputUtilities::CalculateShearCapacityValues(mStressVector, rOutput.begin(), r_properties);
    } else if (r_properties.Has(rVariable)) {
        // Map initial material property to gauss points, as required for the output
        std::fill_n(rOutput.begin(), number_of_integration_points, r_properties.GetValue(rVariable));
    } else {
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] = mConstitutiveLawVector[integration_point]->GetValue(
                rVariable, rOutput[integration_point]);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto number_of_integration_points =
        this->GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);
    rOutput.resize(number_of_integration_points);

    if (rVariable == FLUID_FLUX_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        const auto strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, this->GetStressStatePolicy().GetVoigtSize());

        const auto fluid_fluxes = GeoTransportEquationUtilities::CalculateFluidFluxes<TDim, TNumNodes>(
            this->GetGeometry(), this->GetIntegrationMethod(), this->GetProperties(), mRetentionLawVector,
            GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(
                strain_vectors, this->GetProperties()));
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            GeoElementUtilities::FillArray1dOutput(rOutput[integration_point], fluid_fluxes[integration_point]);
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
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point].resize(mStressVector[integration_point].size(), false);
            rOutput[integration_point] = mStressVector[integration_point];
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

        for (unsigned int integration_point = 0; integration_point < mConstitutiveLawVector.size();
             ++integration_point) {
            rOutput[integration_point] =
                mStressVector[integration_point] +
                PORE_PRESSURE_SIGN_FACTOR * biot_coefficients[integration_point] *
                    bishop_coefficients[integration_point] * fluid_pressures[integration_point] *
                    this->GetStressStatePolicy().GetVoigtVector();
        }
    } else if (rVariable == ENGINEERING_STRAIN_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < mConstitutiveLawVector.size();
             ++integration_point) {
            noalias(Variables.Np) = row(Variables.NContainer, integration_point);

            Matrix J0;
            Matrix InvJ0;
            this->CalculateDerivativesOnInitialConfiguration(Variables.detJInitialConfiguration, J0,
                                                             InvJ0, Variables.GradNpTInitialConfiguration,
                                                             integration_point);

            // Calculating operator B
            Variables.B = this->CalculateBMatrix(Variables.GradNpTInitialConfiguration, Variables.Np);

            // Compute infinitesimal strain
            Variables.StrainVector =
                StressStrainUtilities::CalculateCauchyStrain(Variables.B, Variables.DisplacementVector);
            rOutput[integration_point].resize(Variables.StrainVector.size(), false);
            rOutput[integration_point] = Variables.StrainVector;
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
        // Map initial material property to integration points, as required for the output
        std::fill_n(rOutput.begin(), mConstitutiveLawVector.size(), r_properties.GetValue(rVariable));
    } else {
        for (unsigned int integration_point = 0; integration_point < mConstitutiveLawVector.size(); ++integration_point)
            rOutput[integration_point] = mConstitutiveLawVector[integration_point]->GetValue(
                rVariable, rOutput[integration_point]);
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
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] =
                MathUtils<double>::StressVectorToTensor(mStressVector[integration_point]);
        }
    } else if (rVariable == TOTAL_STRESS_TENSOR) {
        std::vector<Vector> StressVector;
        this->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < mConstitutiveLawVector.size();
             ++integration_point) {
            rOutput[integration_point] =
                MathUtils<double>::StressVectorToTensor(StressVector[integration_point]);
        }
    } else if (rVariable == ENGINEERING_STRAIN_TENSOR) {
        std::vector<Vector> strain_vector;
        CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < mConstitutiveLawVector.size();
             ++integration_point) {
            rOutput[integration_point] =
                MathUtils<double>::StrainVectorToTensor(strain_vector[integration_point]);
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        std::vector<Vector> strain_vector;
        CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo);

        for (unsigned int integration_point = 0; integration_point < mConstitutiveLawVector.size();
             ++integration_point) {
            rOutput[integration_point] =
                MathUtils<double>::StrainVectorToTensor(strain_vector[integration_point]);
        }
    } else if (rVariable == PERMEABILITY_MATRIX) {
        // If the permeability of the element is a given property
        BoundedMatrix<double, TDim, TDim> PermeabilityMatrix;
        GeoElementUtilities::FillPermeabilityMatrix(PermeabilityMatrix, this->GetProperties());
        std::fill_n(rOutput.begin(), number_of_integration_points, PermeabilityMatrix);
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
    const PropertiesType&                           r_properties = this->GetProperties();
    const GeometryType&                             r_geometry   = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(mThisIntegrationMethod);

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geometry, r_properties, rCurrentProcessInfo);
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

    Matrix reduced_kg_matrix =
        prod(rGradNpt, Matrix(prod(MathUtils<double>::StressVectorToTensor(rStressVector), trans(rGradNpt)))) *
        IntegrationCoefficient;

    Matrix geometric_stiffness_matrix = ZeroMatrix(TNumNodes * TDim, TNumNodes * TDim);
    MathUtils<double>::ExpandAndAddReducedMatrix(geometric_stiffness_matrix, reduced_kg_matrix, TDim);

    // Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, geometric_stiffness_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geometry         = this->GetGeometry();
    const auto          integration_method = this->GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geometry.IntegrationPoints(integration_method);
    const auto N_container = r_geometry.ShapeFunctionsValues(integration_method);

    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        N_container, this->GetPressureSolutionVector());
    const auto degrees_saturation = this->CalculateDegreesOfSaturation(fluid_pressures);

    const auto solid_densities =
        GeoTransportEquationUtilities::CalculateSoilDensities(degrees_saturation, this->GetProperties());

    const auto det_Js_initial_configuration =
        GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(r_geometry, integration_method);

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(integration_points, det_Js_initial_configuration);

    const auto mass_matrix_u = GeoEquationOfMotionUtilities::CalculateMassMatrix(
        r_geometry.WorkingSpaceDimension(), r_geometry.PointsNumber(), integration_points.size(),
        r_geometry.ShapeFunctionsValues(integration_method), solid_densities, integration_coefficients);

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

    const PropertiesType&                           r_properties = this->GetProperties();
    const GeometryType&                             r_geometry   = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geometry, r_properties, rCurrentProcessInfo);

    // Stiffness matrix is needed to calculate Biot coefficient
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    RetentionLaw::Parameters RetentionParameters(r_properties);

    const auto b_matrices = CalculateBMatrices(Variables.DN_DXContainer, Variables.NContainer);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);
    const auto det_Js_initial_configuration = GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(
        r_geometry, this->GetIntegrationMethod());
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
        biot_coefficients, degrees_of_saturation, derivatives_of_saturation, r_properties);
    auto relative_permeability_values = this->CalculateRelativePermeabilityValues(fluid_pressures);
    const auto permeability_update_factors = GetOptionalPermeabilityUpdateFactors(strain_vectors);
    std::ranges::transform(permeability_update_factors, relative_permeability_values,
                           relative_permeability_values.begin(), std::multiplies<>{});

    const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);

    for (unsigned int integration_point = 0; integration_point < IntegrationPoints.size(); ++integration_point) {
        this->CalculateKinematics(Variables, integration_point);
        Variables.B                  = b_matrices[integration_point];
        Variables.F                  = deformation_gradients[integration_point];
        Variables.StrainVector       = strain_vectors[integration_point];
        Variables.ConstitutiveMatrix = constitutive_matrices[integration_point];

        // Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, integration_point);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, integration_point);

        Variables.RelativePermeability = relative_permeability_values[integration_point];
        Variables.BishopCoefficient    = bishop_coefficients[integration_point];

        Variables.BiotCoefficient        = biot_coefficients[integration_point];
        Variables.BiotModulusInverse     = biot_moduli_inverse[integration_point];
        Variables.DegreeOfSaturation     = degrees_of_saturation[integration_point];
        Variables.IntegrationCoefficient = integration_coefficients[integration_point];
        Variables.IntegrationCoefficientInitialConfiguration =
            integration_coefficients_on_initial_configuration[integration_point];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, integration_point);
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
    result.reserve(rFluidPressures.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};
    std::transform(rFluidPressures.begin(), rFluidPressures.end(), mRetentionLawVector.begin(),
                   std::back_inserter(result),
                   [&retention_law_params](auto fluid_pressure, const auto& pRetentionLaw) {
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
    result.reserve(rFluidPressures.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};
    std::transform(rFluidPressures.begin(), rFluidPressures.end(), mRetentionLawVector.begin(),
                   std::back_inserter(result),
                   [&retention_law_params](auto fluid_pressure, const auto& pRetentionLaw) {
        retention_law_params.SetFluidPressure(fluid_pressure);
        return pRetentionLaw->CalculateSaturation(retention_law_params);
    });

    return result;
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

    const GeometryType& r_geometry = this->GetGeometry();
    const IndexType number_of_integration_points = r_geometry.IntegrationPointsNumber(mThisIntegrationMethod);

    // Shape functions
    rVariables.NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

    // Gradient of shape functions and determinant of Jacobian
    rVariables.detJContainer.resize(number_of_integration_points, false);

    r_geometry.ShapeFunctionsIntegrationPointsGradients(
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
    result.reserve(rDN_DXContainer.size());
    for (unsigned int integration_point = 0; integration_point < rDN_DXContainer.size(); ++integration_point) {
        result.push_back(this->CalculateBMatrix(rDN_DXContainer[integration_point],
                                                row(rNContainer, integration_point)));
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                ElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);

    this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);

    if (!rVariables.IgnoreUndrained) {
        const auto permeability_matrix =
            GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
                rVariables.GradNpT, rVariables.DynamicViscosityInverse, rVariables.PermeabilityMatrix,
                rVariables.RelativePermeability, rVariables.IntegrationCoefficient);
        GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, permeability_matrix);

        this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix,
                                                                            const ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto stiffness_matrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrixGPoint<TDim, TNumNodes>(
        rVariables.B, rVariables.ConstitutiveMatrix, rVariables.IntegrationCoefficient);

    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, stiffness_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix,
                                                                           const ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix<TDim, TNumNodes>(
        rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
        rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);
    GeoElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix, coupling_matrix);

    if (!rVariables.IgnoreUndrained) {
        const auto p_coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix<TDim, TNumNodes>(
            rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
            rVariables.BiotCoefficient, rVariables.DegreeOfSaturation, rVariables.IntegrationCoefficient);
        GeoElementUtilities::AssemblePUBlockMatrix(
            rLeftHandSideMatrix,
            PORE_PRESSURE_SIGN_FACTOR * rVariables.VelocityCoefficient * trans(p_coupling_matrix));
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

    const auto coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix<TDim, TNumNodes>(
        rVariables.B, this->GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
        rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);

    const array_1d<double, TNumNodes * TDim> coupling_force = prod(coupling_matrix, rVariables.PressureVector);
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, (-1.0) * coupling_force);

    if (!rVariables.IgnoreUndrained) {
        const auto p_coupling_matrix = GeoTransportEquationUtilities::CalculateCouplingMatrix<TDim, TNumNodes>(
            rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
            rVariables.BiotCoefficient, rVariables.DegreeOfSaturation, rVariables.IntegrationCoefficient);
        const array_1d<double, TNumNodes> coupling_flow =
            PORE_PRESSURE_SIGN_FACTOR * prod(trans(p_coupling_matrix), rVariables.VelocityVector);
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, (-1.0) * coupling_flow);
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
    return RetentionLaw::CalculateRelativePermeabilityValues(
        mRetentionLawVector, this->GetProperties(), rFluidPressures);
}

template <unsigned int TDim, unsigned int TNumNodes>
std::vector<double> UPwSmallStrainElement<TDim, TNumNodes>::CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const
{
    KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

    auto result = std::vector<double>{};
    result.reserve(rFluidPressures.size());
    std::ranges::transform(mRetentionLawVector, rFluidPressures, std::back_inserter(result),
                           [&retention_law_params](const auto& pRetentionLaw, auto FluidPressure) {
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
    const auto number_of_integration_points =
        this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    std::vector<Matrix> result;
    result.reserve(number_of_integration_points);
    for (unsigned int integration_point = 0; integration_point < number_of_integration_points; ++integration_point) {
        result.push_back(CalculateDeformationGradient(integration_point));
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

    const GeometryType& r_geometry = this->GetGeometry();
    VariablesUtilities::GetNodalValues(r_geometry, WATER_PRESSURE, rVariables.PressureVector.begin());
    VariablesUtilities::GetNodalValues(r_geometry, DT_WATER_PRESSURE, rVariables.DtPressureVector.begin());

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNodalDisplacementVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& r_geometry = this->GetGeometry();

    // Nodal variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector,
                                                                 r_geometry, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector, r_geometry, VELOCITY);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeNodalVolumeAccelerationVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(
        rVariables.VolumeAcceleration, this->GetGeometry(), VOLUME_ACCELERATION);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::InitializeProperties(ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto& r_properties = this->GetProperties();

    rVariables.IgnoreUndrained = r_properties[IGNORE_UNDRAINED];
    rVariables.UseHenckyStrain = r_properties.Has(USE_HENCKY_STRAIN) ? r_properties[USE_HENCKY_STRAIN] : false;

    rVariables.ConsiderGeometricStiffness =
        r_properties.Has(CONSIDER_GEOMETRIC_STIFFNESS) ? r_properties[CONSIDER_GEOMETRIC_STIFFNESS] : false;

    rVariables.DynamicViscosityInverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
    GeoElementUtilities::FillPermeabilityMatrix(rVariables.PermeabilityMatrix, r_properties);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::CalculateKinematics(ElementVariables& rVariables,
                                                                 unsigned int IntegrationPointIndex)
{
    KRATOS_TRY

    // Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Np)      = row(rVariables.NContainer, IntegrationPointIndex);
    noalias(rVariables.GradNpT) = rVariables.DN_DXContainer[IntegrationPointIndex];

    rVariables.detJ = rVariables.detJContainer[IntegrationPointIndex];

    Matrix J0, InvJ0;
    this->CalculateDerivativesOnInitialConfiguration(rVariables.detJInitialConfiguration, J0, InvJ0,
                                                     rVariables.GradNpTInitialConfiguration,
                                                     IntegrationPointIndex);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
Vector UPwSmallStrainElement<TDim, TNumNodes>::GetPressureSolutionVector()
{
    Vector result(TNumNodes);
    std::ranges::transform(this->GetGeometry(), result.begin(), [](const auto& node) {
        return node.FastGetSolutionStepValue(WATER_PRESSURE);
    });
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

    for (unsigned int integration_point = 0; integration_point < rDeformationGradients.size(); ++integration_point) {
        // Explicitly convert from `row`'s return type to `Vector` to avoid ending up with a pointer
        // to an implicitly converted object
        const auto shape_function_values = Vector{row(rNuContainer, integration_point)};
        ConstitutiveLawUtilities::SetConstitutiveParameters(
            rConstitutiveParameters, rStrainVectors[integration_point],
            rConstitutiveMatrices[integration_point], shape_function_values,
            rDNu_DXContainer[integration_point], rDeformationGradients[integration_point],
            determinants_of_deformation_gradients[integration_point]);
        rConstitutiveParameters.SetStressVector(rStressVectors[integration_point]);

        mConstitutiveLawVector[integration_point]->CalculateMaterialResponseCauchy(rConstitutiveParameters);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string UPwSmallStrainElement<TDim, TNumNodes>::Info() const
{
    const std::string constitutive_info =
        !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";
    return "U-Pw small strain Element #" + std::to_string(this->Id()) + "\nConstitutive law: " + constitutive_info;
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, UPwBaseElement)
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwSmallStrainElement<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, UPwBaseElement)
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
