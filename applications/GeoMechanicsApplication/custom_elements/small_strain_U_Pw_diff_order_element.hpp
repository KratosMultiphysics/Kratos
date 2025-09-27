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

#pragma once

#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_retention/retention_law.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/output_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geometries/geometry_data.h"
#include "includes/cfd_variables.h"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/kratos_export_api.h"
#include "includes/smart_pointers.h"
#include "includes/ublas_interface.h"
#include "stress_state_policy.h"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace Kratos
{

template <typename...>
constexpr bool dependent_false = false;

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUPwDiffOrderElement : public UPwBaseElement
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallStrainUPwDiffOrderElement);

    using UPwBaseElement::UPwBaseElement;

    static constexpr std::size_t TNumPNodes = []() constexpr {
        if constexpr (TDim == 2) {
            if constexpr (TNumNodes == 6) return 3;   // 2D T6P3
            if constexpr (TNumNodes == 8) return 4;   // 2D Q8P4
            if constexpr (TNumNodes == 9) return 4;   // 2D Q9P4
            if constexpr (TNumNodes == 10) return 6;  // 2D T10P6
            if constexpr (TNumNodes == 15) return 10; // 2D T15P10
        } else if constexpr (TDim == 3) {
            if constexpr (TNumNodes == 10) return 4; // 3D T10P4
            if constexpr (TNumNodes == 20) return 8; // 3D T20P8
            if constexpr (TNumNodes == 27) return 8; // 3D T27P8
        } else {
            static_assert(dependent_false<std::integral_constant<std::size_t, TDim>>,
                          "The number of pressure nodes for the given element is not defined.");
        }
        return 0;
    }();

    SmallStrainUPwDiffOrderElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                   std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : UPwBaseElement(NewId, pGeometry, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
        SetUpPressureGeometryPointer();
        if (TNumPNodes != mpPressureGeometry->PointsNumber()) {
            KRATOS_ERROR << "The number of pressure nodes is not correct. Expected: " << TNumPNodes
                         << " - Given: " << mpPressureGeometry->PointsNumber() << std::endl;
        }
    }

    SmallStrainUPwDiffOrderElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   PropertiesType::Pointer            pProperties,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                   std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : UPwBaseElement(NewId, pGeometry, pProperties, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
        SetUpPressureGeometryPointer();
        if (TNumPNodes != mpPressureGeometry->PointsNumber()) {
            KRATOS_ERROR << "The number of pressure nodes is not correct. Expected: " << TNumPNodes
                         << " - Given: " << mpPressureGeometry->PointsNumber() << std::endl;
        }
    }

    ~SmallStrainUPwDiffOrderElement() override = default;

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<SmallStrainUPwDiffOrderElement>(
            NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
            this->CloneIntegrationCoefficientModifier());
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY

        if (const auto ierr = UPwBaseElement::Check(rCurrentProcessInfo); ierr != 0) return ierr;

        const auto& r_geom     = GetGeometry();
        const auto  element_Id = this->Id();

        CheckUtilities::CheckDomainSize(r_geom.DomainSize(), element_Id);

        // check pressure geometry pointer
        KRATOS_DEBUG_ERROR_IF_NOT(mpPressureGeometry) << "Pressure Geometry is not defined\n";

        const auto&           r_prop = this->GetProperties();
        const CheckProperties check_properties(r_prop, "parameter list", element_Id,
                                               CheckProperties::Bounds::AllExclusive);
        check_properties.CheckAvailability(IGNORE_UNDRAINED);
        if (!r_prop[IGNORE_UNDRAINED]) check_properties.CheckPermeabilityProperties(TDim);

        check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW);
        r_prop[CONSTITUTIVE_LAW]->Check(r_prop, r_geom, rCurrentProcessInfo);
        const auto expected_size = this->GetStressStatePolicy().GetVoigtSize();
        ConstitutiveLawUtilities::CheckStrainSize(r_prop, expected_size, element_Id);
        ConstitutiveLawUtilities::CheckHasStrainMeasure_Infinitesimal(r_prop, element_Id);

        return RetentionLaw::Check(mRetentionLawVector, r_prop, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        ElementVariables variables;
        this->InitializeElementVariables(variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(variables.DNu_DXContainer, variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        const auto determinants_of_deformation_gradients =
            GeoMechanicsMathUtilities::CalculateDeterminants(deformation_gradients);
        const auto strain_vectors = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, variables.DisplacementVector,
            variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());

        const auto number_of_integration_points =
            GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
        for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
            this->CalculateKinematics(variables, g_point);
            variables.B            = b_matrices[g_point];
            variables.F            = deformation_gradients[g_point];
            variables.StrainVector = strain_vectors[g_point];

            ConstitutiveLawUtilities::SetConstitutiveParameters(
                ConstitutiveParameters, variables.StrainVector, variables.ConstitutiveMatrix,
                variables.Nu, variables.DNu_DX, variables.F, determinants_of_deformation_gradients[g_point]);

            // Compute constitutive tensor and/or stresses
            noalias(variables.StressVector) = mStressVector[g_point];
            ConstitutiveParameters.SetStressVector(variables.StressVector);
            mConstitutiveLawVector[g_point]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
            mStateVariablesFinalized[g_point] = mConstitutiveLawVector[g_point]->GetValue(
                STATE_VARIABLES, mStateVariablesFinalized[g_point]);
        }

        // Assign pressure values to the intermediate nodes for post-processing
        if (!GetProperties()[IGNORE_UNDRAINED]) AssignPressureToIntermediateNodes();

        KRATOS_CATCH("")
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const GeometryType& r_geom             = GetGeometry();
        const auto          integration_method = this->GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points =
            r_geom.IntegrationPoints(integration_method);
        const auto Np_container = mpPressureGeometry->ShapeFunctionsValues(integration_method);

        const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
            Np_container, this->GetPressureSolutionVector());
        const auto degrees_saturation = this->CalculateDegreesOfSaturation(fluid_pressures);

        const auto solid_densities =
            GeoTransportEquationUtilities::CalculateSoilDensities(degrees_saturation, GetProperties());

        const auto det_Js_initial_configuration =
            GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(r_geom, integration_method);

        const auto integration_coefficients =
            this->CalculateIntegrationCoefficients(integration_points, det_Js_initial_configuration);

        const auto mass_matrix_u = GeoEquationOfMotionUtilities::CalculateMassMatrix(
            TDim, TNumNodes, integration_points.size(),
            r_geom.ShapeFunctionsValues(integration_method), solid_densities, integration_coefficients);

        rMassMatrix = ZeroMatrix(GetNumberOfDOF(), GetNumberOfDOF());
        GeoElementUtilities::AssembleUUBlockMatrix(rMassMatrix, mass_matrix_u);

        KRATOS_CATCH("")
    }

    void SetValuesOnIntegrationPoints(const Variable<Vector>&    rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo&         rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rVariable == CAUCHY_STRESS_VECTOR) {
            KRATOS_ERROR_IF(rValues.size() != GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod))
                << "Unexpected number of values for "
                   "SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints"
                << std::endl;
            mStressVector.resize(rValues.size());
            std::copy(rValues.begin(), rValues.end(), mStressVector.begin());
        } else {
            KRATOS_ERROR_IF(rValues.size() < mConstitutiveLawVector.size())
                << "Insufficient number of values for "
                   "SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints"
                << std::endl;
            for (unsigned int g_point = 0; g_point < mConstitutiveLawVector.size(); ++g_point) {
                mConstitutiveLawVector[g_point]->SetValue(rVariable, rValues[g_point], rCurrentProcessInfo);
            }
        }

        KRATOS_CATCH("")
    }

    using Element::SetValuesOnIntegrationPoints;

    void CalculateOnIntegrationPoints(const Variable<int>& rVariable,
                                      std::vector<int>&    rValues,
                                      const ProcessInfo&   rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const auto number_of_integration_points =
            GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

        rValues.resize(number_of_integration_points);
        for (auto i = SizeType{0}; i < number_of_integration_points; ++i) {
            rValues[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
        }

        KRATOS_CATCH("")
    }

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const auto& r_geom       = GetGeometry();
        const auto& r_properties = this->GetProperties();
        const auto  number_of_integration_points =
            r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());

        rOutput.resize(number_of_integration_points);

        auto for_each_integration_point = [](std::size_t number_of_points, auto&& body) {
            for (unsigned int g_point = 0; g_point < number_of_points; ++g_point) {
                body(g_point);
            }
        };

        auto fill_from_stress_vector = [this, &rOutput, number_of_integration_points,
                                        &for_each_integration_point](auto func) {
            for_each_integration_point(number_of_integration_points, [this, &rOutput, &func](unsigned int g_point) {
                rOutput[g_point] = func(mStressVector[g_point]);
            });
        };

        auto fill_from_vector = [this, &rOutput, number_of_integration_points, &rCurrentProcessInfo,
                                 &for_each_integration_point](const Variable<Vector>& rTensorVariable, auto func) {
            std::vector<Vector> tensor_vector;
            CalculateOnIntegrationPoints(rTensorVariable, tensor_vector, rCurrentProcessInfo);
            for_each_integration_point(number_of_integration_points,
                                       [&rOutput, &tensor_vector, &func](unsigned int g_point) {
                rOutput[g_point] = func(tensor_vector[g_point]);
            });
        };

        if (rVariable == VON_MISES_STRESS)
            return fill_from_stress_vector(StressStrainUtilities::CalculateVonMisesStress);

        if (rVariable == MEAN_EFFECTIVE_STRESS)
            return fill_from_stress_vector(StressStrainUtilities::CalculateMeanStress);

        if (rVariable == MEAN_STRESS)
            return fill_from_vector(TOTAL_STRESS_VECTOR, StressStrainUtilities::CalculateMeanStress);

        if (rVariable == ENGINEERING_VON_MISES_STRAIN)
            return fill_from_vector(ENGINEERING_STRAIN_VECTOR, StressStrainUtilities::CalculateVonMisesStrain);

        if (rVariable == ENGINEERING_VOLUMETRIC_STRAIN)
            return fill_from_vector(ENGINEERING_STRAIN_VECTOR, StressStrainUtilities::CalculateTrace);

        if (rVariable == GREEN_LAGRANGE_VON_MISES_STRAIN)
            return fill_from_vector(GREEN_LAGRANGE_STRAIN_VECTOR, StressStrainUtilities::CalculateVonMisesStrain);

        if (rVariable == GREEN_LAGRANGE_VOLUMETRIC_STRAIN)
            return fill_from_vector(GREEN_LAGRANGE_STRAIN_VECTOR, StressStrainUtilities::CalculateTrace);

        if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION || rVariable == BISHOP_COEFFICIENT ||
            rVariable == DERIVATIVE_OF_SATURATION || rVariable == RELATIVE_PERMEABILITY) {
            ElementVariables variables;
            this->InitializeElementVariables(variables, rCurrentProcessInfo);

            RetentionLaw::Parameters retention_parameters(r_properties);

            for_each_integration_point(
                number_of_integration_points,
                [this, &variables, &rVariable, &retention_parameters, &rOutput](unsigned int g_point) {
                // Compute Np, GradNpT, B and StrainVector
                this->CalculateKinematics(variables, g_point);

                retention_parameters.SetFluidPressure(GeoTransportEquationUtilities::CalculateFluidPressure(
                    variables.Np, variables.PressureVector));

                rOutput[g_point] = mRetentionLawVector[g_point]->CalculateValue(
                    retention_parameters, rVariable, rOutput[g_point]);
            });
            return;
        }
        if (rVariable == HYDRAULIC_HEAD) {
            constexpr auto numerical_limit = std::numeric_limits<double>::epsilon();

            Vector nodal_hydraulic_head(TNumNodes, 0.0);
            for (unsigned int node = 0; node < TNumNodes; ++node) {
                const auto& r_volume_acceleration =
                    r_geom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
                const double g = norm_2(r_volume_acceleration);
                if (g <= numerical_limit) {
                    continue;
                }

                array_1d<double, 3> node_coordinates = ZeroVector(3);
                noalias(node_coordinates)            = r_geom[node].Coordinates();
                const array_1d<double, 3> volume_acceleration_unit_vector = r_volume_acceleration / g;
                const auto water_pressure = r_geom[node].FastGetSolutionStepValue(WATER_PRESSURE);
                const auto fluid_weight   = g * r_properties[DENSITY_WATER];
                nodal_hydraulic_head[node] = -inner_prod(node_coordinates, volume_acceleration_unit_vector) -
                                             PORE_PRESSURE_SIGN_FACTOR * water_pressure / fluid_weight;
            }

            const auto& r_n_container = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());
            for_each_integration_point(number_of_integration_points,
                                       [&r_n_container, &rOutput, &nodal_hydraulic_head](unsigned int g_point) {
                const auto& r_shape_function = row(r_n_container, g_point);
                rOutput[g_point] = std::inner_product(r_shape_function.begin(), r_shape_function.end(),
                                                      nodal_hydraulic_head.begin(), 0.0);
            });
            return;
        }

        if (rVariable == CONFINED_STIFFNESS || rVariable == SHEAR_STIFFNESS) {
            KRATOS_ERROR_IF(TDim != 2 && TDim != 3) << rVariable.Name() << " can not be retrieved for dim "
                                                    << TDim << " in element: " << this->Id() << std::endl;

            const bool is_confined = (rVariable == CONFINED_STIFFNESS);

            constexpr std::array<std::size_t, 2> confined_indices = {
                static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_XX), static_cast<std::size_t>(INDEX_3D_XX)};

            constexpr std::array<std::size_t, 2> shear_indices = {
                static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_XY), static_cast<std::size_t>(INDEX_3D_XZ)};

            const std::size_t variable_index =
                is_confined ? confined_indices[TDim - 2] : shear_indices[TDim - 2];

            ElementVariables variables;
            this->InitializeElementVariables(variables, rCurrentProcessInfo);

            const auto b_matrices = CalculateBMatrices(variables.DNu_DXContainer, variables.NuContainer);
            const auto deformation_gradients = CalculateDeformationGradients();
            auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
                deformation_gradients, b_matrices, variables.DisplacementVector,
                variables.UseHenckyStrain, this->GetStressStatePolicy().GetVoigtSize());

            ConstitutiveLaw::Parameters constitutive_parameters(r_geom, r_properties, rCurrentProcessInfo);
            constitutive_parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
            constitutive_parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            std::vector<Matrix> constitutive_matrices;
            this->CalculateAnyOfMaterialResponse(deformation_gradients, constitutive_parameters,
                                                 variables.NuContainer, variables.DNu_DXContainer,
                                                 strain_vectors, mStressVector, constitutive_matrices);

            std::ranges::transform(constitutive_matrices, rOutput.begin(),
                                   [variable_index](const Matrix& constitutive_matrix) {
                return constitutive_matrix(variable_index, variable_index);
            });
            return;
        }
        if (r_properties.Has(rVariable)) {
            // Map initial material property to gauss points, as required for the output
            std::fill_n(rOutput.begin(), number_of_integration_points, r_properties.GetValue(rVariable));
            return;
        }
        if (rVariable == GEO_SHEAR_CAPACITY) {
            OutputUtilities::CalculateShearCapacityValues(mStressVector, rOutput.begin(), r_properties);
            return;
        }

        for_each_integration_point(number_of_integration_points, [this, &rOutput, &rVariable](unsigned int g_point) {
            rOutput[g_point] = mConstitutiveLawVector[g_point]->GetValue(rVariable, rOutput[g_point]);
        });

        KRATOS_CATCH("")
    }

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const GeometryType& r_geom = GetGeometry();
        const auto          number_of_integration_points =
            r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        rOutput.resize(number_of_integration_points);

        if (rVariable == CAUCHY_STRESS_VECTOR) {
            for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
                if (rOutput[g_point].size() != mStressVector[g_point].size())
                    rOutput[g_point].resize(mStressVector[g_point].size(), false);

                rOutput[g_point] = mStressVector[g_point];
            }
        } else if (rVariable == TOTAL_STRESS_VECTOR) {
            ElementVariables variables;
            this->InitializeElementVariables(variables, rCurrentProcessInfo);

            ConstitutiveLaw::Parameters ConstitutiveParameters(r_geom, GetProperties(), rCurrentProcessInfo);
            ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
            ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

            const auto b_matrices = CalculateBMatrices(variables.DNu_DXContainer, variables.NuContainer);
            const auto deformation_gradients = CalculateDeformationGradients();
            auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
                deformation_gradients, b_matrices, variables.DisplacementVector,
                variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
            std::vector<Matrix> constitutive_matrices;
            this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                                 variables.NuContainer, variables.DNu_DXContainer,
                                                 strain_vectors, mStressVector, constitutive_matrices);
            const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
                constitutive_matrices, this->GetProperties());
            const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
                variables.NpContainer, variables.PressureVector);
            const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);

            for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
                rOutput[g_point] = mStressVector[g_point] +
                                   PORE_PRESSURE_SIGN_FACTOR * biot_coefficients[g_point] *
                                       bishop_coefficients[g_point] * fluid_pressures[g_point] *
                                       GetStressStatePolicy().GetVoigtVector();
            }
        } else if (rVariable == ENGINEERING_STRAIN_VECTOR) {
            ElementVariables variables;
            this->InitializeElementVariables(variables, rCurrentProcessInfo);

            for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
                noalias(variables.Nu) = row(variables.NuContainer, g_point);

                Matrix J0;
                Matrix InvJ0;
                this->CalculateDerivativesOnInitialConfiguration(
                    variables.detJInitialConfiguration, J0, InvJ0, variables.DNu_DXInitialConfiguration, g_point);

                // Calculating operator B
                variables.B = this->CalculateBMatrix(variables.DNu_DXInitialConfiguration, variables.Nu);

                // Compute infinitesimal strain
                variables.StrainVector = StressStrainUtilities::CalculateCauchyStrain(
                    variables.B, variables.DisplacementVector);

                if (rOutput[g_point].size() != variables.StrainVector.size())
                    rOutput[g_point].resize(variables.StrainVector.size(), false);

                rOutput[g_point] = variables.StrainVector;
            }
        } else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
            ElementVariables variables;
            this->InitializeElementVariables(variables, rCurrentProcessInfo);

            const auto b_matrices = CalculateBMatrices(variables.DNu_DXContainer, variables.NuContainer);
            const auto deformation_gradients = CalculateDeformationGradients();
            rOutput                          = StressStrainUtilities::CalculateStrains(
                deformation_gradients, b_matrices, variables.DisplacementVector,
                variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
        } else {
            for (unsigned int i = 0; i < number_of_integration_points; ++i)
                rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }

        KRATOS_CATCH("")
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const auto number_of_integration_points =
            GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
        rOutput.resize(number_of_integration_points);

        if (rVariable == FLUID_FLUX_VECTOR) {
            ElementVariables variables;
            this->InitializeElementVariables(variables, rCurrentProcessInfo);

            const auto b_matrices = CalculateBMatrices(variables.DNu_DXContainer, variables.NuContainer);
            const auto deformation_gradients = CalculateDeformationGradients();
            const auto strain_vectors        = StressStrainUtilities::CalculateStrains(
                deformation_gradients, b_matrices, variables.DisplacementVector,
                variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
            auto relative_permeability_values =
                CalculateRelativePermeabilityValues(GeoTransportEquationUtilities::CalculateFluidPressures(
                    variables.NpContainer, variables.PressureVector));
            const auto permeability_update_factors =
                GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(strain_vectors,
                                                                                  GetProperties());
            std::ranges::transform(relative_permeability_values, permeability_update_factors,
                                   relative_permeability_values.begin(), std::multiplies<>{});

            // Loop over integration points
            for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
                // compute element kinematics (Np, gradNpT, |J|, B, strains)
                this->CalculateKinematics(variables, g_point);
                variables.B = b_matrices[g_point];

                const auto body_acceleration =
                    CalculateBodyAcceleration(variables.Nu, variables.BodyAcceleration);
                const auto relative_permeability = relative_permeability_values[g_point];

                Vector grad_pressure_term(TDim);
                noalias(grad_pressure_term) = prod(trans(variables.DNp_DX), variables.PressureVector);
                noalias(grad_pressure_term) +=
                    PORE_PRESSURE_SIGN_FACTOR * GetProperties()[DENSITY_WATER] * body_acceleration;

                // Compute fluid flux vector q [L/T]
                rOutput[g_point].clear();
                const auto fluid_flux = PORE_PRESSURE_SIGN_FACTOR *
                                        variables.DynamicViscosityInverse * relative_permeability *
                                        prod(variables.IntrinsicPermeability, grad_pressure_term);
                std::copy_n(fluid_flux.begin(), TDim, rOutput[g_point].begin());
            }
        }

        KRATOS_CATCH("")
    }

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const auto number_of_integration_points =
            GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        rOutput.resize(number_of_integration_points);

        if (rVariable == CAUCHY_STRESS_TENSOR) {
            std::vector<Vector> StressVector;
            this->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

            for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
                rOutput[g_point] = MathUtils<double>::StressVectorToTensor(StressVector[g_point]);
            }
        } else if (rVariable == TOTAL_STRESS_TENSOR) {
            std::vector<Vector> StressVector;
            this->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

            for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
                rOutput[g_point] = MathUtils<double>::StressVectorToTensor(StressVector[g_point]);
            }
        } else if (rVariable == ENGINEERING_STRAIN_TENSOR) {
            std::vector<Vector> StrainVector;
            CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

            for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
                rOutput[g_point] = MathUtils<double>::StrainVectorToTensor(StrainVector[g_point]);
            }
        } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
            std::vector<Vector> StrainVector;
            CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

            for (unsigned int g_point = 0; g_point < number_of_integration_points; ++g_point) {
                rOutput[g_point] = MathUtils<double>::StrainVectorToTensor(StrainVector[g_point]);
            }
        } else {
            for (unsigned int i = 0; i < number_of_integration_points; ++i) {
                rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
            }
        }

        KRATOS_CATCH("")
    }

    using Element::CalculateOnIntegrationPoints;

    // Turn back information as a string.
    std::string Info() const override
    {
        const std::string constitutive_info =
            !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";
        return "U-Pw small strain different order Element #" + std::to_string(Id()) +
               "\nConstitutive law: " + constitutive_info;
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

protected:

    using UPwBaseElement::mConstitutiveLawVector;
    using UPwBaseElement::mRetentionLawVector;
    using UPwBaseElement::mStateVariablesFinalized;
    using UPwBaseElement::mStressVector;

    struct ElementVariables {
        // variables at all integration points
        Matrix                                    NuContainer;
        Matrix                                    NpContainer;
        GeometryType::ShapeFunctionsGradientsType DNu_DXContainer;
        GeometryType::ShapeFunctionsGradientsType DNp_DXContainer;
        Vector                                    detJuContainer;

        // variables at each integration point
        Vector Nu;     // Contains the displacement shape functions at every node
        Vector Np;     // Contains the pressure shape functions at every node
        Matrix DNu_DX; // Contains the global derivatives of the displacement shape functions
        Matrix DNu_DXInitialConfiguration; // Contains the global derivatives of the displacement shape functions

        Matrix DNp_DX; // Contains the global derivatives of the pressure shape functions
        Matrix B;
        double IntegrationCoefficient;
        double IntegrationCoefficientInitialConfiguration;
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;

        // variables needed for consistency with the general constitutive law
        Matrix F;

        // needed for updated Lagrangian:
        double detJ;                     // displacement
        double detJInitialConfiguration; // displacement

        // Nodal variables
        Vector BodyAcceleration;
        Vector DisplacementVector;
        Vector VelocityVector;
        Vector PressureVector;
        Vector DeltaPressureVector;
        Vector PressureDtVector;

        /// Retention Law parameters
        double DegreeOfSaturation;
        double DerivativeOfSaturation;
        double RelativePermeability;
        double BishopCoefficient;

        // Properties and processinfo variables
        bool IgnoreUndrained;
        bool UseHenckyStrain;
        bool ConsiderGeometricStiffness;

        // stress/flow variables
        double BiotCoefficient;
        double BiotModulusInverse;
        double DynamicViscosityInverse;
        Matrix IntrinsicPermeability;
        double VelocityCoefficient;
        double DtPressureCoefficient;
    };

    void CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const GeometryType& r_geom = GetGeometry();

        // Definition of variables
        ElementVariables variables;
        this->InitializeElementVariables(variables, rCurrentProcessInfo);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(r_geom, GetProperties(), rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        // Loop over integration points
        const GeometryType::IntegrationPointsArrayType& r_integration_points =
            r_geom.IntegrationPoints(this->GetIntegrationMethod());

        const auto b_matrices = CalculateBMatrices(variables.DNu_DXContainer, variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, variables.DisplacementVector,
            variables.UseHenckyStrain, this->GetStressStatePolicy().GetVoigtSize());
        std::vector<Matrix> constitutive_matrices;
        this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                             variables.NuContainer, variables.DNu_DXContainer,
                                             strain_vectors, mStressVector, constitutive_matrices);
        const auto integration_coefficients =
            this->CalculateIntegrationCoefficients(r_integration_points, variables.detJuContainer);

        const auto stiffness_matrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
            b_matrices, constitutive_matrices, integration_coefficients);

        GeoElementUtilities::AssembleUUBlockMatrix(rStiffnessMatrix, stiffness_matrix);

        KRATOS_CATCH("")
    }

    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      bool               CalculateStiffnessMatrixFlag,
                      bool               CalculateResidualVectorFlag) override
    {
        KRATOS_TRY

        const auto&                                     r_prop = this->GetProperties();
        const auto&                                     r_geom = GetGeometry();
        const GeometryType::IntegrationPointsArrayType& r_integration_points =
            r_geom.IntegrationPoints(this->GetIntegrationMethod());

        ConstitutiveLaw::Parameters ConstitutiveParameters(r_geom, r_prop, rCurrentProcessInfo);

        // Stiffness matrix is needed to calculate Biot coefficient
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        if (CalculateResidualVectorFlag)
            ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        ElementVariables variables;
        this->InitializeElementVariables(variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(variables.DNu_DXContainer, variables.NuContainer);
        const auto integration_coefficients =
            this->CalculateIntegrationCoefficients(r_integration_points, variables.detJuContainer);

        const auto det_Js_initial_configuration = GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(
            r_geom, this->GetIntegrationMethod());

        const auto integration_coefficients_on_initial_configuration =
            this->CalculateIntegrationCoefficients(r_integration_points, det_Js_initial_configuration);

        const auto deformation_gradients = CalculateDeformationGradients();
        auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, variables.DisplacementVector,
            variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
        std::vector<Matrix> constitutive_matrices;
        this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                             variables.NuContainer, variables.DNu_DXContainer,
                                             strain_vectors, mStressVector, constitutive_matrices);
        const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
            constitutive_matrices, this->GetProperties());
        const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
            variables.NpContainer, variables.PressureVector);
        const auto degrees_of_saturation     = CalculateDegreesOfSaturation(fluid_pressures);
        const auto derivatives_of_saturation = CalculateDerivativesOfSaturation(fluid_pressures);
        const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
            biot_coefficients, degrees_of_saturation, derivatives_of_saturation, r_prop);
        auto relative_permeability_values = CalculateRelativePermeabilityValues(fluid_pressures);
        const auto permeability_update_factors = GetOptionalPermeabilityUpdateFactors(strain_vectors);
        std::ranges::transform(permeability_update_factors, relative_permeability_values,
                               relative_permeability_values.begin(), std::multiplies<>{});

        const auto bishop_coefficients = CalculateBishopCoefficients(fluid_pressures);

        for (unsigned int g_point = 0; g_point < r_integration_points.size(); ++g_point) {
            this->CalculateKinematics(variables, g_point);
            variables.B                  = b_matrices[g_point];
            variables.F                  = deformation_gradients[g_point];
            variables.StrainVector       = strain_vectors[g_point];
            variables.ConstitutiveMatrix = constitutive_matrices[g_point];

            variables.RelativePermeability = relative_permeability_values[g_point];
            variables.BishopCoefficient    = bishop_coefficients[g_point];

            variables.BiotCoefficient        = biot_coefficients[g_point];
            variables.BiotModulusInverse     = biot_moduli_inverse[g_point];
            variables.DegreeOfSaturation     = degrees_of_saturation[g_point];
            variables.IntegrationCoefficient = integration_coefficients[g_point];

            variables.IntegrationCoefficientInitialConfiguration =
                integration_coefficients_on_initial_configuration[g_point];

            // Contributions to the left hand side
            if (CalculateStiffnessMatrixFlag)
                this->CalculateAndAddLHS(rLeftHandSideMatrix, variables);

            // Contributions to the right hand side
            if (CalculateResidualVectorFlag)
                this->CalculateAndAddRHS(rRightHandSideVector, variables, g_point);
        }

        KRATOS_CATCH("")
    }

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const GeometryType& r_geom  = GetGeometry();
        const SizeType num_g_points = r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());

        // variables at all integration points
        rVariables.NuContainer.resize(num_g_points, TNumNodes, false);
        rVariables.NuContainer = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());

        rVariables.NpContainer.resize(num_g_points, TNumPNodes, false);
        rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues(this->GetIntegrationMethod());

        rVariables.Nu.resize(TNumNodes, false);
        rVariables.Np.resize(TNumPNodes, false);

        rVariables.DNu_DXContainer.resize(num_g_points, false);
        for (SizeType i = 0; i < num_g_points; ++i)
            ((rVariables.DNu_DXContainer)[i]).resize(TNumNodes, TDim, false);
        rVariables.DNu_DX.resize(TNumNodes, TDim, false);
        rVariables.DNu_DXInitialConfiguration.resize(TNumNodes, TDim, false);
        rVariables.detJuContainer.resize(num_g_points, false);
        r_geom.ShapeFunctionsIntegrationPointsGradients(
            rVariables.DNu_DXContainer, rVariables.detJuContainer, this->GetIntegrationMethod());

        (rVariables.DNp_DXContainer).resize(num_g_points, false);
        for (SizeType i = 0; i < num_g_points; ++i)
            ((rVariables.DNp_DXContainer)[i]).resize(TNumPNodes, TDim, false);
        (rVariables.DNp_DX).resize(TNumPNodes, TDim, false);
        Vector detJpContainer = ZeroVector(num_g_points);
        mpPressureGeometry->ShapeFunctionsIntegrationPointsGradients(
            rVariables.DNp_DXContainer, detJpContainer, this->GetIntegrationMethod());

        // variables computed at each integration point
        const SizeType VoigtSize = this->GetStressStatePolicy().GetVoigtSize();

        rVariables.B.resize(VoigtSize, TNumNodes * TDim, false);
        noalias(rVariables.B) = ZeroMatrix(VoigtSize, TNumNodes * TDim);

        rVariables.StrainVector.resize(VoigtSize, false);
        rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);

        rVariables.StressVector.resize(VoigtSize, false);

        // Needed parameters for consistency with the general constitutive law
        rVariables.F.resize(TDim, TDim, false);
        noalias(rVariables.F) = identity_matrix<double>(TDim);

        // Nodal variables
        this->InitializeNodalVariables(rVariables);

        // Properties variables
        this->InitializeProperties(rVariables);

        // ProcessInfo variables
        rVariables.VelocityCoefficient   = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
        rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

        // Retention law
        rVariables.DegreeOfSaturation   = 1.0;
        rVariables.RelativePermeability = 1.0;
        rVariables.BishopCoefficient    = 1.0;

        KRATOS_CATCH("")
    }

    void InitializeNodalVariables(ElementVariables& rVariables)
    {
        KRATOS_TRY

        const GeometryType& r_geom = GetGeometry();

        Vector BodyAccelerationAux = ZeroVector(3);
        rVariables.BodyAcceleration.resize(TNumNodes * TDim, false);
        rVariables.DisplacementVector.resize(TNumNodes * TDim, false);
        rVariables.VelocityVector.resize(TNumNodes * TDim, false);

        for (SizeType i = 0; i < TNumNodes; ++i) {
            SizeType Local_i    = i * TDim;
            BodyAccelerationAux = r_geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

            rVariables.BodyAcceleration[Local_i] = BodyAccelerationAux[0];
            rVariables.DisplacementVector[Local_i] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
            rVariables.VelocityVector[Local_i] = r_geom[i].FastGetSolutionStepValue(VELOCITY_X);

            rVariables.BodyAcceleration[Local_i + 1] = BodyAccelerationAux[1];
            rVariables.DisplacementVector[Local_i + 1] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
            rVariables.VelocityVector[Local_i + 1] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Y);

            if constexpr (TDim > 2) {
                rVariables.BodyAcceleration[Local_i + 2] = BodyAccelerationAux[2];
                rVariables.DisplacementVector[Local_i + 2] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
                rVariables.VelocityVector[Local_i + 2] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Z);
            }
        }

        rVariables.PressureVector.resize(TNumPNodes, false);
        rVariables.PressureDtVector.resize(TNumPNodes, false);
        rVariables.DeltaPressureVector.resize(TNumPNodes, false);
        const auto& r_p_geometry = *mpPressureGeometry;
        for (SizeType i = 0; i < TNumPNodes; ++i) {
            rVariables.PressureVector[i] = r_p_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE);
            rVariables.PressureDtVector[i] = r_p_geometry[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
            rVariables.DeltaPressureVector[i] =
                r_p_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE) -
                r_p_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE, 1);
        }

        KRATOS_CATCH("")
    }

    void InitializeProperties(ElementVariables& rVariables)
    {
        KRATOS_TRY

        const PropertiesType& r_properties = this->GetProperties();

        rVariables.IgnoreUndrained = r_properties[IGNORE_UNDRAINED];
        rVariables.UseHenckyStrain =
            r_properties.Has(USE_HENCKY_STRAIN) ? r_properties[USE_HENCKY_STRAIN] : false;

        rVariables.ConsiderGeometricStiffness = r_properties.Has(CONSIDER_GEOMETRIC_STIFFNESS)
                                                    ? r_properties[CONSIDER_GEOMETRIC_STIFFNESS]
                                                    : false;

        rVariables.DynamicViscosityInverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
        // Setting the intrinsic permeability matrix
        rVariables.IntrinsicPermeability = GeoElementUtilities::FillPermeabilityMatrix(r_properties, TDim);

        KRATOS_CATCH("")
    }

    virtual void CalculateKinematics(ElementVariables& rVariables, unsigned int GPoint)
    {
        KRATOS_TRY

        // Setting the vector of shape functions and the matrix of the shape functions global gradients
        noalias(rVariables.Nu) = row(rVariables.NuContainer, GPoint);
        noalias(rVariables.Np) = row(rVariables.NpContainer, GPoint);

        noalias(rVariables.DNu_DX) = rVariables.DNu_DXContainer[GPoint];
        noalias(rVariables.DNp_DX) = rVariables.DNp_DXContainer[GPoint];

        rVariables.detJ = rVariables.detJuContainer[GPoint];

        Matrix J0;
        Matrix InvJ0;
        this->CalculateDerivativesOnInitialConfiguration(rVariables.detJInitialConfiguration, J0, InvJ0,
                                                         rVariables.DNu_DXInitialConfiguration, GPoint);

        KRATOS_CATCH("")
    }

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables) const
    {
        KRATOS_TRY

        this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);

        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);

        if (!rVariables.IgnoreUndrained) {
            const auto permeability_matrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix(
                rVariables.DNp_DX, rVariables.DynamicViscosityInverse, rVariables.IntrinsicPermeability,
                rVariables.RelativePermeability, rVariables.IntegrationCoefficient);
            GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, permeability_matrix);

            this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);
        }

        KRATOS_CATCH("")
    }

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables) const
    {
        KRATOS_TRY

        const auto dofs_per_node = rVariables.B.size2();
        Matrix     stiffness_matrix(dofs_per_node, dofs_per_node);

        GeoEquationOfMotionUtilities::CalculateStiffnessMatrixGPoint(
            stiffness_matrix, rVariables.B, rVariables.ConstitutiveMatrix, rVariables.IntegrationCoefficient);

        GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, stiffness_matrix);

        KRATOS_CATCH("")
    }

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables) const
    {
        KRATOS_TRY

        BoundedMatrix<double, TDim * TNumNodes, TNumPNodes> coupling_matrix;
        GeoTransportEquationUtilities::CalculateCouplingMatrix(
            coupling_matrix, rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
            rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);
        GeoElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix, coupling_matrix);

        if (!rVariables.IgnoreUndrained) {
            GeoTransportEquationUtilities::CalculateCouplingMatrix(
                coupling_matrix, rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
                rVariables.BiotCoefficient, rVariables.DegreeOfSaturation, rVariables.IntegrationCoefficient);
            GeoElementUtilities::AssemblePUBlockMatrix(
                rLeftHandSideMatrix,
                PORE_PRESSURE_SIGN_FACTOR * rVariables.VelocityCoefficient * trans(coupling_matrix));
        }

        KRATOS_CATCH("")
    }

    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables) const
    {
        KRATOS_TRY

        const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
            rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

        GeoElementUtilities::AssemblePPBlockMatrix(
            rLeftHandSideMatrix, compressibility_matrix * rVariables.DtPressureCoefficient);

        KRATOS_CATCH("")
    }

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint)
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

    void CalculateAndAddStiffnessForce(VectorType&             rRightHandSideVector,
                                       const ElementVariables& rVariables,
                                       unsigned int            GPoint)
    {
        KRATOS_TRY

        Vector stiffness_force =
            -1.0 * prod(trans(rVariables.B), mStressVector[GPoint]) * rVariables.IntegrationCoefficient;
        GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, stiffness_force);

        KRATOS_CATCH("")
    }

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
    {
        KRATOS_TRY

        const auto soil_density = GeoTransportEquationUtilities::CalculateSoilDensity(
            rVariables.DegreeOfSaturation, this->GetProperties());

        const auto body_acceleration = CalculateBodyAcceleration(rVariables.Nu, rVariables.BodyAcceleration);

        SizeType Index;
        for (SizeType i = 0; i < TNumNodes; ++i) {
            Index = i * TDim;
            for (SizeType idim = 0; idim < TDim; ++idim) {
                rRightHandSideVector[Index + idim] +=
                    rVariables.Nu[i] * soil_density * body_acceleration[idim] *
                    rVariables.IntegrationCoefficientInitialConfiguration;
            }
        }

        KRATOS_CATCH("")
    }

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, const ElementVariables& rVariables) const
    {
        KRATOS_TRY

        BoundedMatrix<double, TDim * TNumNodes, TNumPNodes> coupling_matrix;
        GeoTransportEquationUtilities::CalculateCouplingMatrix(
            coupling_matrix, rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
            rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);
        const Vector coupling_force = prod((-1.0) * coupling_matrix, rVariables.PressureVector);
        GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, coupling_force);

        if (!rVariables.IgnoreUndrained) {
            GeoTransportEquationUtilities::CalculateCouplingMatrix(
                coupling_matrix, rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
                rVariables.BiotCoefficient, rVariables.DegreeOfSaturation, rVariables.IntegrationCoefficient);
            const Vector coupling_flow =
                PORE_PRESSURE_SIGN_FACTOR * prod(trans((-1.0) * coupling_matrix), rVariables.VelocityVector);
            GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, coupling_flow);
        }

        KRATOS_CATCH("")
    }

    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, const ElementVariables& rVariables) const
    {
        KRATOS_TRY

        Matrix compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
            rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);
        Vector compressibility_flow = -prod(compressibility_matrix, rVariables.PressureDtVector);
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, compressibility_flow);

        KRATOS_CATCH("")
    }

    [[nodiscard]] std::vector<double> CalculateRelativePermeabilityValues(const std::vector<double>& rFluidPressures) const
    {
        KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

        auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

        auto result = std::vector<double>{};
        result.reserve(mRetentionLawVector.size());
        std::ranges::transform(mRetentionLawVector, rFluidPressures, std::back_inserter(result),
                               [&retention_law_params](const auto& pRetentionLaw, auto FluidPressure) {
            retention_law_params.SetFluidPressure(FluidPressure);
            return pRetentionLaw->CalculateRelativePermeability(retention_law_params);
        });
        return result;
    }

    [[nodiscard]] std::vector<double> CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const
    {
        KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

        auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

        auto result = std::vector<double>{};
        result.reserve(mRetentionLawVector.size());
        std::ranges::transform(mRetentionLawVector, rFluidPressures, std::back_inserter(result),
                               [&retention_law_params](const auto& pRetentionLaw, auto FluidPressure) {
            retention_law_params.SetFluidPressure(FluidPressure);
            return pRetentionLaw->CalculateBishopCoefficient(retention_law_params);
        });
        return result;
    }

    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, const ElementVariables& rVariables) const
    {
        KRATOS_TRY

        const Matrix permeability_matrix =
            -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse * rVariables.RelativePermeability *
            prod(rVariables.DNp_DX, Matrix(prod(rVariables.IntrinsicPermeability, trans(rVariables.DNp_DX)))) *
            rVariables.IntegrationCoefficient;
        const Vector permeability_flow = -prod(permeability_matrix, rVariables.PressureVector);
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, permeability_flow);

        KRATOS_CATCH("")
    }

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, const ElementVariables& rVariables)
    {
        KRATOS_TRY

        const Matrix grad_Np_T_perm =
            rVariables.DynamicViscosityInverse * rVariables.BishopCoefficient *
            GetProperties()[DENSITY_WATER] * rVariables.RelativePermeability *
            prod(rVariables.DNp_DX, rVariables.IntrinsicPermeability) * rVariables.IntegrationCoefficient;

        Vector body_acceleration = ZeroVector(TDim);

        SizeType index = 0;
        for (SizeType i = 0; i < TNumNodes; ++i) {
            for (SizeType idim = 0; idim < TDim; ++idim) {
                body_acceleration[idim] += rVariables.Nu[i] * rVariables.BodyAcceleration[index];
                index++;
            }
        }

        const Vector fluid_body_flow = prod(grad_Np_T_perm, body_acceleration);
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, fluid_body_flow);

        KRATOS_CATCH("")
    }

    Matrix CalculateBMatrix(const Matrix& rDN_DX, const Vector& rN) const
    {
        return this->GetStressStatePolicy().CalculateBMatrix(rDN_DX, rN, this->GetGeometry());
    }

    std::vector<Matrix> CalculateBMatrices(const GeometryType::ShapeFunctionsGradientsType& rDN_DXContainer,
                                           const Matrix& rNContainer) const
    {
        std::vector<Matrix> result;
        result.reserve(rDN_DXContainer.size());
        for (unsigned int g_point = 0; g_point < rDN_DXContainer.size(); ++g_point) {
            result.push_back(this->CalculateBMatrix(rDN_DXContainer[g_point], row(rNContainer, g_point)));
        }

        return result;
    }

    void AssignPressureToIntermediateNodes();

    virtual Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
    {
        return this->GetStressStatePolicy().CalculateGreenLagrangeStrain(rDeformationGradient);
    }

    Matrix CalculateDeformationGradient(unsigned int GPoint) const
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

    std::vector<Matrix> CalculateDeformationGradients() const
    {
        const auto number_of_integration_points =
            this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        std::vector<Matrix> result;
        result.reserve(number_of_integration_points);
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            result.push_back(CalculateDeformationGradient(integration_point));
        }

        return result;
    }

    ///
    /// \brief This function calculates the constitutive matrices, stresses and strains depending on the
    ///        constitutive parameters. Note that depending on the settings in the rConstitutiveParameters
    ///        the function could calculate the stress, the constitutive matrix, the strains, or a combination.
    ///        In our elements we generally always calculate the constitutive matrix and sometimes the stress.
    ///
    void CalculateAnyOfMaterialResponse(const std::vector<Matrix>&   rDeformationGradients,
                                        ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                        const Matrix&                rNuContainer,
                                        const GeometryType::ShapeFunctionsGradientsType& rDNu_DXContainer,
                                        std::vector<Vector>& rStrainVectors,
                                        std::vector<Vector>& rStressVectors,
                                        std::vector<Matrix>& rConstitutiveMatrices)
    {
        const SizeType voigt_size = TDim == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN;

        if (rStrainVectors.size() != rDeformationGradients.size()) {
            rStrainVectors.resize(rDeformationGradients.size());
            std::fill(rStrainVectors.begin(), rStrainVectors.end(), ZeroVector(voigt_size));
        }
        if (rStressVectors.size() != rDeformationGradients.size()) {
            rStressVectors.resize(rDeformationGradients.size());
            std::fill(rStressVectors.begin(), rStressVectors.end(), ZeroVector(voigt_size));
        }
        if (rConstitutiveMatrices.size() != rDeformationGradients.size()) {
            rConstitutiveMatrices.resize(rDeformationGradients.size());
            std::fill(rConstitutiveMatrices.begin(), rConstitutiveMatrices.end(),
                      ZeroMatrix(voigt_size, voigt_size));
        }

        const auto determinants_of_deformation_gradients =
            GeoMechanicsMathUtilities::CalculateDeterminants(rDeformationGradients);

        for (unsigned int g_point = 0; g_point < rDeformationGradients.size(); ++g_point) {
            // Explicitly convert from `row`'s return type to `Vector` to avoid ending up with a
            // pointer to an implicitly converted object
            const auto shape_function_values = Vector{row(rNuContainer, g_point)};
            ConstitutiveLawUtilities::SetConstitutiveParameters(
                rConstitutiveParameters, rStrainVectors[g_point], rConstitutiveMatrices[g_point],
                shape_function_values, rDNu_DXContainer[g_point], rDeformationGradients[g_point],
                determinants_of_deformation_gradients[g_point]);
            rConstitutiveParameters.SetStressVector(rStressVectors[g_point]);

            mConstitutiveLawVector[g_point]->CalculateMaterialResponseCauchy(rConstitutiveParameters);
        }
    }

    [[nodiscard]] Vector GetPressureSolutionVector() const
    {
        Vector result(TNumPNodes);
        std::transform(mpPressureGeometry->begin(), mpPressureGeometry->end(), result.begin(),
                       [](const auto& node) { return node.FastGetSolutionStepValue(WATER_PRESSURE); });
        return result;
    }

    [[nodiscard]] std::vector<double> CalculateDegreesOfSaturation(const std::vector<double>& rFluidPressures)
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

    [[nodiscard]] std::vector<double> CalculateDerivativesOfSaturation(const std::vector<double>& rFluidPressures)
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

    [[nodiscard]] virtual std::vector<double> GetOptionalPermeabilityUpdateFactors(const std::vector<Vector>& rStrainVectors) const
    {
        return GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(rStrainVectors,
                                                                                 GetProperties());
    }

    [[nodiscard]] SizeType GetNumberOfDOF() const override { return TNumNodes * TDim + TNumPNodes; }

private:
    GeometryType::Pointer mpPressureGeometry;

    [[nodiscard]] DofsVectorType GetDofs() const override
    {
        return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), *mpPressureGeometry, TDim);
    }

    /**
     * @brief Sets the up the pressure geometry pointer object
     * This function sets the pointer for the auxiliary geometry for the pressure problem
     * The pressure geometry pointer is set according to the element geometry number of nodes and dimension
     */
    void SetUpPressureGeometryPointer();

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, UPwBaseElement)
        rSerializer.save("PressureGeometry", mpPressureGeometry);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, UPwBaseElement)
        rSerializer.load("PressureGeometry", mpPressureGeometry);
    }

    // Private Operations

    template <class TValueType>
    inline void ThreadSafeNodeWrite(NodeType& rNode, const Variable<TValueType>& Var, const TValueType Value)
    {
        rNode.SetLock();
        rNode.FastGetSolutionStepValue(Var) = Value;
        rNode.UnSetLock();
    }

    Vector CalculateBodyAcceleration(Vector& rNu, Vector rBodyAcceleration) const
    {
        Vector body_acceleration = ZeroVector(TDim);
        for (SizeType i = 0; i < TNumNodes; ++i) {
            if constexpr (TDim == 2) {
                body_acceleration[0] += rNu[i] * rBodyAcceleration[i * 2 + 0];
                body_acceleration[1] += rNu[i] * rBodyAcceleration[i * 2 + 1];
            } else if constexpr (TDim == 3) {
                body_acceleration[0] += rNu[i] * rBodyAcceleration[i * 3 + 0];
                body_acceleration[1] += rNu[i] * rBodyAcceleration[i * 3 + 1];
                body_acceleration[2] += rNu[i] * rBodyAcceleration[i * 3 + 2];
            }
        }
        return body_acceleration;
    }
}; // Class SmallStrainUPwDiffOrderElement

} // namespace Kratos
