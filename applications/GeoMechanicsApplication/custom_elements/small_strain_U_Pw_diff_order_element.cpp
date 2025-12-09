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

// Project includes
#include "containers/array_1d.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "includes/cfd_variables.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/node_utilities.h"
#include "custom_utilities/output_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "stress_state_policy.h"

#include <numeric>

namespace Kratos
{
Element::Pointer SmallStrainUPwDiffOrderElement::Create(IndexType               NewId,
                                                        NodesArrayType const&   ThisNodes,
                                                        PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer SmallStrainUPwDiffOrderElement::Create(IndexType               NewId,
                                                        GeometryType::Pointer   pGeom,
                                                        PropertiesType::Pointer pProperties) const
{
    return make_intrusive<SmallStrainUPwDiffOrderElement>(NewId, pGeom, pProperties,
                                                          this->GetStressStatePolicy().Clone(),
                                                          this->CloneIntegrationCoefficientModifier());
}

int SmallStrainUPwDiffOrderElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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
    if (!r_prop[IGNORE_UNDRAINED])
        check_properties.CheckPermeabilityProperties(r_geom.WorkingSpaceDimension());

    check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW);
    r_prop[CONSTITUTIVE_LAW]->Check(r_prop, r_geom, rCurrentProcessInfo);
    const auto expected_size = this->GetStressStatePolicy().GetVoigtSize();
    ConstitutiveLawUtilities::CheckStrainSize(r_prop, expected_size, element_Id);
    ConstitutiveLawUtilities::CheckHasStrainMeasure_Infinitesimal(r_prop, element_Id);

    return RetentionLaw::Check(mRetentionLawVector, r_prop, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
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
        r_geom.WorkingSpaceDimension(), r_geom.PointsNumber(), integration_points.size(),
        r_geom.ShapeFunctionsValues(integration_method), solid_densities, integration_coefficients);

    rMassMatrix = ZeroMatrix(GetNumberOfDOF(), GetNumberOfDOF());
    GeoElementUtilities::AssembleUUBlockMatrix(rMassMatrix, mass_matrix_u);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(deformation_gradients);
    const auto strain_vectors = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        GetStressStatePolicy().GetVoigtSize());

    const auto number_of_integration_points = GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
    for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
        this->ExtractShapeFunctionDataAtIntegrationPoint(Variables, GPoint);
        Variables.B            = b_matrices[GPoint];
        Variables.F            = deformation_gradients[GPoint];
        Variables.StrainVector = strain_vectors[GPoint];

        ConstitutiveLawUtilities::SetConstitutiveParameters(
            ConstitutiveParameters, Variables.StrainVector, Variables.ConstitutiveMatrix, Variables.Nu,
            Variables.DNu_DX, Variables.F, determinants_of_deformation_gradients[GPoint]);

        // Compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
        mStateVariablesFinalized[GPoint] =
            mConstitutiveLawVector[GPoint]->GetValue(STATE_VARIABLES, mStateVariablesFinalized[GPoint]);
    }

    // Assign pressure values to the intermediate nodes for post-processing
    if (!GetProperties()[IGNORE_UNDRAINED]) AssignPressureToIntermediateNodes();

    KRATOS_CATCH("")
}

Vector SmallStrainUPwDiffOrderElement::GetPressures(const size_t n_nodes) const
{
    const auto& r_geom = GetGeometry();
    Vector      pressure(n_nodes);
    std::transform(r_geom.begin(), r_geom.begin() + n_nodes, pressure.begin(),
                   [](const auto& node) { return node.FastGetSolutionStepValue(WATER_PRESSURE); });
    return pressure;
}

void set_arithmetic_average_pressure(Geometry<Node>&                               rGeometry,
                                     const Vector&                                 rPressure,
                                     const std::vector<std::pair<size_t, size_t>>& rIndexPpairs,
                                     size_t DestinationOffset = 0)
{
    for (size_t i = 0; const auto& [first_index, second_index] : rIndexPpairs) {
        NodeUtilities::ThreadSafeNodeWrite(rGeometry[DestinationOffset + i], WATER_PRESSURE,
                                           0.5 * (rPressure[first_index] + rPressure[second_index]));
        ++i;
    }
}

void set_arithmetic_average_pressure(Geometry<Node>& rGeometry,
                                     const Vector&   rPressure,
                                     const std::vector<std::tuple<size_t, size_t, size_t, size_t>>& rIndices,
                                     size_t DestinationOffset = 0)
{
    for (size_t i = 0; const auto& [first_index, second_index, third_index, fourth_index] : rIndices) {
        NodeUtilities::ThreadSafeNodeWrite(rGeometry[DestinationOffset + i], WATER_PRESSURE,
                                           0.25 * (rPressure[first_index] + rPressure[second_index] +
                                                   rPressure[third_index] + rPressure[fourth_index]));
        ++i;
    }
}

void SmallStrainUPwDiffOrderElement::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType&  r_geom      = GetGeometry();
    const SizeType num_u_nodes = r_geom.PointsNumber();
    const SizeType n_dim       = r_geom.WorkingSpaceDimension();

    switch (num_u_nodes) {
    case 6: // 2D T6P3
    {
        const Vector                                 pressure = GetPressures(3);
        const std::vector<std::pair<size_t, size_t>> pairs    = {{0, 1}, {1, 2}, {2, 0}};
        set_arithmetic_average_pressure(r_geom, pressure, pairs, 3);
        break;
    }
    case 8: // 2D Q8P4
    {
        const Vector                                 pressure = GetPressures(4);
        const std::vector<std::pair<size_t, size_t>> pairs    = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        set_arithmetic_average_pressure(r_geom, pressure, pairs, 4);
        break;
    }
    case 9: // 2D Q9P4
    {
        const Vector                                 pressure = GetPressures(4);
        const std::vector<std::pair<size_t, size_t>> pairs    = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        set_arithmetic_average_pressure(r_geom, pressure, pairs, 4);
        const std::vector<std::tuple<size_t, size_t, size_t, size_t>>& indices = {{0, 1, 2, 3}};
        set_arithmetic_average_pressure(r_geom, pressure, indices, 8);
        break;
    }
    case 10: // 3D T10P4  //2D T10P6
    {
        if (n_dim == 3) {
            const Vector                                 pressure = GetPressures(4);
            const std::vector<std::pair<size_t, size_t>> pairs    = {{0, 1}, {1, 2}, {2, 0},
                                                                     {0, 3}, {1, 3}, {2, 3}};
            set_arithmetic_average_pressure(r_geom, pressure, pairs, 4);
        } else if (n_dim == 2) {
            constexpr double c1 = 1.0 / 9.0;
            const Vector     p  = GetPressures(6);
            NodeUtilities::ThreadSafeNodeWrite(r_geom[3], WATER_PRESSURE,
                                               (2.0 * p[0] - p[1] + 8.0 * p[3]) * c1);
            NodeUtilities::ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE,
                                               (2.0 * p[1] - p[0] + 8.0 * p[3]) * c1);
            NodeUtilities::ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE,
                                               (2.0 * p[1] - p[2] + 8.0 * p[4]) * c1);
            NodeUtilities::ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE,
                                               (2.0 * p[2] - p[1] + 8.0 * p[4]) * c1);
            NodeUtilities::ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE,
                                               (2.0 * p[2] - p[0] + 8.0 * p[5]) * c1);
            NodeUtilities::ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE,
                                               (2.0 * p[0] - p[2] + 8.0 * p[5]) * c1);
            NodeUtilities::ThreadSafeNodeWrite(
                r_geom[9], WATER_PRESSURE, (4.0 * (p[3] + p[4] + p[5]) - (p[0] + p[1] + p[2])) * c1);
        }
        break;
    }
    case 15: // 2D T15P10
    {
        constexpr double c1 = 0.0390625;
        const Vector     p  = GetPressures(10);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[3], WATER_PRESSURE,
                                           (3.0 * p[0] + p[1] + 27.0 * p[3] - 5.4 * p[4]) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE,
                                           (14.4 * (p[3] + p[4]) - 1.6 * (p[0] + p[1])) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE,
                                           (3.0 * p[1] + p[0] + 27.0 * p[4] - 5.4 * p[3]) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE,
                                           (3.0 * p[1] + p[2] + 27.0 * p[5] - 5.4 * p[6]) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE,
                                           (14.4 * (p[5] + p[6]) - 1.6 * (p[1] + p[2])) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE,
                                           (3.0 * p[2] + p[1] + 27.0 * p[6] - 5.4 * p[5]) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[9], WATER_PRESSURE,
                                           (3.0 * p[2] + p[0] + 27.0 * p[7] - 5.4 * p[8]) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[10], WATER_PRESSURE,
                                           (14.4 * (p[7] + p[8]) - 1.6 * (p[0] + p[2])) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[11], WATER_PRESSURE,
                                           (3.0 * p[0] + p[2] + 27.0 * p[8] - 5.4 * p[7]) * c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[12], WATER_PRESSURE,
                                           (p[1] + p[2] + 7.2 * (p[3] + p[8]) - 3.6 * (p[4] + p[7]) -
                                            1.8 * (p[5] + p[6]) + 21.6 * p[9] - 1.6 * p[0]) *
                                               c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[13], WATER_PRESSURE,
                                           (p[0] + p[2] + 7.2 * (p[4] + p[5]) - 3.6 * (p[3] + p[6]) -
                                            1.8 * (p[7] + p[8]) + 21.6 * p[9] - 1.6 * p[1]) *
                                               c1);
        NodeUtilities::ThreadSafeNodeWrite(r_geom[14], WATER_PRESSURE,
                                           (p[0] + p[1] + 7.2 * (p[6] + p[7]) - 3.6 * (p[5] + p[8]) -
                                            1.8 * (p[3] + p[4]) + 21.6 * p[9] - 1.6 * p[2]) *
                                               c1);
        break;
    }
    case 20: // 3D H20P8
    {
        const Vector                                 pressure = GetPressures(8);
        const std::vector<std::pair<size_t, size_t>> pairs =
            // edges -- bottom
            {{0, 1},
             {1, 2},
             {2, 3},
             {3, 0},
             // edges -- middle
             {4, 0},
             {5, 1},
             {6, 2},
             {7, 3},
             // edges -- top
             {4, 5},
             {5, 6},
             {6, 7},
             {7, 4}};
        set_arithmetic_average_pressure(r_geom, pressure, pairs, 8);
        break;
    }
    case 27: // 3D H27P8
    {
        const Vector                                 pressure = GetPressures(8);
        const std::vector<std::pair<size_t, size_t>> pairs =
            // edges -- bottom
            {{0, 1},
             {1, 2},
             {2, 3},
             {3, 0},
             // edges -- middle
             {4, 0},
             {5, 1},
             {6, 2},
             {7, 3},
             // edges -- top
             {4, 5},
             {5, 6},
             {6, 7},
             {7, 0}};
        set_arithmetic_average_pressure(r_geom, pressure, pairs, 8);
        // face centers
        const std::vector<std::tuple<size_t, size_t, size_t, size_t>>& indices = {
            {0, 1, 2, 3}, {0, 1, 4, 5}, {1, 2, 5, 6}, {2, 3, 6, 7}, {3, 0, 7, 4}, {4, 5, 6, 7}};
        set_arithmetic_average_pressure(r_geom, pressure, indices, 20);
        // element center
        NodeUtilities::ThreadSafeNodeWrite(r_geom[26], WATER_PRESSURE,
                                           0.125 * (pressure[0] + pressure[1] + pressure[2] + pressure[3] +
                                                    pressure[4] + pressure[5] + pressure[6] + pressure[7]));
        break;
    }
    default:
        KRATOS_ERROR << "Unexpected geometry type for different order "
                        "interpolation element"
                     << this->Id() << std::endl;
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                  const std::vector<Vector>& rValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
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
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            mConstitutiveLawVector[GPoint]->SetValue(rVariable, rValues[GPoint], rCurrentProcessInfo);
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
                                                                  std::vector<int>&    rValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
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

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                  std::vector<double>&    rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geom       = GetGeometry();
    const auto& r_properties = this->GetProperties();
    const auto number_of_integration_points = r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());

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

        RetentionLaw::Parameters RetentionParameters(GetProperties());

        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            // Compute Np, GradNpT, B and StrainVector
            this->ExtractShapeFunctionDataAtIntegrationPoint(Variables, GPoint);

            RetentionParameters.SetFluidPressure(GeoTransportEquationUtilities::CalculateFluidPressure(
                Variables.Np, Variables.PressureVector));

            rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateValue(
                RetentionParameters, rVariable, rOutput[GPoint]);
        }
    } else if (rVariable == HYDRAULIC_HEAD) {
        constexpr auto numerical_limit = std::numeric_limits<double>::epsilon();
        const auto&    r_prop          = this->GetProperties();

        // Defining the shape functions, the Jacobian and the shape functions local gradients Containers
        const Matrix&  n_container = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());
        const SizeType num_u_nodes = r_geom.PointsNumber();

        // Defining necessary variables
        Vector nodal_hydraulic_head = ZeroVector(num_u_nodes);
        for (unsigned int node = 0; node < num_u_nodes; ++node) {
            Vector NodeVolumeAcceleration(3);
            noalias(NodeVolumeAcceleration) = r_geom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
            const double g = norm_2(NodeVolumeAcceleration);
            if (g > numerical_limit) {
                const auto fluid_weight = g * r_prop[DENSITY_WATER];

                Vector node_coordinates(3);
                noalias(node_coordinates) = r_geom[node].Coordinates();
                Vector node_volume_acceleration_unit_vector(3);
                noalias(node_volume_acceleration_unit_vector) = NodeVolumeAcceleration / g;

                const auto water_pressure = r_geom[node].FastGetSolutionStepValue(WATER_PRESSURE);
                nodal_hydraulic_head[node] =
                    -inner_prod(node_coordinates, node_volume_acceleration_unit_vector) -
                    PORE_PRESSURE_SIGN_FACTOR * water_pressure / fluid_weight;
            } else {
                nodal_hydraulic_head[node] = 0.0;
            }
        }

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            const auto& shape_function_values = row(n_container, integration_point);
            rOutput[integration_point] =
                std::inner_product(shape_function_values.begin(), shape_function_values.end(),
                                   nodal_hydraulic_head.begin(), 0.0);
        }
    } else if (rVariable == CONFINED_STIFFNESS || rVariable == SHEAR_STIFFNESS) {
        KRATOS_ERROR_IF(r_geom.WorkingSpaceDimension() != 2 && r_geom.WorkingSpaceDimension() != 3)
            << rVariable.Name() << " can not be retrieved for dim "
            << r_geom.WorkingSpaceDimension() << " in element: " << this->Id() << std::endl;
        size_t variable_index = 0;
        if (rVariable == CONFINED_STIFFNESS) {
            variable_index = r_geom.WorkingSpaceDimension() == 2 ? static_cast<size_t>(INDEX_2D_PLANE_STRAIN_XX)
                                                                 : static_cast<size_t>(INDEX_3D_XX);
        } else {
            variable_index = r_geom.WorkingSpaceDimension() == 2 ? static_cast<size_t>(INDEX_2D_PLANE_STRAIN_XY)
                                                                 : static_cast<size_t>(INDEX_3D_XZ);
        }

        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, this->GetStressStatePolicy().GetVoigtSize());

        ConstitutiveLaw::Parameters ConstitutiveParameters(r_geom, r_properties, rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        std::vector<Matrix> constitutive_matrices;
        this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                             Variables.NuContainer, Variables.DNu_DXContainer,
                                             strain_vectors, mStressVector, constitutive_matrices);

        std::transform(constitutive_matrices.begin(), constitutive_matrices.end(), rOutput.begin(),
                       [variable_index](const Matrix& constitutive_matrix) {
            return constitutive_matrix(variable_index, variable_index);
        });
    } else if (r_properties.Has(rVariable)) {
        // Map initial material property to gauss points, as required for the output
        std::fill_n(rOutput.begin(), number_of_integration_points, r_properties.GetValue(rVariable));
    } else if (rVariable == GEO_SHEAR_CAPACITY) {
        OutputUtilities::CalculateShearCapacityValues(mStressVector, rOutput.begin(), GetProperties());
    } else {
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            rOutput[integration_point] = mConstitutiveLawVector[integration_point]->GetValue(
                rVariable, rOutput[integration_point]);
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                  std::vector<array_1d<double, 3>>& rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const auto number_of_integration_points = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    rOutput.resize(number_of_integration_points);

    if (rVariable == FLUID_FLUX_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        const auto strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
        auto relative_permeability_values = RetentionLaw::CalculateRelativePermeabilityValues(
            mRetentionLawVector, this->GetProperties(),
            GeoTransportEquationUtilities::CalculateFluidPressures(Variables.NpContainer, Variables.PressureVector));
        const auto permeability_update_factors =
            GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(strain_vectors, GetProperties());
        std::transform(relative_permeability_values.cbegin(), relative_permeability_values.cend(),
                       permeability_update_factors.cbegin(), relative_permeability_values.begin(),
                       std::multiplies<>{});

        // Loop over integration points
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        for (unsigned int g_point = 0; g_point < mConstitutiveLawVector.size(); ++g_point) {
            // compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->ExtractShapeFunctionDataAtIntegrationPoint(Variables, g_point);
            Variables.B = b_matrices[g_point];

            Vector   body_acceleration = ZeroVector(dimension);
            SizeType Index             = 0;
            for (SizeType i = 0; i < r_geometry.PointsNumber(); ++i) {
                for (unsigned int idim = 0; idim < dimension; ++idim) {
                    body_acceleration[idim] += Variables.Nu[i] * Variables.BodyAcceleration[Index];
                    ++Index;
                }
            }

            const auto relative_permeability = relative_permeability_values[g_point];

            Vector grad_pressure_term(dimension);
            noalias(grad_pressure_term) = prod(trans(Variables.DNp_DX), Variables.PressureVector);
            noalias(grad_pressure_term) +=
                PORE_PRESSURE_SIGN_FACTOR * GetProperties()[DENSITY_WATER] * body_acceleration;

            // Compute fluid flux vector q [L/T]
            rOutput[g_point].clear();
            const auto fluid_flux = PORE_PRESSURE_SIGN_FACTOR * Variables.DynamicViscosityInverse *
                                    relative_permeability *
                                    prod(Variables.IntrinsicPermeability, grad_pressure_term);
            std::copy_n(fluid_flux.begin(), dimension, rOutput[g_point].begin());
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                  std::vector<Vector>&    rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geom = GetGeometry();
    const auto number_of_integration_points = r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());
    rOutput.resize(number_of_integration_points);

    if (rVariable == CAUCHY_STRESS_VECTOR) {
        for (unsigned int GPoint = 0; GPoint < number_of_integration_points; ++GPoint) {
            if (rOutput[GPoint].size() != mStressVector[GPoint].size())
                rOutput[GPoint].resize(mStressVector[GPoint].size(), false);

            rOutput[GPoint] = mStressVector[GPoint];
        }
    } else if (rVariable == TOTAL_STRESS_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        ConstitutiveLaw::Parameters ConstitutiveParameters(r_geom, GetProperties(), rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
        std::vector<Matrix> constitutive_matrices;
        this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                             Variables.NuContainer, Variables.DNu_DXContainer,
                                             strain_vectors, mStressVector, constitutive_matrices);
        const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
            constitutive_matrices, this->GetProperties());
        const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
            Variables.NpContainer, Variables.PressureVector);
        const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] =
                mStressVector[GPoint] + PORE_PRESSURE_SIGN_FACTOR * biot_coefficients[GPoint] *
                                            bishop_coefficients[GPoint] * fluid_pressures[GPoint] *
                                            GetStressStatePolicy().GetVoigtVector();
        }
    } else if (rVariable == ENGINEERING_STRAIN_VECTOR) {
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            noalias(Variables.Nu) = row(Variables.NuContainer, GPoint);

            Matrix J0;
            Matrix InvJ0;
            double detJInitialConfiguration;
            Matrix DNu_DXInitialConfiguration;
            this->CalculateDerivativesOnInitialConfiguration(detJInitialConfiguration, J0, InvJ0,
                                                             DNu_DXInitialConfiguration, GPoint);

            // Calculating operator B
            Variables.B = this->CalculateBMatrix(DNu_DXInitialConfiguration, Variables.Nu);

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

        const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
        const auto deformation_gradients = CalculateDeformationGradients();
        rOutput                          = StressStrainUtilities::CalculateStrains(
            deformation_gradients, b_matrices, Variables.DisplacementVector,
            Variables.UseHenckyStrain, GetStressStatePolicy().GetVoigtSize());
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                  std::vector<Matrix>&    rOutput,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rOutput.resize(GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()));

    if (rVariable == CAUCHY_STRESS_TENSOR) {
        std::vector<Vector> StressVector;
        this->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
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
    } else {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            rOutput[i] = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::Calculate(const Variable<Vector>& rVariable,
                                               Vector&                 rOutput,
                                               const ProcessInfo&      rCurrentProcessInfo)
{
    KRATOS_ERROR_IF_NOT(rVariable == INTERNAL_FORCES_VECTOR || rVariable == EXTERNAL_FORCES_VECTOR)
        << "Variable " << rVariable.Name() << " is unknown for element with Id " << this->GetId() << ".";

    rOutput = Vector(this->GetNumberOfDOF(), 0.0);

    const PropertiesType&                           r_prop = this->GetProperties();
    const GeometryType&                             r_geom = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geom, r_prop, rCurrentProcessInfo);

    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(r_integration_points, Variables.detJuContainer);

    const auto det_Js_initial_configuration = GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(
        r_geom, this->GetIntegrationMethod());

    const auto deformation_gradients = CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NuContainer, Variables.DNu_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);
    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Variables.NpContainer, Variables.PressureVector);
    const auto degrees_of_saturation = CalculateDegreesOfSaturation(fluid_pressures);

    auto relative_permeability_values = RetentionLaw::CalculateRelativePermeabilityValues(
        mRetentionLawVector, this->GetProperties(), fluid_pressures);
    const auto permeability_update_factors = GetOptionalPermeabilityUpdateFactors(strain_vectors);
    std::ranges::transform(permeability_update_factors, relative_permeability_values,
                           relative_permeability_values.begin(), std::multiplies<>{});
    const auto bishop_coefficients = CalculateBishopCoefficients(fluid_pressures);

    if (rVariable == INTERNAL_FORCES_VECTOR) {
        const auto derivatives_of_saturation = CalculateDerivativesOfSaturation(fluid_pressures);
        const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
            constitutive_matrices, this->GetProperties());
        const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
            biot_coefficients, degrees_of_saturation, derivatives_of_saturation, r_prop);
        rOutput = CalculateInternalForces(Variables, b_matrices, integration_coefficients,
                                          biot_coefficients, degrees_of_saturation, biot_moduli_inverse,
                                          relative_permeability_values, bishop_coefficients);
    } else if (rVariable == EXTERNAL_FORCES_VECTOR) {
        const auto integration_coefficients_on_initial_configuration =
            this->CalculateIntegrationCoefficients(r_integration_points, det_Js_initial_configuration);
        rOutput = CalculateExternalForces(
            Variables, integration_coefficients, integration_coefficients_on_initial_configuration,
            degrees_of_saturation, relative_permeability_values, bishop_coefficients);
    }
}

Vector SmallStrainUPwDiffOrderElement::CalculateInternalForces(ElementVariables& rVariables,
                                                               const std::vector<Matrix>& rBMatrices,
                                                               const std::vector<double>& rIntegrationCoefficients,
                                                               const std::vector<double>& rBiotCoefficients,
                                                               const std::vector<double>& rDegreesOfSaturation,
                                                               const std::vector<double>& rBiotModuliInverse,
                                                               const std::vector<double>& rRelativePermeabilityValues,
                                                               const std::vector<double>& rBishopCoefficients) const
{
    Vector result(this->GetNumberOfDOF(), 0.0);
    for (unsigned int integration_point = 0; integration_point < rIntegrationCoefficients.size(); ++integration_point) {
        rVariables.B                      = rBMatrices[integration_point];
        rVariables.IntegrationCoefficient = rIntegrationCoefficients[integration_point];

        this->CalculateAndAddStiffnessForce(result, rVariables, integration_point);
    }

    for (unsigned int integration_point = 0; integration_point < rIntegrationCoefficients.size(); ++integration_point) {
        rVariables.B                      = rBMatrices[integration_point];
        rVariables.BishopCoefficient      = rBishopCoefficients[integration_point];
        rVariables.BiotCoefficient        = rBiotCoefficients[integration_point];
        rVariables.DegreeOfSaturation     = rDegreesOfSaturation[integration_point];
        rVariables.IntegrationCoefficient = rIntegrationCoefficients[integration_point];
        noalias(rVariables.Np)            = row(rVariables.NpContainer, integration_point);

        this->CalculateAndAddCouplingTerms(result, rVariables);
    }
    if (!rVariables.IgnoreUndrained) {
        for (unsigned int integration_point = 0;
             integration_point < rIntegrationCoefficients.size(); ++integration_point) {
            noalias(rVariables.Np)            = row(rVariables.NpContainer, integration_point);
            rVariables.BiotModulusInverse     = rBiotModuliInverse[integration_point];
            rVariables.IntegrationCoefficient = rIntegrationCoefficients[integration_point];

            CalculateAndAddCompressibilityFlow(result, rVariables);
        }
        for (unsigned int integration_point = 0;
             integration_point < rIntegrationCoefficients.size(); ++integration_point) {
            noalias(rVariables.DNp_DX)        = rVariables.DNp_DXContainer[integration_point];
            rVariables.RelativePermeability   = rRelativePermeabilityValues[integration_point];
            rVariables.IntegrationCoefficient = rIntegrationCoefficients[integration_point];

            CalculateAndAddPermeabilityFlow(result, rVariables);
        }
    }

    return result;
}

Vector SmallStrainUPwDiffOrderElement::CalculateExternalForces(
    ElementVariables&          rVariables,
    const std::vector<double>& rIntegrationCoefficients,
    const std::vector<double>& rIntegrationCoefficientsOnInitialConfiguration,
    const std::vector<double>& rDegreesOfSaturation,
    const std::vector<double>& rRelativePermeabilityValues,
    const std::vector<double>& rBishopCoefficients) const
{
    Vector result = ZeroVector(this->GetNumberOfDOF());
    for (unsigned int integration_point = 0; integration_point < rIntegrationCoefficients.size(); ++integration_point) {
        noalias(rVariables.Nu)        = row(rVariables.NuContainer, integration_point);
        rVariables.DegreeOfSaturation = rDegreesOfSaturation[integration_point];
        rVariables.IntegrationCoefficientInitialConfiguration =
            rIntegrationCoefficientsOnInitialConfiguration[integration_point];
        this->CalculateAndAddMixBodyForce(result, rVariables);
    }
    if (!rVariables.IgnoreUndrained) {
        for (unsigned int integration_point = 0;
             integration_point < rIntegrationCoefficients.size(); ++integration_point) {
            noalias(rVariables.Nu)            = row(rVariables.NuContainer, integration_point);
            noalias(rVariables.DNp_DX)        = rVariables.DNp_DXContainer[integration_point];
            rVariables.RelativePermeability   = rRelativePermeabilityValues[integration_point];
            rVariables.BishopCoefficient      = rBishopCoefficients[integration_point];
            rVariables.IntegrationCoefficient = rIntegrationCoefficients[integration_point];

            this->CalculateAndAddFluidBodyFlow(result, rVariables);
        }
    }

    return result;
}

void SmallStrainUPwDiffOrderElement::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                  VectorType&        rRightHandSideVector,
                                                  const ProcessInfo& rCurrentProcessInfo,
                                                  bool               CalculateStiffnessMatrixFlag,
                                                  bool               CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const PropertiesType&                           r_prop = this->GetProperties();
    const GeometryType&                             r_geom = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geom, r_prop, rCurrentProcessInfo);

    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(r_integration_points, Variables.detJuContainer);

    const auto det_Js_initial_configuration = GeoEquationOfMotionUtilities::CalculateDetJsInitialConfiguration(
        r_geom, this->GetIntegrationMethod());

    const auto integration_coefficients_on_initial_configuration =
        this->CalculateIntegrationCoefficients(r_integration_points, det_Js_initial_configuration);

    const auto deformation_gradients = CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NuContainer, Variables.DNu_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);
    const auto biot_coefficients = GeoTransportEquationUtilities::CalculateBiotCoefficients(
        constitutive_matrices, this->GetProperties());
    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Variables.NpContainer, Variables.PressureVector);
    const auto degrees_of_saturation     = CalculateDegreesOfSaturation(fluid_pressures);
    const auto derivatives_of_saturation = CalculateDerivativesOfSaturation(fluid_pressures);
    const auto biot_moduli_inverse = GeoTransportEquationUtilities::CalculateInverseBiotModuli(
        biot_coefficients, degrees_of_saturation, derivatives_of_saturation, r_prop);
    auto relative_permeability_values = RetentionLaw::CalculateRelativePermeabilityValues(
        mRetentionLawVector, this->GetProperties(), fluid_pressures);
    const auto permeability_update_factors = GetOptionalPermeabilityUpdateFactors(strain_vectors);
    std::ranges::transform(permeability_update_factors, relative_permeability_values,
                           relative_permeability_values.begin(), std::multiplies<>{});

    const auto bishop_coefficients = CalculateBishopCoefficients(fluid_pressures);

    if (CalculateStiffnessMatrixFlag) {
        for (unsigned int GPoint = 0; GPoint < r_integration_points.size(); ++GPoint) {
            this->ExtractShapeFunctionDataAtIntegrationPoint(Variables, GPoint);
            Variables.B                  = b_matrices[GPoint];
            Variables.F                  = deformation_gradients[GPoint];
            Variables.StrainVector       = strain_vectors[GPoint];
            Variables.ConstitutiveMatrix = constitutive_matrices[GPoint];

            Variables.RelativePermeability = relative_permeability_values[GPoint];
            Variables.BishopCoefficient    = bishop_coefficients[GPoint];

            Variables.BiotCoefficient        = biot_coefficients[GPoint];
            Variables.BiotModulusInverse     = biot_moduli_inverse[GPoint];
            Variables.DegreeOfSaturation     = degrees_of_saturation[GPoint];
            Variables.IntegrationCoefficient = integration_coefficients[GPoint];

            Variables.IntegrationCoefficientInitialConfiguration =
                integration_coefficients_on_initial_configuration[GPoint];

            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
        }
    }

    if (CalculateResidualVectorFlag) {
        const auto internal_forces = CalculateInternalForces(
            Variables, b_matrices, integration_coefficients, biot_coefficients, degrees_of_saturation,
            biot_moduli_inverse, relative_permeability_values, bishop_coefficients);

        const auto external_forces = CalculateExternalForces(
            Variables, integration_coefficients, integration_coefficients_on_initial_configuration,
            degrees_of_saturation, relative_permeability_values, bishop_coefficients);
        rRightHandSideVector = external_forces - internal_forces;
    }
    KRATOS_CATCH("")
}

std::vector<double> SmallStrainUPwDiffOrderElement::GetOptionalPermeabilityUpdateFactors(const std::vector<Vector>& rStrainVectors) const
{
    return GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(rStrainVectors, GetProperties());
}

std::vector<double> SmallStrainUPwDiffOrderElement::CalculateDerivativesOfSaturation(const std::vector<double>& rFluidPressures)
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

std::vector<double> SmallStrainUPwDiffOrderElement::CalculateDegreesOfSaturation(const std::vector<double>& rFluidPressures)
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

void SmallStrainUPwDiffOrderElement::CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix,
                                                                      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geom = GetGeometry();

    // Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(r_geom, GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Loop over integration points
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());

    const auto b_matrices = CalculateBMatrices(Variables.DNu_DXContainer, Variables.NuContainer);
    const auto deformation_gradients = CalculateDeformationGradients();
    auto       strain_vectors        = StressStrainUtilities::CalculateStrains(
        deformation_gradients, b_matrices, Variables.DisplacementVector, Variables.UseHenckyStrain,
        this->GetStressStatePolicy().GetVoigtSize());
    std::vector<Matrix> constitutive_matrices;
    this->CalculateAnyOfMaterialResponse(deformation_gradients, ConstitutiveParameters,
                                         Variables.NuContainer, Variables.DNu_DXContainer,
                                         strain_vectors, mStressVector, constitutive_matrices);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(r_integration_points, Variables.detJuContainer);

    const auto stiffness_matrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        b_matrices, constitutive_matrices, integration_coefficients);

    GeoElementUtilities::AssembleUUBlockMatrix(rStiffnessMatrix, stiffness_matrix);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::InitializeElementVariables(ElementVariables& rVariables,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geom       = GetGeometry();
    const SizeType      num_u_nodes  = r_geom.PointsNumber();
    const SizeType      num_p_nodes  = mpPressureGeometry->PointsNumber();
    const SizeType      num_g_points = r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());
    const SizeType      n_dim        = r_geom.WorkingSpaceDimension();

    // Variables at all integration points
    rVariables.NuContainer.resize(num_g_points, num_u_nodes, false);
    rVariables.NuContainer = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());

    rVariables.NpContainer.resize(num_g_points, num_p_nodes, false);
    rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues(this->GetIntegrationMethod());

    rVariables.Nu.resize(num_u_nodes, false);
    rVariables.Np.resize(num_p_nodes, false);

    rVariables.DNu_DXContainer.resize(num_g_points, false);
    for (SizeType i = 0; i < num_g_points; ++i)
        ((rVariables.DNu_DXContainer)[i]).resize(num_u_nodes, n_dim, false);
    rVariables.DNu_DX.resize(num_u_nodes, n_dim, false);
    rVariables.detJuContainer.resize(num_g_points, false);

    try {
        r_geom.ShapeFunctionsIntegrationPointsGradients(
            rVariables.DNu_DXContainer, rVariables.detJuContainer, this->GetIntegrationMethod());
    } catch (Kratos::Exception& e) {
        KRATOS_INFO("Original error message") << e.what() << std::endl;
#ifdef KRATOS_COMPILED_IN_WINDOWS
        KRATOS_INFO("Error in calculation of dNu/dx. Most probably the element is "
                    "distorted. Element ID: ")
            << this->Id() << std::endl;
#endif
        KRATOS_ERROR << "In calculation of dNu/dx. Most probably the element "
                        "is distorted. Element ID: "
                     << this->Id() << std::endl;
    }

    (rVariables.DNp_DXContainer).resize(num_g_points, false);
    for (SizeType i = 0; i < num_g_points; ++i)
        ((rVariables.DNp_DXContainer)[i]).resize(num_p_nodes, n_dim, false);
    (rVariables.DNp_DX).resize(num_p_nodes, n_dim, false);
    Vector detJpContainer = ZeroVector(num_g_points);

    try {
        mpPressureGeometry->ShapeFunctionsIntegrationPointsGradients(
            rVariables.DNp_DXContainer, detJpContainer, this->GetIntegrationMethod());
    } catch (Kratos::Exception& e) {
        KRATOS_INFO("Original error message") << e.what() << std::endl;
#ifdef KRATOS_COMPILED_IN_WINDOWS
        KRATOS_INFO("Error in calculation of dNp/dx. Most probably the element is "
                    "distorted. Element ID: ")
            << this->Id() << std::endl;
#endif
        KRATOS_ERROR << "In calculation of dNp/dx. Most probably the element "
                        "is distorted. Element ID: "
                     << this->Id() << std::endl;
    }

    // Variables computed at each integration point
    const SizeType VoigtSize = this->GetStressStatePolicy().GetVoigtSize();

    rVariables.B.resize(VoigtSize, num_u_nodes * n_dim, false);
    noalias(rVariables.B) = ZeroMatrix(VoigtSize, num_u_nodes * n_dim);

    rVariables.StrainVector.resize(VoigtSize, false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);

    rVariables.StressVector.resize(VoigtSize, false);

    // Needed parameters for consistency with the general constitutive law
    rVariables.F.resize(n_dim, n_dim, false);
    noalias(rVariables.F) = identity_matrix<double>(n_dim);

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

void SmallStrainUPwDiffOrderElement::InitializeNodalVariables(ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& r_geom      = GetGeometry();
    const SizeType      n_dim       = r_geom.WorkingSpaceDimension();
    const SizeType      num_u_nodes = r_geom.PointsNumber();
    const SizeType      num_p_nodes = mpPressureGeometry->PointsNumber();

    Vector BodyAccelerationAux = ZeroVector(3);
    rVariables.BodyAcceleration.resize(num_u_nodes * n_dim, false);
    rVariables.DisplacementVector.resize(num_u_nodes * n_dim, false);
    rVariables.VelocityVector.resize(num_u_nodes * n_dim, false);

    for (SizeType i = 0; i < num_u_nodes; ++i) {
        SizeType Local_i    = i * n_dim;
        BodyAccelerationAux = r_geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        rVariables.BodyAcceleration[Local_i]   = BodyAccelerationAux[0];
        rVariables.DisplacementVector[Local_i] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
        rVariables.VelocityVector[Local_i]     = r_geom[i].FastGetSolutionStepValue(VELOCITY_X);

        rVariables.BodyAcceleration[Local_i + 1] = BodyAccelerationAux[1];
        rVariables.DisplacementVector[Local_i + 1] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        rVariables.VelocityVector[Local_i + 1] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Y);

        if (n_dim > 2) {
            rVariables.BodyAcceleration[Local_i + 2] = BodyAccelerationAux[2];
            rVariables.DisplacementVector[Local_i + 2] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
            rVariables.VelocityVector[Local_i + 2] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Z);
        }
    }

    rVariables.PressureVector.resize(num_p_nodes, false);
    rVariables.PressureDtVector.resize(num_p_nodes, false);
    rVariables.DeltaPressureVector.resize(num_p_nodes, false);
    const auto& r_p_geometry = *mpPressureGeometry;
    for (SizeType i = 0; i < num_p_nodes; ++i) {
        rVariables.PressureVector[i] = r_p_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.PressureDtVector[i] = r_p_geometry[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
        rVariables.DeltaPressureVector[i] = r_p_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE) -
                                            r_p_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE, 1);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::InitializeProperties(ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto& r_properties = this->GetProperties();

    rVariables.IgnoreUndrained = r_properties[IGNORE_UNDRAINED];
    rVariables.UseHenckyStrain = r_properties.Has(USE_HENCKY_STRAIN) ? r_properties[USE_HENCKY_STRAIN] : false;

    rVariables.ConsiderGeometricStiffness =
        r_properties.Has(CONSIDER_GEOMETRIC_STIFFNESS) ? r_properties[CONSIDER_GEOMETRIC_STIFFNESS] : false;

    rVariables.DynamicViscosityInverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
    // Setting the intrinsic permeability matrix
    rVariables.IntrinsicPermeability =
        GeoElementUtilities::FillPermeabilityMatrix(r_properties, GetGeometry().WorkingSpaceDimension());

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::ExtractShapeFunctionDataAtIntegrationPoint(ElementVariables& rVariables,
                                                                                unsigned int GPoint)
{
    KRATOS_TRY

    noalias(rVariables.Nu) = row(rVariables.NuContainer, GPoint);
    noalias(rVariables.Np) = row(rVariables.NpContainer, GPoint);

    noalias(rVariables.DNu_DX) = rVariables.DNu_DXContainer[GPoint];
    noalias(rVariables.DNp_DX) = rVariables.DNp_DXContainer[GPoint];

    rVariables.detJ = rVariables.detJuContainer[GPoint];

    KRATOS_CATCH("")
}

Matrix SmallStrainUPwDiffOrderElement::CalculateBMatrix(const Matrix& rDN_DX, const Vector& rN) const
{
    return this->GetStressStatePolicy().CalculateBMatrix(rDN_DX, rN, this->GetGeometry());
}

std::vector<Matrix> SmallStrainUPwDiffOrderElement::CalculateBMatrices(
    const GeometryType::ShapeFunctionsGradientsType& rDN_DXContainer, const Matrix& rNContainer) const
{
    std::vector<Matrix> result;
    result.reserve(rDN_DXContainer.size());
    for (unsigned int GPoint = 0; GPoint < rDN_DXContainer.size(); ++GPoint) {
        result.push_back(this->CalculateBMatrix(rDN_DXContainer[GPoint], row(rNContainer, GPoint)));
    }

    return result;
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddLHS(MatrixType&             rLeftHandSideMatrix,
                                                        const ElementVariables& rVariables) const
{
    KRATOS_TRY

    CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);

    this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);

    if (!rVariables.IgnoreUndrained) {
        const auto permeability_matrix = GeoTransportEquationUtilities::CalculatePermeabilityMatrix(
            rVariables.DNp_DX, rVariables.DynamicViscosityInverse, rVariables.IntrinsicPermeability,
            rVariables.RelativePermeability, rVariables.IntegrationCoefficient);
        GeoElementUtilities::AssemblePPBlockMatrix(rLeftHandSideMatrix, permeability_matrix);

        CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix,
                                                                    const ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto number_of_dofs = rVariables.B.size2();
    Matrix     stiffness_matrix(number_of_dofs, number_of_dofs);

    GeoEquationOfMotionUtilities::CalculateStiffnessMatrixGPoint(
        stiffness_matrix, rVariables.B, rVariables.ConstitutiveMatrix, rVariables.IntegrationCoefficient);

    GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, stiffness_matrix);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix,
                                                                   const ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix coupling_matrix(this->GetGeometry().WorkingSpaceDimension() * this->GetGeometry().PointsNumber(),
                           mpPressureGeometry->PointsNumber(), 0.0);
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

void SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                                          const ElementVariables& rVariables)
{
    KRATOS_TRY

    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    GeoElementUtilities::AssemblePPBlockMatrix(
        rLeftHandSideMatrix, compressibility_matrix * rVariables.DtPressureCoefficient);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector,
                                                                   const ElementVariables& rVariables,
                                                                   unsigned int GPoint) const
{
    KRATOS_TRY

    Vector stiffness_force = prod(trans(rVariables.B), mStressVector[GPoint]) * rVariables.IntegrationCoefficient;
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, stiffness_force);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector,
                                                                 ElementVariables& rVariables) const
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      Dim       = rGeom.WorkingSpaceDimension();
    const SizeType      NumUNodes = rGeom.PointsNumber();

    const auto soil_density = GeoTransportEquationUtilities::CalculateSoilDensity(
        rVariables.DegreeOfSaturation, this->GetProperties());

    Vector   body_acceleration = ZeroVector(Dim);
    SizeType Index             = 0;
    for (SizeType i = 0; i < NumUNodes; ++i) {
        for (SizeType idim = 0; idim < Dim; ++idim) {
            body_acceleration[idim] += rVariables.Nu[i] * rVariables.BodyAcceleration[Index];
            ++Index;
        }
    }

    for (SizeType i = 0; i < NumUNodes; ++i) {
        Index = i * Dim;
        for (SizeType idim = 0; idim < Dim; ++idim) {
            rRightHandSideVector[Index + idim] += rVariables.Nu[i] * soil_density * body_acceleration[idim] *
                                                  rVariables.IntegrationCoefficientInitialConfiguration;
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector,
                                                                  const ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix coupling_matrix(this->GetGeometry().WorkingSpaceDimension() * this->GetGeometry().PointsNumber(),
                           mpPressureGeometry->PointsNumber(), 0.0);
    GeoTransportEquationUtilities::CalculateCouplingMatrix(
        coupling_matrix, rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
        rVariables.BiotCoefficient, rVariables.BishopCoefficient, rVariables.IntegrationCoefficient);
    const Vector coupling_force = prod(coupling_matrix, rVariables.PressureVector);
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, coupling_force);

    if (!rVariables.IgnoreUndrained) {
        GeoTransportEquationUtilities::CalculateCouplingMatrix(
            coupling_matrix, rVariables.B, GetStressStatePolicy().GetVoigtVector(), rVariables.Np,
            rVariables.BiotCoefficient, rVariables.DegreeOfSaturation, rVariables.IntegrationCoefficient);
        const Vector coupling_flow =
            PORE_PRESSURE_SIGN_FACTOR * prod(trans(coupling_matrix), rVariables.VelocityVector);
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, coupling_flow);
    }

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector,
                                                                        const ElementVariables& rVariables)
{
    KRATOS_TRY

    Matrix compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);
    Vector compressibility_flow = prod(compressibility_matrix, rVariables.PressureDtVector);
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, compressibility_flow);

    KRATOS_CATCH("")
}

std::vector<double> SmallStrainUPwDiffOrderElement::CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const
{
    KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

    auto result = std::vector<double>{};
    result.reserve(mRetentionLawVector.size());
    std::transform(mRetentionLawVector.begin(), mRetentionLawVector.end(), rFluidPressures.begin(),
                   std::back_inserter(result),
                   [&retention_law_params](const auto& pRetentionLaw, auto FluidPressure) {
        retention_law_params.SetFluidPressure(FluidPressure);
        return pRetentionLaw->CalculateBishopCoefficient(retention_law_params);
    });
    return result;
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                                                     const ElementVariables& rVariables)
{
    KRATOS_TRY

    const Matrix permeability_matrix =
        -PORE_PRESSURE_SIGN_FACTOR * rVariables.DynamicViscosityInverse * rVariables.RelativePermeability *
        prod(rVariables.DNp_DX, Matrix(prod(rVariables.IntrinsicPermeability, trans(rVariables.DNp_DX)))) *
        rVariables.IntegrationCoefficient;
    const Vector permeability_flow = prod(permeability_matrix, rVariables.PressureVector);
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, permeability_flow);

    KRATOS_CATCH("")
}

void SmallStrainUPwDiffOrderElement::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                                                  const ElementVariables& rVariables) const
{
    KRATOS_TRY

    const Matrix grad_Np_T_perm = rVariables.DynamicViscosityInverse * rVariables.BishopCoefficient *
                                  GetProperties()[DENSITY_WATER] * rVariables.RelativePermeability *
                                  prod(rVariables.DNp_DX, rVariables.IntrinsicPermeability) *
                                  rVariables.IntegrationCoefficient;

    const GeometryType& r_geom      = GetGeometry();
    const SizeType      dimension   = r_geom.WorkingSpaceDimension();
    const SizeType      num_U_nodes = r_geom.PointsNumber();

    Vector body_acceleration = ZeroVector(dimension);

    SizeType index = 0;
    for (SizeType i = 0; i < num_U_nodes; ++i) {
        for (SizeType idim = 0; idim < dimension; ++idim) {
            body_acceleration[idim] += rVariables.Nu[i] * rVariables.BodyAcceleration[index];
            index++;
        }
    }

    const Vector fluid_body_flow = prod(grad_Np_T_perm, body_acceleration);
    GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector, fluid_body_flow);

    KRATOS_CATCH("")
}

Vector SmallStrainUPwDiffOrderElement::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
{
    return this->GetStressStatePolicy().CalculateGreenLagrangeStrain(rDeformationGradient);
}

Matrix SmallStrainUPwDiffOrderElement::CalculateDeformationGradient(unsigned int GPoint) const
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

std::vector<Matrix> SmallStrainUPwDiffOrderElement::CalculateDeformationGradients() const
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

SizeType SmallStrainUPwDiffOrderElement::GetNumberOfDOF() const
{
    return GetGeometry().PointsNumber() * GetGeometry().WorkingSpaceDimension() +
           mpPressureGeometry->PointsNumber();
}

Element::DofsVectorType SmallStrainUPwDiffOrderElement::GetDofs() const
{
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), *mpPressureGeometry,
                                                      GetGeometry().WorkingSpaceDimension());
}

void SmallStrainUPwDiffOrderElement::SetUpPressureGeometryPointer()
{
    const auto& r_geometry        = GetGeometry();
    const auto  number_of_U_nodes = r_geometry.PointsNumber();
    const auto  dimension         = r_geometry.WorkingSpaceDimension();
    switch (number_of_U_nodes) {
    case 6: // 2D T6P3
        mpPressureGeometry = make_shared<Triangle2D3<Node>>(r_geometry(0), r_geometry(1), r_geometry(2));
        break;
    case 8: // 2D Q8P4
        mpPressureGeometry = make_shared<Quadrilateral2D4<Node>>(r_geometry(0), r_geometry(1),
                                                                 r_geometry(2), r_geometry(3));
        break;
    case 9: // 2D Q9P4
        mpPressureGeometry = make_shared<Quadrilateral2D4<Node>>(r_geometry(0), r_geometry(1),
                                                                 r_geometry(2), r_geometry(3));
        break;
    case 10:
        if (dimension == 3) // 3D T10P4
            mpPressureGeometry = make_shared<Tetrahedra3D4<Node>>(r_geometry(0), r_geometry(1),
                                                                  r_geometry(2), r_geometry(3));
        else if (dimension == 2) // 2D T10P6
            mpPressureGeometry = make_shared<Triangle2D6<Node>>(
                r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3), r_geometry(4), r_geometry(5));
        break;
    case 15: // 2D T15P10
        mpPressureGeometry = make_shared<Triangle2D10<Node>>(
            r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3), r_geometry(4),
            r_geometry(5), r_geometry(6), r_geometry(7), r_geometry(8), r_geometry(9));
        break;
    case 20: // 3D H20P8
        mpPressureGeometry =
            make_shared<Hexahedra3D8<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3),
                                            r_geometry(4), r_geometry(5), r_geometry(6), r_geometry(7));
        break;
    case 27: // 3D H27P8
        mpPressureGeometry =
            make_shared<Hexahedra3D8<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3),
                                            r_geometry(4), r_geometry(5), r_geometry(6), r_geometry(7));
        break;
    default:
        KRATOS_ERROR << "Unexpected geometry type for different order interpolation element "
                     << this->Id() << std::endl;
    }
}

void SmallStrainUPwDiffOrderElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, UPwBaseElement)
    rSerializer.save("PressureGeometry", mpPressureGeometry);
}

void SmallStrainUPwDiffOrderElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, UPwBaseElement)
    rSerializer.load("PressureGeometry", mpPressureGeometry);
}

Vector SmallStrainUPwDiffOrderElement::GetPressureSolutionVector() const
{
    Vector result(mpPressureGeometry->PointsNumber());
    std::transform(mpPressureGeometry->begin(), mpPressureGeometry->end(), result.begin(),
                   [](const auto& node) { return node.FastGetSolutionStepValue(WATER_PRESSURE); });
    return result;
}

void SmallStrainUPwDiffOrderElement::CalculateAnyOfMaterialResponse(
    const std::vector<Matrix>&                       rDeformationGradients,
    ConstitutiveLaw::Parameters&                     rConstitutiveParameters,
    const Matrix&                                    rNuContainer,
    const GeometryType::ShapeFunctionsGradientsType& rDNu_DXContainer,
    std::vector<Vector>&                             rStrainVectors,
    std::vector<Vector>&                             rStressVectors,
    std::vector<Matrix>&                             rConstitutiveMatrices)
{
    const SizeType voigt_size =
        GetGeometry().WorkingSpaceDimension() == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN;

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
        std::fill(rConstitutiveMatrices.begin(), rConstitutiveMatrices.end(), ZeroMatrix(voigt_size, voigt_size));
    }

    const auto determinants_of_deformation_gradients =
        GeoMechanicsMathUtilities::CalculateDeterminants(rDeformationGradients);

    for (unsigned int GPoint = 0; GPoint < rDeformationGradients.size(); ++GPoint) {
        // Explicitly convert from `row`'s return type to `Vector` to avoid ending up with a pointer
        // to an implicitly converted object
        const auto shape_function_values = Vector{row(rNuContainer, GPoint)};
        ConstitutiveLawUtilities::SetConstitutiveParameters(
            rConstitutiveParameters, rStrainVectors[GPoint], rConstitutiveMatrices[GPoint],
            shape_function_values, rDNu_DXContainer[GPoint], rDeformationGradients[GPoint],
            determinants_of_deformation_gradients[GPoint]);
        rConstitutiveParameters.SetStressVector(rStressVectors[GPoint]);

        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(rConstitutiveParameters);
    }
}

std::string SmallStrainUPwDiffOrderElement::Info() const
{
    const std::string constitutive_info =
        !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";
    return "U-Pw small strain different order Element #" + std::to_string(Id()) +
           "\nConstitutive law: " + constitutive_info;
}

void SmallStrainUPwDiffOrderElement::PrintInfo(std::ostream& rOStream) const { rOStream << Info(); }

} // Namespace Kratos
