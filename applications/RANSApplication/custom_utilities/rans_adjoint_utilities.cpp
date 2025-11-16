//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <tuple>

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/model_part.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_adjoint_utilities.h"

namespace Kratos
{

double RansAdjointUtilities::CalculateVectorNormDerivative(
    const double VectorNorm,
    const array_1d<double, 3>& rVector,
    const array_1d<double, 3>& rVectorDerivative)
{
    if (VectorNorm > 0.0) {
        return inner_prod(rVector, rVectorDerivative) / VectorNorm;
    } else {
        return 0.0;
    }
}

array_1d<double, 3> RansAdjointUtilities::CalculateUnitVectorDerivative(
    const double VectorMagnitude,
    const array_1d<double, 3>& rUnitVector,
    const array_1d<double, 3>& rVectorDerivative)
{
    if (VectorMagnitude > 0.0) {
        return (rVectorDerivative - rUnitVector * inner_prod(rUnitVector, rVectorDerivative)) /
               VectorMagnitude;
    } else {
        return ZeroVector(3);
    }
}

double RansAdjointUtilities::CalculateWallHeightConditionDerivative(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentElementGeometry,
    const IndexType DirectionIndex,
    const array_1d<double, 3>& rUnitNormal,
    const array_1d<double, 3>& rUnitNormalDerivative)
{
    const auto& condition_center = rConditionGeometry.Center();
    const auto& parent_center = rParentElementGeometry.Center();

    array_1d<double, 3> condition_center_derivative = ZeroVector(3);
    condition_center_derivative[DirectionIndex] = 1.0 / rConditionGeometry.PointsNumber();

    array_1d<double, 3> parent_center_derivative = ZeroVector(3);
    parent_center_derivative[DirectionIndex] = 1.0 / rParentElementGeometry.PointsNumber();

    return inner_prod(condition_center - parent_center, rUnitNormalDerivative) +
           inner_prod(condition_center_derivative - parent_center_derivative, rUnitNormal);
}

double RansAdjointUtilities::CalculateWallHeightParentElementDerivative(
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentElementGeometry,
    const IndexType DirectionIndex,
    const array_1d<double, 3>& rUnitNormal,
    const array_1d<double, 3>& rUnitNormalDerivative)
{
    const auto& condition_center = rConditionGeometry.Center();
    const auto& parent_center = rParentElementGeometry.Center();

    array_1d<double, 3> parent_center_derivative = ZeroVector(3);
    parent_center_derivative[DirectionIndex] = 1.0 / rParentElementGeometry.PointsNumber();

    return inner_prod(condition_center - parent_center, rUnitNormalDerivative) -
           inner_prod(parent_center_derivative, rUnitNormal);
}

void RansAdjointUtilities::CalculateYPlusAndUtauDerivative(
    double& rYPlusDerivative,
    double& rUTauDerivative,
    const double YPlus,
    const double WallVelocity,
    const double WallVelocityDerivative,
    const double WallHeight,
    const double WallHeightDerivative,
    const double KinematicViscosity,
    const double Kappa,
    const double Beta,
    const double YPlusLimit)
{
    KRATOS_TRY

    const double u_tau = YPlus * KinematicViscosity / WallHeight;

    if (YPlus > YPlusLimit) {
        // compute logarithmic wall law derivatives
        rUTauDerivative = (WallVelocityDerivative -
                           u_tau * WallHeightDerivative / (Kappa * WallHeight)) /
                          (1 / Kappa + WallVelocity / u_tau);
    } else {
        // compute linear wall law derivatives
        rUTauDerivative =
            (WallVelocityDerivative / WallHeight -
             WallVelocity * WallHeightDerivative / std::pow(WallHeight, 2)) *
            KinematicViscosity / (2 * u_tau);
    }

    rYPlusDerivative =
        (rUTauDerivative * WallHeight + u_tau * WallHeightDerivative) / KinematicViscosity;

    KRATOS_CATCH("");
}

template<>
double RansAdjointUtilities::GeometricalDerivatives<2, 2>::DomainSizeDerivative(
    const GeometryType& rGeometry,
    const IndexType NodeIndex,
    const IndexType DirectionIndex)
{
    const double lx = rGeometry[0].X() - rGeometry[1].X();
    const double lx_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 0);

    const double ly = rGeometry[0].Y() - rGeometry[1].Y();
    const double ly_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 1);

    const double length = lx * lx + ly * ly;
    const double length_derivative = 2.0 * lx * lx_derivative + 2.0 * ly * ly_derivative;

    const double domain_size_derivative = 0.5 * length_derivative / std::sqrt(length);

    return domain_size_derivative;
}

template<>
double RansAdjointUtilities::GeometricalDerivatives<3, 3>::DomainSizeDerivative(
    const GeometryType& rGeometry,
    const IndexType NodeIndex,
    const IndexType DirectionIndex)
{
    const array_1d<double, 3>& a1 = rGeometry[0] - rGeometry[1];
    array_1d<double, 3> a1_derivative = ZeroVector(3);
    a1_derivative[DirectionIndex] = (NodeIndex == 0) - (NodeIndex == 1);

    const double a = norm_2(a1);
    const double a_derivative = CalculateVectorNormDerivative(a, a1, a1_derivative);

    const array_1d<double, 3>& b1 = rGeometry[1] - rGeometry[2];
    array_1d<double, 3> b1_derivative = ZeroVector(3);
    b1_derivative[DirectionIndex] = (NodeIndex == 1) - (NodeIndex == 2);

    const double b = norm_2(b1);
    const double b_derivative = CalculateVectorNormDerivative(b, b1, b1_derivative);

    const array_1d<double, 3>& c1 = rGeometry[2] - rGeometry[0];
    array_1d<double, 3> c1_derivative = ZeroVector(3);
    c1_derivative[DirectionIndex] = (NodeIndex == 2) - (NodeIndex == 0);

    const double c = norm_2(c1);
    const double c_derivative = CalculateVectorNormDerivative(c, c1, c1_derivative);

    const double s = (a+b+c)/2.0;
    const double s_derivative = (a_derivative+b_derivative+c_derivative)/2.0;

    const double value = std::sqrt(s*(s-a)*(s-b)*(s-c));
    double value_derivative = 0.0;

    value_derivative += s_derivative*(s-a)*(s-b)*(s-c);
    value_derivative += s*(s_derivative-a_derivative)*(s-b)*(s-c);
    value_derivative += s*(s-a)*(s_derivative-b_derivative)*(s-c);
    value_derivative += s*(s-a)*(s-b)*(s_derivative-c_derivative);

    value_derivative *= 0.5 / value;

    return value_derivative;
}

void RansAdjointUtilities::CopyAdjointSolutionToNonHistorical(ModelPart& rModelPart)
{
    KRATOS_TRY

    block_for_each(rModelPart.Nodes(), [](ModelPart::NodeType& rNode) {
        rNode.SetValue(ADJOINT_FLUID_VECTOR_1, rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1));
        rNode.SetValue(ADJOINT_FLUID_VECTOR_2, rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2));
        rNode.SetValue(ADJOINT_FLUID_VECTOR_3, rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3));
        rNode.SetValue(AUX_ADJOINT_FLUID_VECTOR_1, rNode.FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1));
        rNode.SetValue(ADJOINT_FLUID_SCALAR_1, rNode.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1));
        rNode.SetValue(RANS_SCALAR_1_ADJOINT_1, rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_1));
        rNode.SetValue(RANS_SCALAR_1_ADJOINT_2, rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_2));
        rNode.SetValue(RANS_SCALAR_1_ADJOINT_3, rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3));
        rNode.SetValue(RANS_AUX_ADJOINT_SCALAR_1, rNode.FastGetSolutionStepValue(RANS_AUX_ADJOINT_SCALAR_1));
        rNode.SetValue(RANS_SCALAR_2_ADJOINT_1, rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_1));
        rNode.SetValue(RANS_SCALAR_2_ADJOINT_2, rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_2));
        rNode.SetValue(RANS_SCALAR_2_ADJOINT_3, rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_3));
        rNode.SetValue(RANS_AUX_ADJOINT_SCALAR_2, rNode.FastGetSolutionStepValue(RANS_AUX_ADJOINT_SCALAR_2));
        rNode.SetValue(SHAPE_SENSITIVITY, rNode.FastGetSolutionStepValue(SHAPE_SENSITIVITY));
    });

    KRATOS_CATCH("");
}

void RansAdjointUtilities::RescaleAdjointSolution(ModelPart& rModelPart)
{
    KRATOS_TRY

    using MultipleReduction =
        CombinedReduction<SumReduction<double>, SumReduction<double>, SumReduction<double>,
                          SumReduction<double>, SumReduction<double>, SumReduction<double>,
                          SumReduction<double>, SumReduction<double>, SumReduction<double>,
                          SumReduction<double>, SumReduction<double>, SumReduction<double>,
                          SumReduction<double>, SumReduction<double>>;

    // now calculate adjoint energy of each adjoint variable
    double adjoint_stable_fluid_vector_1_energy, adjoint_stable_fluid_vector_2_energy,
        adjoint_stable_fluid_vector_3_energy, aux_adjoint_stable_fluid_vector_1_energy,
        adjoint_stable_fluid_scalar_1_energy, rans_scalar_1_adjoint_stable_1_energy,
        rans_scalar_1_adjoint_stable_2_energy, rans_scalar_1_adjoint_stable_3_energy,
        rans_aux_adjoint_stable_scalar_1_energy, rans_scalar_2_adjoint_stable_1_energy,
        rans_scalar_2_adjoint_stable_2_energy, rans_scalar_2_adjoint_stable_3_energy,
        rans_aux_adjoint_stable_scalar_2_energy, shape_sensitivity_stable_energy;

    std::tie(adjoint_stable_fluid_vector_1_energy, adjoint_stable_fluid_vector_2_energy,
             adjoint_stable_fluid_vector_3_energy, aux_adjoint_stable_fluid_vector_1_energy,
             adjoint_stable_fluid_scalar_1_energy, rans_scalar_1_adjoint_stable_1_energy,
             rans_scalar_1_adjoint_stable_2_energy, rans_scalar_1_adjoint_stable_3_energy,
             rans_aux_adjoint_stable_scalar_1_energy, rans_scalar_2_adjoint_stable_1_energy,
             rans_scalar_2_adjoint_stable_2_energy, rans_scalar_2_adjoint_stable_3_energy,
             rans_aux_adjoint_stable_scalar_2_energy, shape_sensitivity_stable_energy) =
        block_for_each<MultipleReduction>(rModelPart.Nodes(), [](const ModelPart::NodeType& rNode) {
            return std::make_tuple(
                std::pow(norm_2(rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1)), 2),
                std::pow(norm_2(rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2)), 2),
                std::pow(norm_2(rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3)), 2),
                std::pow(norm_2(rNode.FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1)), 2),
                std::pow(rNode.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1), 2),
                std::pow(rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_1), 2),
                std::pow(rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_2), 2),
                std::pow(rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3), 2),
                std::pow(rNode.FastGetSolutionStepValue(RANS_AUX_ADJOINT_SCALAR_1), 2),
                std::pow(rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_1), 2),
                std::pow(rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_2), 2),
                std::pow(rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_3), 2),
                std::pow(rNode.FastGetSolutionStepValue(RANS_AUX_ADJOINT_SCALAR_2), 2),
                std::pow(norm_2(rNode.FastGetSolutionStepValue(SHAPE_SENSITIVITY)), 2));
        });

    // here do sum all to support MPI
    // TODO:

    // now calculate adjoint energy of each adjoint variable
    double adjoint_refined_fluid_vector_1_energy, adjoint_refined_fluid_vector_2_energy,
        adjoint_refined_fluid_vector_3_energy, aux_adjoint_refined_fluid_vector_1_energy,
        adjoint_refined_fluid_scalar_1_energy, rans_scalar_1_adjoint_refined_1_energy,
        rans_scalar_1_adjoint_refined_2_energy, rans_scalar_1_adjoint_refined_3_energy,
        rans_aux_adjoint_refined_scalar_1_energy, rans_scalar_2_adjoint_refined_1_energy,
        rans_scalar_2_adjoint_refined_2_energy, rans_scalar_2_adjoint_refined_3_energy,
        rans_aux_adjoint_refined_scalar_2_energy, shape_sensitivity_refined_energy;

    std::tie(adjoint_refined_fluid_vector_1_energy, adjoint_refined_fluid_vector_2_energy,
             adjoint_refined_fluid_vector_3_energy, aux_adjoint_refined_fluid_vector_1_energy,
             adjoint_refined_fluid_scalar_1_energy, rans_scalar_1_adjoint_refined_1_energy,
             rans_scalar_1_adjoint_refined_2_energy, rans_scalar_1_adjoint_refined_3_energy,
             rans_aux_adjoint_refined_scalar_1_energy, rans_scalar_2_adjoint_refined_1_energy,
             rans_scalar_2_adjoint_refined_2_energy, rans_scalar_2_adjoint_refined_3_energy,
             rans_aux_adjoint_refined_scalar_2_energy, shape_sensitivity_refined_energy) =
        block_for_each<MultipleReduction>(rModelPart.Nodes(), [](const ModelPart::NodeType& rNode) {
            return std::make_tuple(
                std::pow(norm_2(rNode.GetValue(ADJOINT_FLUID_VECTOR_1)), 2),
                std::pow(norm_2(rNode.GetValue(ADJOINT_FLUID_VECTOR_2)), 2),
                std::pow(norm_2(rNode.GetValue(ADJOINT_FLUID_VECTOR_3)), 2),
                std::pow(norm_2(rNode.GetValue(AUX_ADJOINT_FLUID_VECTOR_1)), 2),
                std::pow(rNode.GetValue(ADJOINT_FLUID_SCALAR_1), 2),
                std::pow(rNode.GetValue(RANS_SCALAR_1_ADJOINT_1), 2),
                std::pow(rNode.GetValue(RANS_SCALAR_1_ADJOINT_2), 2),
                std::pow(rNode.GetValue(RANS_SCALAR_1_ADJOINT_3), 2),
                std::pow(rNode.GetValue(RANS_AUX_ADJOINT_SCALAR_1), 2),
                std::pow(rNode.GetValue(RANS_SCALAR_2_ADJOINT_1), 2),
                std::pow(rNode.GetValue(RANS_SCALAR_2_ADJOINT_2), 2),
                std::pow(rNode.GetValue(RANS_SCALAR_2_ADJOINT_3), 2),
                std::pow(rNode.GetValue(RANS_AUX_ADJOINT_SCALAR_2), 2),
                std::pow(norm_2(rNode.GetValue(SHAPE_SENSITIVITY)), 2));
        });

    // here do sum all to support MPI
    // TODO:

    // now calculate the ratios
    const double ratio_adjoint_refined_fluid_vector_1_energy = adjoint_refined_fluid_vector_1_energy / adjoint_stable_fluid_vector_1_energy;
    const double ratio_adjoint_refined_fluid_vector_2_energy = adjoint_refined_fluid_vector_2_energy / adjoint_stable_fluid_vector_2_energy;
    const double ratio_adjoint_refined_fluid_vector_3_energy = adjoint_refined_fluid_vector_3_energy / adjoint_stable_fluid_vector_3_energy;
    const double ratio_aux_adjoint_refined_fluid_vector_1_energy = aux_adjoint_refined_fluid_vector_1_energy / aux_adjoint_stable_fluid_vector_1_energy;
    const double ratio_adjoint_refined_fluid_scalar_1_energy = adjoint_refined_fluid_scalar_1_energy / adjoint_stable_fluid_scalar_1_energy;
    const double ratio_rans_scalar_1_adjoint_refined_1_energy = rans_scalar_1_adjoint_refined_1_energy / rans_scalar_1_adjoint_stable_1_energy;
    const double ratio_rans_scalar_1_adjoint_refined_2_energy = rans_scalar_1_adjoint_refined_2_energy / rans_scalar_1_adjoint_stable_2_energy;
    const double ratio_rans_scalar_1_adjoint_refined_3_energy = rans_scalar_1_adjoint_refined_3_energy / rans_scalar_1_adjoint_stable_3_energy;
    const double ratio_rans_aux_adjoint_refined_scalar_1_energy = rans_aux_adjoint_refined_scalar_1_energy / rans_aux_adjoint_stable_scalar_1_energy;
    const double ratio_rans_scalar_2_adjoint_refined_1_energy = rans_scalar_2_adjoint_refined_1_energy / rans_scalar_2_adjoint_stable_1_energy;
    const double ratio_rans_scalar_2_adjoint_refined_2_energy = rans_scalar_2_adjoint_refined_2_energy / rans_scalar_2_adjoint_stable_2_energy;
    const double ratio_rans_scalar_2_adjoint_refined_3_energy = rans_scalar_2_adjoint_refined_3_energy / rans_scalar_2_adjoint_stable_3_energy;
    const double ratio_rans_aux_adjoint_refined_scalar_2_energy = rans_aux_adjoint_refined_scalar_2_energy / rans_aux_adjoint_stable_scalar_2_energy;
    // const double ratio_shape_sensitivity_refined_energy = shape_sensitivity_refined_energy / shape_sensitivity_stable_energy;

    // now rescale the adjoint solution
    block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1) = rNode.GetValue(ADJOINT_FLUID_VECTOR_1) / ratio_adjoint_refined_fluid_vector_1_energy;
        rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2) = rNode.GetValue(ADJOINT_FLUID_VECTOR_2) / ratio_adjoint_refined_fluid_vector_2_energy;
        rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3) = rNode.GetValue(ADJOINT_FLUID_VECTOR_3) / ratio_adjoint_refined_fluid_vector_3_energy;
        rNode.FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1) = rNode.GetValue(AUX_ADJOINT_FLUID_VECTOR_1) / ratio_aux_adjoint_refined_fluid_vector_1_energy;
        rNode.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1) = rNode.GetValue(ADJOINT_FLUID_SCALAR_1) / ratio_adjoint_refined_fluid_scalar_1_energy;

        rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_1) = rNode.GetValue(RANS_SCALAR_1_ADJOINT_1) / ratio_rans_scalar_1_adjoint_refined_1_energy;
        rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_2) = rNode.GetValue(RANS_SCALAR_1_ADJOINT_2) / ratio_rans_scalar_1_adjoint_refined_2_energy;
        rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3) = rNode.GetValue(RANS_SCALAR_1_ADJOINT_3) / ratio_rans_scalar_1_adjoint_refined_3_energy;
        rNode.FastGetSolutionStepValue(RANS_AUX_ADJOINT_SCALAR_1) = rNode.GetValue(RANS_AUX_ADJOINT_SCALAR_1) / ratio_rans_aux_adjoint_refined_scalar_1_energy;

        rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_1) = rNode.GetValue(RANS_SCALAR_2_ADJOINT_1) / ratio_rans_scalar_2_adjoint_refined_1_energy;
        rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_2) = rNode.GetValue(RANS_SCALAR_2_ADJOINT_2) / ratio_rans_scalar_2_adjoint_refined_2_energy;
        rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_3) = rNode.GetValue(RANS_SCALAR_2_ADJOINT_3) / ratio_rans_scalar_2_adjoint_refined_3_energy;
        rNode.FastGetSolutionStepValue(RANS_AUX_ADJOINT_SCALAR_2) = rNode.GetValue(RANS_AUX_ADJOINT_SCALAR_2) / ratio_rans_aux_adjoint_refined_scalar_2_energy;

        rNode.FastGetSolutionStepValue(SHAPE_SENSITIVITY) = rNode.GetValue(SHAPE_SENSITIVITY);
    });

    KRATOS_CATCH("");
}

void RansAdjointUtilities::RescaleShapeSensitivity(ModelPart& rModelPart)
{
    KRATOS_TRY

    const double stable_solution_max_norm = block_for_each<MaxReduction<double>>(rModelPart.Nodes(), [](const ModelPart::NodeType& rNode){
        return norm_2(rNode.FastGetSolutionStepValue(SHAPE_SENSITIVITY));
    });

    const double refined_solution_max_norm = block_for_each<MaxReduction<double>>(rModelPart.Nodes(), [](const ModelPart::NodeType& rNode){
        return norm_2(rNode.GetValue(SHAPE_SENSITIVITY));
    });

    const auto& max_list = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(std::vector<double>{stable_solution_max_norm, refined_solution_max_norm});

    const double rescale_ratio = max_list[0] / max_list[1];

    block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode){
        rNode.FastGetSolutionStepValue(SHAPE_SENSITIVITY) = rNode.GetValue(SHAPE_SENSITIVITY) * rescale_ratio;
    });


    KRATOS_CATCH("");
}

void RansAdjointUtilities::CalculateTransientReponseFunctionInterpolationError(
    ModelPart& rModelPart,
    const double Gamma,
    const double DeltaTime)
{
    KRATOS_TRY

    if (rModelPart.NumberOfElements() == 0) {
        return;
    }

    if (rModelPart.NumberOfNodes() == 0) {
        return;
    }

    const auto& r_front_element = rModelPart.Elements().front();

    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().front().Has(RESPONSE_FUNCTION_INTERPOLATION_ERROR))
        << "RESPONSE_FUNCTION_INTERPOLATION_ERROR is not calculated properly in nodes in "
        << rModelPart.Name() << ".\n";

    KRATOS_ERROR_IF_NOT(r_front_element.Has(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1))
        << "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1 is not calculated properly in elements in "
        << rModelPart.Name() << ".\n";

    KRATOS_ERROR_IF_NOT(r_front_element.Has(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2))
        << "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2 is not calculated properly in elements in "
        << rModelPart.Name() << ".\n";

    KRATOS_ERROR_IF_NOT(r_front_element.Has(RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE))
        << "RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE is not calculated properly in elements in "
        << rModelPart.Name() << ".\n";

    const IndexType number_of_dofs_per_node = r_front_element.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1).size();

    KRATOS_ERROR_IF(r_front_element.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2).size() != number_of_dofs_per_node)
        << "Size mismatch between RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1 and RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2 vectors in " << rModelPart.Name() << "."
        << " [ RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1.size() = " << number_of_dofs_per_node
        << ", RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2.size() = " << r_front_element.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2).size()
        << " ].\n";

    KRATOS_ERROR_IF(r_front_element.GetValue(RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE).size() != number_of_dofs_per_node)
        << "Size mismatch between RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1 and RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE vectors in " << rModelPart.Name() << "."
        << " [ RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1.size() = " << number_of_dofs_per_node
        << ", RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE.size() = " << r_front_element.GetValue(RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE).size()
        << " ].\n";

    const IndexType number_of_elements = rModelPart.GetCommunicator().GlobalNumberOfElements();

    const double coeff_1 = std::abs((Gamma - 1.0) / Gamma);
    const double coeff_2 = std::abs(1.0 / (Gamma * DeltaTime));
    const double coeff_3 = 1.0 / number_of_elements;
    const double coeff_4 = 1.0 / r_front_element.GetGeometry().size();

    const ZeroVector zero(number_of_dofs_per_node);
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode) {
        rNode.SetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY, zero);
    });

    block_for_each(rModelPart.Elements(), Vector(), [&](ElementType& rElement, Vector& rTLS) {
        const double domain_size = rElement.GetGeometry().DomainSize();
        // get a^{n,h}_m
        const auto& r_auxiliary_values_1 = rElement.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1);
        // get b^{n,h}_m
        const auto& r_auxiliary_values_2 = rElement.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2);

        // these currently hold previous time step errors
        auto& r_response_function_dofs_interpolation_error  = rElement.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR);
        auto& r_response_function_dofs_interpolation_error_rate  = rElement.GetValue(RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE);

        KRATOS_DEBUG_ERROR_IF(r_auxiliary_values_1.size() != number_of_dofs_per_node)
            << "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1.size() is "
               "not equal to number of dofs per node in element data container "
               "with id "
            << rElement.Id() << " in " << rModelPart.Name()
            << " [ number of dofs per node = " << number_of_dofs_per_node
            << ", RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1.size() = "
            << r_auxiliary_values_1.size() << " ].\n";

        KRATOS_DEBUG_ERROR_IF(r_auxiliary_values_2.size() != number_of_dofs_per_node)
            << "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2.size() is "
               "not equal to number of dofs per node in element data container "
               "with id "
            << rElement.Id() << " in " << rModelPart.Name()
            << " [ number of dofs per node = " << number_of_dofs_per_node
            << ", RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2.size() = "
            << r_auxiliary_values_2.size() << " ].\n";

        KRATOS_DEBUG_ERROR_IF(r_response_function_dofs_interpolation_error_rate.size() != number_of_dofs_per_node)
            << "RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE.size() is "
               "not equal to number of dofs per node in element data container "
               "with id "
            << rElement.Id() << " in " << rModelPart.Name()
            << " [ number of dofs per node = " << number_of_dofs_per_node
            << ", RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE.size() = "
            << r_response_function_dofs_interpolation_error_rate.size() << " ].\n";

        if (rTLS.size() != number_of_dofs_per_node) {
            rTLS.resize(number_of_dofs_per_node, false);
        }

        for (IndexType i = 0; i < number_of_dofs_per_node; ++i) {
            // now calculate current time step interpolation
            //**************************************************************************************
            rTLS[i] = (coeff_3 +
                       (coeff_1 * r_response_function_dofs_interpolation_error_rate[i] +
                        coeff_2 * r_response_function_dofs_interpolation_error[i]) *
                           r_auxiliary_values_2[i]) *
                      domain_size / std::max(r_auxiliary_values_1[i], 1e-30);

            // update elemental allowed response function interpolation error rate
            r_response_function_dofs_interpolation_error_rate[i] =
                r_response_function_dofs_interpolation_error_rate[i] * coeff_1 +
                coeff_2 * (r_response_function_dofs_interpolation_error[i] + rTLS[i]);

            r_response_function_dofs_interpolation_error[i] = rTLS[i];
            //**************************************************************************************
        }

        // distribute elemental response function interpolation error to nodes
        for (auto& r_node : rElement.GetGeometry()) {
            r_node.SetLock();
            noalias(r_node.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY)) += rTLS * coeff_4;
            r_node.UnSetLock();
        }
    });

    rModelPart.GetCommunicator().AssembleNonHistoricalData(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY);

    block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
        auto& aggregated_interpolation_errors = rNode.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR);
        const auto& current_interpolation_error = rNode.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY);

        for (IndexType i = 0; i < number_of_dofs_per_node; ++i) {
            aggregated_interpolation_errors[i] = 1.0 / (1.0 / aggregated_interpolation_errors[i] + 1.0 / current_interpolation_error[i]);
        }
    });

    KRATOS_CATCH("");
}

} // namespace Kratos