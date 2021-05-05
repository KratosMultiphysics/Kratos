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
    const double ratio_shape_sensitivity_refined_energy = shape_sensitivity_refined_energy / shape_sensitivity_stable_energy;

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

} // namespace Kratos