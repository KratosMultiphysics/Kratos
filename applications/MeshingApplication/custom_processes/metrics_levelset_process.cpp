// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// Project includes
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/metrics_math_utils.h"
#include "custom_processes/metrics_levelset_process.h"

namespace Kratos
{
template<SizeType TDim>
ComputeLevelSetSolMetricProcess<TDim>::ComputeLevelSetSolMetricProcess(
        ModelPart& rThisModelPart,
        const Variable<array_1d<double,3>> rVariableGradient,
        Parameters ThisParameters
        ):mThisModelPart(rThisModelPart),
          mVariableGradient(rVariableGradient)
{
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                         : 0.1,
        "enforce_current"                      : true,
        "anisotropy_remeshing"                 : true,
        "anisotropy_parameters":
        {
            "reference_variable_name"              : "DISTANCE",
            "hmin_over_hmax_anisotropic_ratio"      : 1.0,
            "boundary_layer_max_distance"           : 1.0,
            "interpolation"                         : "Linear"
        }
    })" );
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mEnforceCurrent = ThisParameters["enforce_current"].GetBool();

    // In case we have isotropic remeshing (default values)
    if (ThisParameters["anisotropy_remeshing"].GetBool() == false) {
        mRatioReferenceVariable = default_parameters["anisotropy_parameters"]["reference_variable_name"].GetString();
        mAnisotropicRatio = default_parameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
        mBoundLayer = default_parameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
        mInterpolation = ConvertInter(default_parameters["anisotropy_parameters"]["interpolation"].GetString());
    } else {
        mRatioReferenceVariable = ThisParameters["anisotropy_parameters"]["reference_variable_name"].GetString();
        mAnisotropicRatio = ThisParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
        mBoundLayer = ThisParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
        mInterpolation = ConvertInter(ThisParameters["anisotropy_parameters"]["interpolation"].GetString());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void ComputeLevelSetSolMetricProcess<TDim>::Execute()
{
    // Iterate in the nodes
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    const int num_nodes = nodes_array.end() - nodes_array.begin();

    // Some checks
    VariableUtils().CheckVariableExists(mVariableGradient, nodes_array);
    VariableUtils().CheckVariableExists(NODAL_H, nodes_array);

    // Ratio reference variable
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mRatioReferenceVariable)) << "Variable " << mRatioReferenceVariable << " is not a double variable" << std::endl;
    const auto& reference_var = KratosComponents<Variable<double>>::Get(mRatioReferenceVariable);

    // Tensor variable definition
    const Variable<TensorArrayType>& tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_"+std::to_string(TDim)+"D");

    // Setting metric in case not defined
    if (!nodes_array.begin()->Has(tensor_variable)) {
        // Declaring auxiliar vector
        const TensorArrayType aux_zero_vector = ZeroVector(3 * (TDim - 1));
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = nodes_array.begin() + i;
            it_node->SetValue(tensor_variable, aux_zero_vector);
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i)  {
        auto it_node = nodes_array.begin() + i;

        array_1d<double, 3>& gradient_value = it_node->FastGetSolutionStepValue(mVariableGradient);

        // Isotropic by default
        double ratio = 1.0;

        double element_size = mMinSize;
        KRATOS_DEBUG_ERROR_IF_NOT(it_node->SolutionStepsDataHas(NODAL_H)) << "ERROR:: NODAL_H not defined for node " << it_node->Id();
        const double nodal_h = it_node->FastGetSolutionStepValue(NODAL_H);
        if (it_node->SolutionStepsDataHas(reference_var)) {
            const double ratio_reference = it_node->FastGetSolutionStepValue(reference_var);
            ratio = CalculateAnisotropicRatio(ratio_reference, mAnisotropicRatio, mBoundLayer, mInterpolation);
            if (((element_size > nodal_h) && (mEnforceCurrent)) || (std::abs(ratio_reference) > mBoundLayer))
                element_size = nodal_h;
        } else {
            if (((element_size > nodal_h) && (mEnforceCurrent)))
                element_size = nodal_h;
        }

        // For postprocess pourposes
        it_node->SetValue(ANISOTROPIC_RATIO, ratio);

        const double tolerance = 1.0e-12;
        const double norm_gradient_value = norm_2(gradient_value);
        if (norm_gradient_value > tolerance)
            gradient_value /= norm_gradient_value;

        // We compute the metric
        KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(tensor_variable)) << "METRIC_TENSOR_" + std::to_string(TDim) + "D  not defined for node " << it_node->Id() << std::endl;
        TensorArrayType& metric = it_node->GetValue(tensor_variable);

        const double norm_metric = norm_2(metric);
        if (norm_metric > 0.0) { // NOTE: This means we combine differents metrics, at the same time means that the metric should be reseted each time
            const TensorArrayType& old_metric = it_node->GetValue(tensor_variable);
            const TensorArrayType& new_metric = ComputeLevelSetMetricTensor(gradient_value, ratio, element_size);

            metric = MetricsMathUtils<TDim>::IntersectMetrics(old_metric, new_metric);
        } else {
            metric = ComputeLevelSetMetricTensor(gradient_value, ratio, element_size);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double, 3> ComputeLevelSetSolMetricProcess<2>::ComputeLevelSetMetricTensor(
    const array_1d<double, 3>& GradientValue,
    const double Ratio,
    const double ElementSize
    )
{
    array_1d<double, 3> metric;

    const double coeff_0 = 1.0/(ElementSize * ElementSize);
    const double coeff_1 = coeff_0/(Ratio * Ratio);

    const double v0v0 = GradientValue[0]*GradientValue[0];
    const double v0v1 = GradientValue[0]*GradientValue[1];
    const double v1v1 = GradientValue[1]*GradientValue[1];

    metric[0] = coeff_0*(1.0 - v0v0) + coeff_1*v0v0;
    metric[1] = coeff_0*(1.0 - v1v1) + coeff_1*v1v1;
    metric[2] = coeff_0*(    - v0v1) + coeff_1*v0v1;

    return metric;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double, 6> ComputeLevelSetSolMetricProcess<3>::ComputeLevelSetMetricTensor(
    const array_1d<double, 3>& GradientValue,
    const double Ratio,
    const double ElementSize
    )
{
    array_1d<double, 6> metric;

    const double coeff_0 = 1.0/(ElementSize * ElementSize);
    const double coeff_1 = coeff_0/(Ratio * Ratio);

    const double v0v0 = GradientValue[0]*GradientValue[0];
    const double v0v1 = GradientValue[0]*GradientValue[1];
    const double v0v2 = GradientValue[0]*GradientValue[2];
    const double v1v1 = GradientValue[1]*GradientValue[1];
    const double v1v2 = GradientValue[1]*GradientValue[2];
    const double v2v2 = GradientValue[2]*GradientValue[2];

    metric[0] = coeff_0*(1.0 - v0v0) + coeff_1*v0v0;
    metric[1] = coeff_0*(1.0 - v1v1) + coeff_1*v1v1;
    metric[2] = coeff_0*(1.0 - v2v2) + coeff_1*v2v2;
    metric[3] = coeff_0*(    - v0v1) + coeff_1*v0v1;
    metric[4] = coeff_0*(    - v1v2) + coeff_1*v1v2;
    metric[5] = coeff_0*(    - v0v2) + coeff_1*v0v2;

    return metric;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
double ComputeLevelSetSolMetricProcess<TDim>::CalculateAnisotropicRatio(
    const double Distance,
    const double AnisotropicRatio,
    const double BoundLayer,
    const Interpolation& rInterpolation
    )
{
    const double tolerance = 1.0e-12;
    double ratio = 1.0; // NOTE: Isotropic mesh
    if (AnisotropicRatio < 1.0) {
        if (std::abs(Distance) <= BoundLayer) {
            if (rInterpolation == Interpolation::CONSTANT)
                ratio = AnisotropicRatio;
            else if (rInterpolation == Interpolation::LINEAR)
                ratio = AnisotropicRatio + (std::abs(Distance)/BoundLayer) * (1.0 - AnisotropicRatio);
            else if (rInterpolation == Interpolation::EXPONENTIAL) {
                ratio = - std::log(std::abs(Distance)/BoundLayer) * AnisotropicRatio + tolerance;
                if (ratio > 1.0) ratio = 1.0;
            }
        }
    }

    return ratio;
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeLevelSetSolMetricProcess<2>;
template class ComputeLevelSetSolMetricProcess<3>;

};// namespace Kratos.
