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
        const Variable<array_1d<double,3>>& rVariableGradient,
        Parameters ThisParameters
        ):mThisModelPart(rThisModelPart),
          mVariableGradient(rVariableGradient)
{
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                         : 0.1,
        "maximal_size"                         : 1.0,
        "sizing_parameters":
        {
            "reference_variable_name"          : "DISTANCE",
            "boundary_layer_max_distance"      : 1.0,
            "interpolation"                    : "constant"
        },
        "enforce_current"                      : true,
        "anisotropy_remeshing"                 : true,
        "anisotropy_parameters":
        {
            "reference_variable_name"              : "DISTANCE",
            "hmin_over_hmax_anisotropic_ratio"      : 1.0,
            "boundary_layer_max_distance"           : 1.0,
            "interpolation"                         : "linear"
        }
    })" );
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();
    mSizeReferenceVariable = ThisParameters["sizing_parameters"]["reference_variable_name"].GetString();
    mSizeBoundLayer = ThisParameters["sizing_parameters"]["boundary_layer_max_distance"].GetDouble();
    mSizeInterpolation = ConvertInter(ThisParameters["sizing_parameters"]["interpolation"].GetString());
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
    NodesArrayType& r_nodes_array = mThisModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    const int num_nodes = r_nodes_array.end() - it_node_begin;

    // Some checks
    VariableUtils().CheckVariableExists(mVariableGradient, r_nodes_array);
    for (auto& i_node : r_nodes_array)
        KRATOS_ERROR_IF_NOT(i_node.Has(NODAL_H)) << "NODAL_H must be computed" << std::endl;

    // Ratio reference variable
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mRatioReferenceVariable)) << "Variable " << mRatioReferenceVariable << " is not a double variable" << std::endl;
    const auto& ani_reference_var = KratosComponents<Variable<double>>::Get(mRatioReferenceVariable);

    // Size reference variable
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mSizeReferenceVariable)) << "Variable " << mSizeReferenceVariable << " is not a double variable" << std::endl;
    const auto& r_size_reference_var = KratosComponents<Variable<double>>::Get(mSizeReferenceVariable);

    // Tensor variable definition
    const Variable<TensorArrayType>& r_tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_" +   std::to_string(TDim) + "D");

    // Setting r_metric in case not defined
    if (!it_node_begin->Has(r_tensor_variable)) {
        // Declaring auxiliar vector
        const TensorArrayType aux_zero_vector = ZeroVector(3 * (TDim - 1));
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;
            it_node->SetValue(r_tensor_variable, aux_zero_vector);
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i)  {
        auto it_node = it_node_begin + i;

        array_1d<double, 3>& r_gradient_value = it_node->FastGetSolutionStepValue(mVariableGradient);

        // Isotropic by default
        double ratio = 1.0;
        if (it_node->SolutionStepsDataHas(ani_reference_var)) {
            const double ratio_reference = it_node->FastGetSolutionStepValue(ani_reference_var);
            ratio = CalculateAnisotropicRatio(ratio_reference);
        }

        // MinSize by default
        double element_size = mMinSize;
        const double nodal_h = it_node->GetValue(NODAL_H);
        if (it_node->SolutionStepsDataHas(r_size_reference_var)) {
            const double size_reference = it_node->FastGetSolutionStepValue(r_size_reference_var);
            element_size = CalculateElementSize(size_reference, nodal_h);
            if (((element_size > nodal_h) && (mEnforceCurrent)) || (std::abs(size_reference) > mSizeBoundLayer))
                element_size = nodal_h;
        } else {
            if (((element_size > nodal_h) && (mEnforceCurrent)))
                element_size = nodal_h;
        }

        // For postprocess pourposes
        it_node->SetValue(ANISOTROPIC_RATIO, ratio);

        const double tolerance = 1.0e-12;
        const double norm_gradient_value = norm_2(r_gradient_value);
        if (norm_gradient_value > tolerance)
            r_gradient_value /= norm_gradient_value;

        // We compute the metric
        KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(r_tensor_variable)) << "METRIC_TENSOR_" + std::to_string(TDim) + "D  not defined for node " << it_node->Id() << std::endl;
        TensorArrayType& r_metric = it_node->GetValue(r_tensor_variable);

        const double norm_metric = norm_2(r_metric);
        if (norm_metric > 0.0) { // NOTE: This means we combine differents metrics, at the same time means that the r_metric should be reseted each time
            const TensorArrayType& r_old_metric = it_node->GetValue(r_tensor_variable);
            const TensorArrayType new_metric = ComputeLevelSetMetricTensor(r_gradient_value, ratio, element_size);

            noalias(r_metric) = MetricsMathUtils<TDim>::IntersectMetrics(r_old_metric, new_metric);
        } else {
            noalias(r_metric) = ComputeLevelSetMetricTensor(r_gradient_value, ratio, element_size);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double, 3> ComputeLevelSetSolMetricProcess<2>::ComputeLevelSetMetricTensor(
    const array_1d<double, 3>& rGradientValue,
    const double Ratio,
    const double ElementSize
    )
{
    array_1d<double, 3> r_metric;

    const double coeff_0 = 1.0/(ElementSize * ElementSize);
    const double coeff_1 = coeff_0/(Ratio * Ratio);

    const double v0v0 = rGradientValue[0]*rGradientValue[0];
    const double v0v1 = rGradientValue[0]*rGradientValue[1];
    const double v1v1 = rGradientValue[1]*rGradientValue[1];

    r_metric[0] = coeff_0*(1.0 - v0v0) + coeff_1*v0v0;
    r_metric[1] = coeff_0*(1.0 - v1v1) + coeff_1*v1v1;
    r_metric[2] = coeff_0*(    - v0v1) + coeff_1*v0v1;

    return r_metric;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double, 6> ComputeLevelSetSolMetricProcess<3>::ComputeLevelSetMetricTensor(
    const array_1d<double, 3>& rGradientValue,
    const double Ratio,
    const double ElementSize
    )
{
    array_1d<double, 6> r_metric;

    const double coeff_0 = 1.0/(ElementSize * ElementSize);
    const double coeff_1 = coeff_0/(Ratio * Ratio);

    const double v0v0 = rGradientValue[0]*rGradientValue[0];
    const double v0v1 = rGradientValue[0]*rGradientValue[1];
    const double v0v2 = rGradientValue[0]*rGradientValue[2];
    const double v1v1 = rGradientValue[1]*rGradientValue[1];
    const double v1v2 = rGradientValue[1]*rGradientValue[2];
    const double v2v2 = rGradientValue[2]*rGradientValue[2];

    r_metric[0] = coeff_0*(1.0 - v0v0) + coeff_1*v0v0;
    r_metric[1] = coeff_0*(1.0 - v1v1) + coeff_1*v1v1;
    r_metric[2] = coeff_0*(1.0 - v2v2) + coeff_1*v2v2;
    r_metric[3] = coeff_0*(    - v0v1) + coeff_1*v0v1;
    r_metric[4] = coeff_0*(    - v1v2) + coeff_1*v1v2;
    r_metric[5] = coeff_0*(    - v0v2) + coeff_1*v0v2;

    return r_metric;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
double ComputeLevelSetSolMetricProcess<TDim>::CalculateAnisotropicRatio(
    const double Distance
    )
{
    const double tolerance = 1.0e-12;
    double ratio = 1.0; // NOTE: Isotropic mesh
    if (mAnisotropicRatio < 1.0) {
        if (std::abs(Distance) <= mBoundLayer) {
            if (mInterpolation == Interpolation::CONSTANT)
                ratio = mAnisotropicRatio;
            else if (mInterpolation == Interpolation::LINEAR)
                ratio = mAnisotropicRatio + (std::abs(Distance)/mBoundLayer) * (1.0 - mAnisotropicRatio);
            else if (mInterpolation == Interpolation::EXPONENTIAL) {
                ratio = - std::log(std::abs(Distance)/mBoundLayer) * mAnisotropicRatio + tolerance;
                if (ratio > 1.0) ratio = 1.0;
            }
        }
    }

    return ratio;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
double ComputeLevelSetSolMetricProcess<TDim>::CalculateElementSize(
    const double Distance,
    const double NodalH
    )
{
    double size = NodalH;
    if (std::abs(Distance) <= mSizeBoundLayer) {
        if (mSizeInterpolation == Interpolation::CONSTANT)
            size = mMinSize;
        else if (mSizeInterpolation == Interpolation::LINEAR)
            size = mMinSize + (std::abs(Distance)/mSizeBoundLayer) * (mMaxSize-mMinSize);
        else if (mSizeInterpolation == Interpolation::EXPONENTIAL) {
            size = - std::log(1-std::abs(Distance)/mSizeBoundLayer) * (mMaxSize-mMinSize) + mMinSize;
            if (size > mMaxSize) size = mMaxSize;
        }
    }


    return size;
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeLevelSetSolMetricProcess<2>;
template class ComputeLevelSetSolMetricProcess<3>;

};// namespace Kratos.
