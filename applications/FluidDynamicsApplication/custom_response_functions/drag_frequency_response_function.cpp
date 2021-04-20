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
#include <functional>
#include <cmath>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "drag_frequency_response_function.h"

namespace Kratos
{

template <unsigned int TDim>
DragFrequencyResponseFunction<TDim>::DragFrequencyResponseFunction(
    Parameters Settings,
    ModelPart& rModelPart)
    : BaseType(rModelPart)
{
    KRATOS_TRY;

    Parameters default_settings(R"(
    {
        "structure_model_part_name" : "PLEASE_SPECIFY_STRUCTURE_MODEL_PART",
        "drag_direction"            : [1.0, 0.0, 0.0],
        "is_real_component"         : true,
        "frequency_bin_index"       : 1,
        "window_size"               : 300,
        "total_number_of_time_steps": 1000,
        "echo_level"                : 0
        "
    })");

    Settings.ValidateAndAssignDefaults(default_settings);

    mStructureModelPartName = Settings["structure_model_part_name"].GetString();

    if (Settings["drag_direction"].IsArray() == false ||
        Settings["drag_direction"].size() != 3) {
        KRATOS_ERROR << "Invalid \"drag_direction\"." << std::endl;
    }

    for (unsigned int d = 0; d < TDim; ++d)
        mDragDirection[d] = Settings["drag_direction"][d].GetDouble();

    if (std::abs(norm_2(mDragDirection) - 1.0) > 1e-3) {
        const double magnitude = norm_2(mDragDirection);
        if (magnitude == 0.0)
            KRATOS_ERROR << "\"drag_direction\" is zero." << std::endl;

        KRATOS_WARNING("DragFrequencyResponseFunction")
            << "Non unit magnitude in \"drag_direction\"." << std::endl;
        KRATOS_WARNING("DragFrequencyResponseFunction") << "Normalizing ..." << std::endl;

        for (unsigned int d = 0; d < TDim; ++d)
            mDragDirection[d] /= magnitude;
    }

    mEchoLevel = Settings["echo_level"].GetInt();

    mFrequencyBinIndex = Settings["frequency_bin_index"].GetInt();
    KRATOS_ERROR_IF(mFrequencyBinIndex < 0)
        << "Frequency bin index should be greater than or equal to zero. [ "
           "frequency_bin_index = "
        << mFrequencyBinIndex << " ].\n";

    mTotalNumberOfTimeSteps = Settings["total_number_of_time_steps"].GetInt();
    KRATOS_ERROR_IF(mTotalNumberOfTimeSteps <= 0)
        << "Total number of time steps should be greater than to zero. [ "
           "total_number_of_time_steps = "
        << mTotalNumberOfTimeSteps << " ].\n";

    mIsRealComponentRequested = Settings["is_real_component"].GetBool();
    if (mIsRealComponentRequested) {
        mComponentFunction = [](double x) { return std::cos(x); };
    } else {
        mComponentFunction = [](double x) { return std::sin(x); };
    }

    mWindowSize = Settings["window_size"].GetInt();
    KRATOS_ERROR_IF(mWindowSize <= 0)
        << "Window size should be greater than to zero. [ "
        << "window_size = " << mWindowSize << " ].\n";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void DragFrequencyResponseFunction<TDim>::Initialize()
{
    KRATOS_TRY;

    BaseType::Check();

    VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Nodes());
    VariableUtils().SetFlag(STRUCTURE, true,  mrModelPart.GetSubModelPart(mStructureModelPartName).Nodes());

    const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];

    KRATOS_INFO_IF("DragFrequencyResponseFunction", mEchoLevel > 0) << "Summary: \n"
        << "\tMain model part name      : " << mrModelPart.Name() << "\n"
        << "\tStructural model part name: " << mStructureModelPartName << "\n"
        << "\tDrag direction            : " << mDragDirection << "\n"
        << "\tWindow size               : " << mWindowSize << "\n"
        << "\tWindow duration           : " << mWindowSize * delta_time << " s\n"
        << "\tEvaluated frequency bin   : " << mFrequencyBinIndex << "\n"
        << "\tEvaluated frequency       : " << mFrequencyBinIndex / (delta_time * mTotalNumberOfTimeSteps) << " Hz \n"
        << "\tMaximum possible frequency: " << 1.0 / (delta_time * 2.0) << " Hz \n"
        << "\tTotal number of time steps: " << mTotalNumberOfTimeSteps << "\n"
        << "\tComponent type            : " << (mIsRealComponentRequested ? "real\n" : "imaginary\n");

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void DragFrequencyResponseFunction<TDim>::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDragFrequencyContribution(
        rResidualGradient, rAdjointElement.GetGeometry().Points(), rProcessInfo[STEP], rResponseGradient);
}


template <unsigned int TDim>
void DragFrequencyResponseFunction<TDim>::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDragFrequencyContribution(
        rResidualGradient, rAdjointElement.GetGeometry().Points(), rProcessInfo[STEP], rResponseGradient);
}


template <unsigned int TDim>
void DragFrequencyResponseFunction<TDim>::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDragFrequencyContribution(rResidualGradient,
                                        rAdjointElement.GetGeometry().Points(),
                                        rProcessInfo[STEP], rResponseGradient);
}

template <unsigned int TDim>
void DragFrequencyResponseFunction<TDim>::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    CalculateDragFrequencyContribution(rSensitivityMatrix,
                                        rAdjointElement.GetGeometry().Points(),
                                        rProcessInfo[STEP], rSensitivityGradient);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
double DragFrequencyResponseFunction<TDim>::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY;

    KRATOS_ERROR
        << "DragFrequencyResponseFunction::CalculateValue(ModelPart& "
            "rModelPart) is not implemented!!!\n";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void DragFrequencyResponseFunction<TDim>::CalculateDragFrequencyContribution(
    const Matrix& rDerivativesOfResidual,
    const Element::NodesArrayType& rNodes,
    const int Step,
    Vector& rDerivativesOfDrag) const
{
    KRATOS_TRY

    const int offsetted_step = Step - mTotalNumberOfTimeSteps + mWindowSize;
    if (offsetted_step >= 0) {
        // calculate raw drag sensitivities
        BaseType::CalculateDragContribution(rDerivativesOfResidual, rNodes, rDerivativesOfDrag);

        const double component_coefficient = mComponentFunction(2 * M_PI * Step * mFrequencyBinIndex / mTotalNumberOfTimeSteps);
        const double windowing_value = 0.5 * (1.0 - std::cos(2.0 * M_PI * offsetted_step / mWindowSize));

        rDerivativesOfDrag *= (component_coefficient * windowing_value);
    } else {
        if (rDerivativesOfDrag.size() != rDerivativesOfResidual.size1()) {
            rDerivativesOfDrag.resize(rDerivativesOfResidual.size1(), false);
        }

        rDerivativesOfDrag.clear();
    }

    KRATOS_CATCH("");
}

// template instantiations
template class DragFrequencyResponseFunction<2>;
template class DragFrequencyResponseFunction<3>;

} /* namespace Kratos.*/