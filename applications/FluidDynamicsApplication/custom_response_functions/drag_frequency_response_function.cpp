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
        "component_type"            : "",
        "frequency_bin_index"       : -1,
        "window_time_length"        : -1.0,
        "echo_level"                : 0
    })");

    Settings.ValidateAndAssignDefaults(default_settings);

    mStructureModelPartName = Settings["structure_model_part_name"].GetString();

    mStartTime = 0.0;

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

    const auto& component_type = Settings["component_type"].GetString();
    if (component_type == "real") {
        mIsRealComponentRequested = true;
    } else if (component_type == "imag") {
        mIsRealComponentRequested = false;
    } else {
        KRATOS_ERROR << "Unsupported component type requested [ component_type = \"" << component_type << "\" ]. Supported component types are:\n   \"real\"\n   \"imag\"\n";
    }

    mWindowingLength = Settings["window_time_length"].GetDouble();
    KRATOS_ERROR_IF(mWindowingLength <= 0)
        << "Window time length should be greater than to zero. [ "
        << "window_time_length = " << mWindowingLength << " ].\n";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void DragFrequencyResponseFunction<TDim>::InitializeSolutionStep()
{
    KRATOS_TRY

    const auto& r_process_info = mrModelPart.GetProcessInfo();
    const double current_time = r_process_info[TIME];

    if (!mIsInitialized) {
        mIsInitialized = true;

        // since transient adjoints are run in backwards in time
        // in here it is assumed that always transient simulations start with time = 0.0s
        const double delta_time = -1.0 * r_process_info[DELTA_TIME];

        mpFluidFFTUtilities = new FluidFFTUtilities(current_time, mWindowingLength, delta_time);

        KRATOS_INFO_IF("DragFrequencyResponseFunction", mEchoLevel > 0) << "Summary: \n"
            << "\tMain model part name      : " << mrModelPart.Name() << "\n"
            << "\tStructural model part name: " << mStructureModelPartName << "\n"
            << "\tDrag direction            : " << mDragDirection << "\n"
            << "\tTime step length          : " << delta_time << "s\n"
            << "\tWindowing duration        : " << mWindowingLength << "s\n"
            << "\tTotal duration            : " << current_time << "s\n"
            << "\tFrequency resolution      : " << mpFluidFFTUtilities->GetFrequencyResolution() << "Hz\n"
            << "\tEvaluated frequency bin   : " << mFrequencyBinIndex << "\n"
            << "\tEvaluated frequency       : " << mpFluidFFTUtilities->GetFrequency(mFrequencyBinIndex) << "Hz \n"
            << "\tMaximum possible frequency: " << mpFluidFFTUtilities->GetMaximumFrequency() << "Hz \n"
            << "\tComponent type            : " << (mIsRealComponentRequested ? "real\n" : "imaginary\n");
    }

    BaseType::InitializeSolutionStep();

    // Fix coefficients for all elements since they are only time dependent
    mIsWithinWindowingRange = mpFluidFFTUtilities->IsWithinWindowingRange(current_time);

    if (mIsRealComponentRequested) {
        mCurrentTimeStepCoefficient = mpFluidFFTUtilities->CalculateHannWindowCoefficient(current_time) * mpFluidFFTUtilities->CalculateFFTRealCoefficient(mFrequencyBinIndex, current_time);
    } else {
        mCurrentTimeStepCoefficient = mpFluidFFTUtilities->CalculateHannWindowCoefficient(current_time) * mpFluidFFTUtilities->CalculateFFTImagCoefficient(mFrequencyBinIndex, current_time);
    }

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
        rResidualGradient, rAdjointElement.GetGeometry().Points(), rProcessInfo, rResponseGradient);
}


template <unsigned int TDim>
void DragFrequencyResponseFunction<TDim>::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    CalculateDragFrequencyContribution(
        rResidualGradient, rAdjointElement.GetGeometry().Points(), rProcessInfo, rResponseGradient);
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
                                        rProcessInfo, rResponseGradient);
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
                                        rProcessInfo, rSensitivityGradient);

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
    const ProcessInfo& rProcessInfo,
    Vector& rDerivativesOfDrag) const
{
    KRATOS_TRY

    if (mIsWithinWindowingRange) {
        // calculate raw drag sensitivities
        BaseType::CalculateDragContribution(rDerivativesOfResidual, rNodes, rDerivativesOfDrag, rProcessInfo);

        // apply hann windowing and fft component coefficient
        rDerivativesOfDrag *= mCurrentTimeStepCoefficient;

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