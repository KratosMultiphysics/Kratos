//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "smooth_clamper.h"

namespace Kratos {

SmoothClamper::SmoothClamper(
    const double Min,
    const double Max)
    : mMin(Min),
      mMax(Max),
      mDelta(mMax - mMin)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mMin >= mMax) << "Min must be lower than Max";

    KRATOS_CATCH("");
}

double SmoothClamper::ProjectForward(const double X) const
{
    const double x_tilde = std::clamp((X - mMin) / mDelta, 0.0, 1.0);
    return mMin + x_tilde * x_tilde * (3.0 - 2.0 * x_tilde) * mDelta;
}

double SmoothClamper::CalculateForwardProjectionGradient(const double X) const
{
    const double x_tilde = std::clamp((X - mMin) / mDelta, 0.0, 1.0);
    return 6 * x_tilde - 6 * x_tilde * x_tilde;
}

double SmoothClamper::ProjectBackward(const double Y) const
{
    double x_tilde;
    if (Y < mMin) {
        x_tilde = 0;
    } else if (Y > mMax) {
        x_tilde = 1.0;
    } else {
        const double y = (Y - mMin) / mDelta;
        x_tilde = 0.5 - std::sin(std::asin(1.0 - 2.0 * y) / 3.0);
    }
    return mMin + x_tilde * mDelta;
}

TensorAdaptor<double>::Pointer SmoothClamper::ProjectForward(const TensorAdaptor<double>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);

    auto p_output = rInput.Clone();
    auto input_data_view = rInput.ViewData();
    auto output_data_view = p_output->ViewData();

    IndexPartition<IndexType>(rInput.Size()).for_each([this, &output_data_view, &input_data_view](const auto Index) {
        output_data_view[Index]= this->ProjectForward(input_data_view[Index]);
    });

    return p_output;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SmoothClamper::CalculateForwardProjectionGradient(const TensorAdaptor<double>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);
    // y = 3x^2 - 2x^3
    // dy/dx = 3.2.x - 2.3.x^2

    auto p_output = rInput.Clone();
    auto input_data_view = rInput.ViewData();
    auto output_data_view = p_output->ViewData();

    IndexPartition<IndexType>(rInput.Size()).for_each([this, &output_data_view, &input_data_view](const auto Index) {
        output_data_view[Index]= this->CalculateForwardProjectionGradient(input_data_view[Index]);
    });

    return p_output;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SmoothClamper::ProjectBackward(const TensorAdaptor<double>& rInput) const
{
    KRATOS_TRY

    auto p_output = rInput.Clone();
    auto input_data_view = rInput.ViewData();
    auto output_data_view = p_output->ViewData();

    IndexPartition<IndexType>(rInput.Size()).for_each([this, &output_data_view, &input_data_view](const auto Index) {
        output_data_view[Index]= this->ProjectBackward(input_data_view[Index]);
    });


    return p_output;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/