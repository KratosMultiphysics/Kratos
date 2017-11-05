#include "numerical_diffusion.h"

namespace Kratos
{
    NumericalDiffusion::NUMERICAL_DIFFUSION_METHOD NumericalDiffusion::mNumericalDiffusionMethod = NumericalDiffusion::NUMERICAL_DIFFUSION_METHOD::SINGULAR_VALUE_PRESSURE_COUPLED;
    double NumericalDiffusion::mBeta = 0.0;
}
