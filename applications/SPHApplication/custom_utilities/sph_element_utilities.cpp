#pragma once

#include "custom_utilities/sph_element_utilities.h"

namespace Kratos
{

double SPHElementUtilities::GetRayleighAlpha(const Properties& rProperties, const ProcessInfo& rProcessInfo)
{
    if (rProperties.Has(RAYLEIGH_ALPHA)) {
        return rProperties[RAYLEIGH_ALPHA];
    } else if (rProcessInfo.Has(RAYLEIGH_ALPHA)) {
        return rProcessInfo[RAYLEIGH_ALPHA];
    }

    return 0.0;
}

double SPHElementUtilities::GetRayleighBeta(const Properties& rProperties, const ProcessInfo& rProcessInfo)
{
    if (rProperties.Has(RAYLEIGH_BETA)) {
        return rProperties[RAYLEIGH_BETA];
    } else if (rProcessInfo.Has(RAYLEIGH_BETA)) {
        return rProcessInfo[RAYLEIGH_BETA];
    }

    return 0.0;
}

void SPHElementUtilities::CalculateRayleighDampingMatrix(
    Element& rElement,
    MatrixType& rDampingMatrix,
    const ProcessInfo& rProcessInfo,
    const SizeType mat_size
)
{
    KRATOS_TRY
    // Rayleigh Damping Matrix: alpha*M + beta*K

    const double alpha = GetRayleighAlpha(rElement.GetProperties(), rProcessInfo);
    const double beta = GetRayleighBeta(rElement.GetProperties(), rProcessInfo);

    if (std::abs(alpha) < 1E-12 && std::abs(beta) < 1E-12) {
        if (rDampingMatrix.size1() != mat_size || rDampingMatrix.size2() != mat_size) {
            rDampingMatrix.resize(mat_size, mat_size, false);
        }
        noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);
    } else if (std::abs(alpha) > 1E-12 && std::abs(beta) < 1E-12) {
        rElement.CalculateMassMatrix(rDampingMatrix, rProcessInfo);
        rDampingMatrix *= alpha;
    } else if (std::abs(alpha) < 1E-12 && std::abs(beta) > 1E-12) {
        rElement.CalculateLeftHandSide(rDampingMatrix, rProcessInfo); 
        rDampingMatrix *= beta;
    } else {
        rElement.CalculateLeftHandSide(rDampingMatrix, rProcessInfo); 
        rDampingMatrix *= beta;

        Matrix MassMatrix;
        rElement.CalculateMassMatrix(MassMatrix, rProcessInfo);
        noalias(rDampingMatrix) += alpha  * MassMatrix;
    }
    KRATOS_CATCH("")
}


}