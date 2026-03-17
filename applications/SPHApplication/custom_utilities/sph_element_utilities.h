
#pragma once 

#include <cmath>
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "sph_application_variables.h"

/**
 * @class SPHElementUtilities
 * @brief 
 * @details The methods are static, so it can be called without constructing the class
 */

namespace Kratos
{
class SPHElementUtilities
{
public:

    using SizeType = std::size_t;
    using MatrixType = Matrix;
    using VectorType = Vector;

    /**
     * @brief These following function are copied from structural mechanics application
     * @details For more info go to structural_mechanics/custom_utilities/structural_mechanical_element_utilities.h
     */
    static double GetRayleighAlpha(const Properties& rProperties, const ProcessInfo& rProcessInfo);
    static double GetRayleighBeta(const Properties& rProperties, const ProcessInfo& rProcessInfo);
    static void CalculateRayleighDampingMatrix(Element& rElement, MatrixType& rDampingMatrix, const ProcessInfo& rProcessInfo, const SizeType mat_size);

};

}