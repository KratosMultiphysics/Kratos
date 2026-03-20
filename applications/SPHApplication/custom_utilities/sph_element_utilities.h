
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

    static void GetLocalBodyForces(Element& rElement, VectorType& body_force);


};

}