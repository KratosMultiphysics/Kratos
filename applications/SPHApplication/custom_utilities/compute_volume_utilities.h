
#pragma once 

#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "sph_application_variables.h"

/**
 * @class ComputeVolumeUtilities
 * @brief 
 * @details The methods are static, so it can be called without constructing the class
 */

namespace Kratos
{
class ComputeVolumeUtilities
{
public:

    static void CalculateBoundaryNetNormal(Geometry<Node>& rGeom);

    static void CheckBoundaryNetNormal(ModelPart::NodesContainerType& rNodes);


};

}