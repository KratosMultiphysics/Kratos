
#pragma once 

#include <cmath>
#include "includes/model_part.h"



/**
 * @class NeighboursSearchUtilities
 * @brief This class includes some utilities necessaries for the computation of particle neighbours
 * @details The methods are static, so it can be called without constructing the class
 */

namespace Kratos
{

class NeighboursSearchUtilities
{
public:

static double ComputeSmoothingLength(const ModelPart& rThisModelPart, double Coeff);

static double ComputeInterparticleMinDist(const ModelPart& rThisModelPart);

};

}
