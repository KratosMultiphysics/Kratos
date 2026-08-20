//  ____  ____  _   _                   _ _           _   _             
// / ___||  _ \| | | | __ _ _ __  _ __ | (_) ___ __ _| |_(_) ___  _ __  
// \___ \| |_) | |_| |/ _` | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_ \ 
//  ___) |  __/|  _  | (_| | |_) | |_) | | | (_| (_| | |_| | (_) | | | |
// |____/|_|   |_| |_|\__,_| .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                         |_|   |_|                                    

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once 

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
    /**
     * @details The computation of the smoothing length and the subsequent search of the neighbours are based on the 
     * global minimum distance between the particles. This means that this class is not suitable for simulations with non-uniform particle distributions. 
     * In that case, the possible solutions are:
     *  - Different smoothing lengths for each particle.
     *  - Selecting a priori the number of closet neighbours to be considered for each particle.   
     */

/**
 * @brief This method computes the smoothing length of the particles in the model part
 * @param Coeff The coefficient to multiply the smoothing length
*/
static double ComputeSmoothingLength(const ModelPart& rThisModelPart, double Coeff);

/**
 * @brief This method computes the minimum distance between particles in the model part
*/
static double ComputeInterparticleMinDist(const ModelPart& rThisModelPart);

/**
 * @brief This method computes the neighbours of the particles in the model part
 * @param Radius The radius to search for neighbours
*/
static void SearchNeighbours(ModelPart& rThisModelPart, double Radius);

};

}
