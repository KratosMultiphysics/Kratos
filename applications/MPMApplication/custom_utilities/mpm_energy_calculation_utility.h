//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/element.h"

namespace Kratos
{

class MPMEnergyCalculationUtility
{
public:

    /**
     * @brief Compute the potential energy of a material point element
     * @param rElement The material point element whose potential energy is to be computed
     * @return The potential energy of the input material point element
     */
    static double CalculatePotentialEnergy(Element& rElement);

    /**
     * @brief Compute the potential energy of a model part
     * @param rModelPart The model part whose potential energy is to be computed
     * @return The potential energy of the input model part
     */
    static double CalculatePotentialEnergy(ModelPart& rModelPart);

    /**
     * @brief Compute the kinetic energy of a material point element
     * @param rElement The material point element whose kinetic energy is to be computed
     * @return The kinetic energy of the input material point element
     */
    static double CalculateKineticEnergy(Element& rElement);

    /**
     * @brief Compute the kinetic energy of a model part
     * @param rModelPart The model part whose kinetic energy is to be computed
     * @return The kinetic energy of the input model part
     */
    static double CalculateKineticEnergy(ModelPart& rModelPart);

    /**
     * @brief Compute the strain energy of a material point element
     * @param rElement The material point element whose strain energy is to be computed
     * @return The strain energy of the input material point element
     */
    static double CalculateStrainEnergy(Element& rElement);

    /**
     * @brief Compute the strain energy of a model part
     * @param rModelPart The model part whose strain energy is to be computed
     * @return The strain energy of the input model part
     */
    static double CalculateStrainEnergy(ModelPart& rModelPart);

    /**
     * @brief Compute the total energy of a material point element
     * @detail The total energy is computed summing the kinetic, strain and potential energy
     * @param rElement The material point element whose total energy is to be computed
     * @return The total energy of the input material point element
     */
    static double CalculateTotalEnergy(Element& rElement);

    /**
     * @brief Compute the total energy of a model part
     * @detail The total energy is computed summing the kinetic, strain and potential energy
     * @param rModelPart The model part whose total energy is to be computed
     * @return The total energy of the input model part
     */
    static double CalculateTotalEnergy(ModelPart& rModelPart);

    /**
     * @brief Compute the kinetic, potential, strain and total energy of a model part
     * @param rModelPart The model part whose total energy is to be computed
     * @param rPotentialEnergy The potential energy of the input model part
     * @param rKineticEnergy The kinetic energy of the input model part
     * @param rStrainEnergy The strain energy of the input model part
     * @param rTotallEnergy The total energy of the input model part
     */
    static void CalculateAllEnergies(
        ModelPart& rModelPart,
        double& rPotentialEnergy,
        double& rKineticEnergy,
        double& rStrainEnergy,
        double& rTotalEnergy
    );

};

} // end namespace Kratos
