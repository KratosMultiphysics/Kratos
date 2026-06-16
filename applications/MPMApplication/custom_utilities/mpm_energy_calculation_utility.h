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

class KRATOS_API(MPM_APPLICATION) MPMEnergyCalculationUtility
{
public:

    /**
     * @brief Compute the potential energy of a material point element
     * @param rElement The material point element whose potential energy is to be computed
     * @return The potential energy of the input material point element
     */
    static double CalculatePotentialEnergy(Element& rElement);

    /**
     * @brief Compute the sum of potential energy of material points in a model part
     * @param rMpmModelPart The model part whose potential energy is to be computed
     * @return The potential energy of the input model part
     */
    static double CalculatePotentialEnergy(ModelPart& rMpmModelPart);

    /**
     * @brief Compute the kinetic energy of a material point element
     * @param rElement The material point element whose kinetic energy is to be computed
     * @return The kinetic energy of the input material point element
     */
    static double CalculateKineticEnergy(Element& rElement);

    /**
     * @brief Compute the sum of  kinetic energy of material points in a model part   
     * @param rMpmModelPart The model part whose kinetic energy is to be computed
     * @return The kinetic energy of the input model part
     */
    static double CalculateKineticEnergy(ModelPart& rMpmModelPart);

    /**
     * @brief Compute the strain energy of a material point element
     * @param rElement The material point element whose strain energy is to be computed
     * @return The strain energy of the input material point element
     */
    static double CalculateStrainEnergy(Element& rElement);

    /**
     * @brief Compute the sum of strain energy of material points in a model part

     * @param rMpmModelPart The model part whose strain energy is to be computed
     * @return The strain energy of the input model part
     */
    static double CalculateStrainEnergy(ModelPart& rMpmModelPart);

    /**
     * @brief Compute the kinetic, potential, strain and total energy of a material point element
     * @param rElement The element whose total energy is to be computed
     * @return rPotentialEnergy The potential energy of the input element
     * @return rKineticEnergy The kinetic energy of the input element
     * @return rStrainEnergy The strain energy of the input element
     * @return rTotallEnergy The total energy of the input element
     */
    static std::tuple<double,double,double,double> CalculateAllEnergies(Element& rElement);

    /**
     * @brief Compute the total kinetic, potential, strain and total energy of material points in a model part
     * @param rMpmModelPart The model part whose total energy is to be computed
     * @param rPotentialEnergy The potential energy of the input model part
     * @param rKineticEnergy The kinetic energy of the input model part
     * @param rStrainEnergy The strain energy of the input model part
     * @param rTotallEnergy The total energy of the input model part
     */
    static std::tuple<double,double,double,double> CalculateAllEnergies(ModelPart& rMpmModelPart);

};

} // end namespace Kratos
