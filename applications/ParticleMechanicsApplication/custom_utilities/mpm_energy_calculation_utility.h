//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_ENERGY_CALCULATION_UTILITY
#define KRATOS_MPM_ENERGY_CALCULATION_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "particle_mechanics_application_variables.h"
#include "containers/model.h"
#include "includes/element.h"

namespace Kratos
{
namespace MPMEnergyCalculationUtility
{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef Element::Pointer ElementPointerType;

    //compute potential energy
    double CalculatePotentialEnergy(Element& rElement);
    double CalculatePotentialEnergy(ModelPart& rModelPart);

     //compute kinetic energy
    double CalculateKineticEnergy(Element& rElement);
    double CalculateKineticEnergy(ModelPart& rModelPart);

     //compute strain energy
    double CalculateStrainEnergy(Element& rElement);
    double CalculateStrainEnergy(ModelPart& rModelPart);

    //compute total energy
    double CalculateTotalEnergy(Element& rElement);
    double CalculateTotalEnergy(ModelPart& rModelPart);

} // end namespace MPMEnergyCalculationUtility

} // end namespace Kratos

#endif // KRATOS_MPM_ENERGY_CALCULATION_UTILITY


