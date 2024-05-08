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
#include "mpm_application_variables.h"
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
    double KRATOS_API(MPM_APPLICATION) CalculatePotentialEnergy(Element& rElement);
    double KRATOS_API(MPM_APPLICATION) CalculatePotentialEnergy(ModelPart& rModelPart);

     //compute kinetic energy
    double KRATOS_API(MPM_APPLICATION) CalculateKineticEnergy(Element& rElement);
    double KRATOS_API(MPM_APPLICATION) CalculateKineticEnergy(ModelPart& rModelPart);

     //compute strain energy
    double KRATOS_API(MPM_APPLICATION) CalculateStrainEnergy(Element& rElement);
    double KRATOS_API(MPM_APPLICATION) CalculateStrainEnergy(ModelPart& rModelPart);

    //compute total energy
    double KRATOS_API(MPM_APPLICATION) CalculateTotalEnergy(Element& rElement);
    double KRATOS_API(MPM_APPLICATION) CalculateTotalEnergy(ModelPart& rModelPart);

} // end namespace MPMEnergyCalculationUtility

} // end namespace Kratos

#endif // KRATOS_MPM_ENERGY_CALCULATION_UTILITY


