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


// System includes

// External includes

// Project includes
#include "utilities/reduction_utilities.h"
#include "mpm_energy_calculation_utility.h"

namespace Kratos
{
    double MPMEnergyCalculationUtility::CalculatePotentialEnergy(Element& rElement)
    {
        const ProcessInfo& process_info = ProcessInfo();

        double mp_potential_energy = 0.0;
        std::vector<double> mp_mass(1);
        rElement.CalculateOnIntegrationPoints(MP_MASS, mp_mass, process_info);

        std::vector<array_1d<double, 3>> mp_volume_acceleration = { ZeroVector(3) };
        rElement.CalculateOnIntegrationPoints(MP_VOLUME_ACCELERATION, mp_volume_acceleration, process_info);

        std::vector<array_1d<double, 3>> mp_coord = { ZeroVector(3) };
        rElement.CalculateOnIntegrationPoints(MP_COORD, mp_coord, process_info);

        for(unsigned int k = 0; k<3; k++)
            mp_potential_energy += mp_mass[0] * std::abs(mp_volume_acceleration[0][k]) * mp_coord[0][k];

        return mp_potential_energy;
    }

    double MPMEnergyCalculationUtility::CalculatePotentialEnergy(ModelPart& rModelPart)
    {
        return block_for_each<SumReduction<double>>(rModelPart.Elements(),
                [&](auto& r_element) { return this->CalculatePotentialEnergy(r_element); });
    }

    double MPMEnergyCalculationUtility::CalculateKineticEnergy(Element& rElement)
    {
        const ProcessInfo& process_info = ProcessInfo();

        double MP_kinetic_energy = 0.0;
        std::vector<double> mp_mass(1);
        rElement.CalculateOnIntegrationPoints(MP_MASS, mp_mass, process_info);

        std::vector<array_1d<double, 3>> mp_velocity = { ZeroVector(3) };
        rElement.CalculateOnIntegrationPoints(MP_VELOCITY, mp_velocity, process_info);

        for(SizeType k = 0; k<3; ++k)
            MP_kinetic_energy   += 0.5 * mp_mass[0] * mp_velocity[0][k] * mp_velocity[0][k];

        return MP_kinetic_energy;
    }

    double MPMEnergyCalculationUtility::CalculateKineticEnergy(ModelPart& rModelPart)
    {
        return block_for_each<SumReduction<double>>(rModelPart.Elements(),
                [&](auto& r_element) { return this->CalculateKineticEnergy(r_element); });
    }

    double MPMEnergyCalculationUtility::CalculateStrainEnergy(Element& rElement)
    {
        const ProcessInfo& process_info = ProcessInfo();

        double mp_strain_energy = 0.0;
        std::vector<double> mp_volume(1);
        rElement.CalculateOnIntegrationPoints(MP_VOLUME, mp_volume, process_info);

        std::vector<Vector> mp_cauchy_stress(1);
        rElement.CalculateOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, mp_cauchy_stress, process_info);

        std::vector<Vector> mp_almansi_strain(1);
        rElement.CalculateOnIntegrationPoints(MP_ALMANSI_STRAIN_VECTOR, mp_almansi_strain, process_info);

        for(SizeType j = 0; j < mp_cauchy_stress[0].size(); ++j)
        {
            mp_strain_energy +=  0.5 * mp_volume[0] * mp_cauchy_stress[0][j] * mp_almansi_strain[0][j];
        }

        return mp_strain_energy;
    }

    double MPMEnergyCalculationUtility::CalculateStrainEnergy(ModelPart& rModelPart)
    {
        return block_for_each<SumReduction<double>>(rModelPart.Elements(),
                [&](auto& r_element) { return this->CalculateStrainEnergy(r_element); });
    }

    double MPMEnergyCalculationUtility::CalculateTotalEnergy(Element& rElement)
    {
        const double r_MP_potential_energy = CalculatePotentialEnergy(rElement);
        const double r_MP_kinetic_energy   = CalculateKineticEnergy(rElement);
        const double r_MP_strain_energy    = CalculateStrainEnergy(rElement);

        return (r_MP_potential_energy + r_MP_kinetic_energy + r_MP_strain_energy);
    }

    double MPMEnergyCalculationUtility::CalculateTotalEnergy(ModelPart& rModelPart)
    {
        return block_for_each<SumReduction<double>>(rModelPart.Elements(),
                [&](auto& r_element) { return this->CalculateTotalEnergy(r_element); });
    }

} // end namespace Kratos
