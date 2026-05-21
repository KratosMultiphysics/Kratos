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
#include "mpm_energy_calculation_utility.h"

namespace Kratos
{
namespace MPMEnergyCalculationUtility
{
    /**
     * @brief Assign and compute potential energy*
     */
    double CalculatePotentialEnergy(Element& rElement)
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

    double CalculatePotentialEnergy(ModelPart& rModelPart)
    {
        double model_part_potential_energy = 0.0;
        for(SizeType i = 0; i < rModelPart.Elements().size(); ++i) {

            auto element_itr = rModelPart.Elements().begin() + i;
            model_part_potential_energy += CalculatePotentialEnergy(**(element_itr.base()));
        }
        return model_part_potential_energy;
    }

     //Assign and compute kinetic energy
    double CalculateKineticEnergy(Element& rElement)
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

    double CalculateKineticEnergy(ModelPart& rModelPart)
    {
        double model_part_kinetic_energy = 0.0;
        for(SizeType i = 0; i < rModelPart.Elements().size(); ++i){

            auto element_itr = rModelPart.Elements().begin() + i;
            model_part_kinetic_energy += CalculateKineticEnergy(**(element_itr.base()));
        }
        return model_part_kinetic_energy;
    }

     //Assign and compute strain energy
    double CalculateStrainEnergy(Element& rElement)
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

    double CalculateStrainEnergy(ModelPart& rModelPart)
    {
        double model_part_strain_energy = 0.0;
        for(SizeType i = 0; i < rModelPart.Elements().size(); ++i){

            auto element_itr = rModelPart.Elements().begin() + i;
            model_part_strain_energy += CalculateStrainEnergy(**(element_itr.base()));
        }
        return model_part_strain_energy;
    }

    /**
     * @brief Assign and compute total energy
     * @details Compute total energy inside material point, summing potential, kinetic and strain energy
     *
     */
    double CalculateTotalEnergy(Element& rElement)
    {
        const double r_MP_potential_energy = CalculatePotentialEnergy(rElement);
        const double r_MP_kinetic_energy   = CalculateKineticEnergy(rElement);
        const double r_MP_strain_energy    = CalculateStrainEnergy(rElement);

        return (r_MP_potential_energy + r_MP_kinetic_energy + r_MP_strain_energy);
    }


    double CalculateTotalEnergy(ModelPart& rModelPart)
    {
        double model_part_total_energy = 0.0;
        for(SizeType i = 0; i < rModelPart.Elements().size(); ++i){

            auto element_itr = rModelPart.Elements().begin() + i;
            model_part_total_energy += CalculateTotalEnergy(**(element_itr.base()));
        }
        return model_part_total_energy;
    }

} // end namespace MPMEnergyCalculationUtility

} // end namespace Kratos


