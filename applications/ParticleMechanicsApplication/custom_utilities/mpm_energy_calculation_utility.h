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

    /**
     * @brief Assign and compute potential energy*
     */
    void CalculatePotentialEnergy(ElementPointerType& rpElement)
    {
        double MP_PotentialEnergy = 0.0;

        for(unsigned int k = 0; k<3; k++)
            MP_PotentialEnergy += rpElement->GetValue(MP_MASS) * std::abs(rpElement->GetValue(MP_VOLUME_ACCELERATION)[k] * rpElement->GetValue(MP_COORD)[k]);

        rpElement->SetValue(MP_POTENTIAL_ENERGY, MP_PotentialEnergy);
    }

    void CalculatePotentialEnergy(ModelPart& rModelPart)
    {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){

            auto element_itr = rModelPart.Elements().begin() + i;
            CalculatePotentialEnergy(*(element_itr.base()));
        }
    }

    /**
     * @brief Assign and compute kinetic energy*
     */
    void CalculateKineticEnergy(ElementPointerType& rpElement)
    {
        double MP_KineticEnergy = 0.0;

        for(unsigned int k = 0; k<3; k++)
            MP_KineticEnergy   += 0.5 * rpElement->GetValue(MP_MASS) * rpElement->GetValue(MP_VELOCITY)[k] * rpElement->GetValue(MP_VELOCITY)[k] ;

        rpElement->SetValue(MP_KINETIC_ENERGY, MP_KineticEnergy);
    }

    void CalculateKineticEnergy(ModelPart& rModelPart)
    {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){

            auto element_itr = rModelPart.Elements().begin() + i;
            CalculateKineticEnergy(*(element_itr.base()));
        }
    }

    /**
     * @brief Assign and compute strain energy*
     */

    void CalculateStrainEnergy(ElementPointerType& rpElement)
    {
        double MP_StrainEnergy = 0.0;

        for(unsigned int j = 0; j < rpElement->GetValue(MP_CAUCHY_STRESS_VECTOR).size(); j++)
        {
            MP_StrainEnergy +=  0.5 * rpElement->GetValue(MP_VOLUME) * rpElement->GetValue(MP_CAUCHY_STRESS_VECTOR)[j] * rpElement->GetValue(MP_ALMANSI_STRAIN_VECTOR)[j];
        }

        rpElement->SetValue(MP_STRAIN_ENERGY, MP_StrainEnergy);
    }

    void CalculateStrainEnergy(ModelPart& rModelPart)
    {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){

            auto element_itr = rModelPart.Elements().begin() + i;
            CalculateStrainEnergy(*(element_itr.base()));
        }
    }

    /**
     * @brief Assign and compute total energy
     * @details Compute total energy inside particle, summing potential, kinetic and strain energy
     *
     */
    void CalculateTotalEnergy(ElementPointerType& rpElement)
    {
        CalculatePotentialEnergy(rpElement);
        CalculateKineticEnergy(rpElement);
        CalculateStrainEnergy(rpElement);

        const double & MP_PotentialEnergy = rpElement->GetValue(MP_POTENTIAL_ENERGY);
        const double & MP_KineticEnergy   = rpElement->GetValue(MP_KINETIC_ENERGY);
        const double & MP_StrainEnergy    = rpElement->GetValue(MP_STRAIN_ENERGY);
        const double MP_TotalEnergy = MP_PotentialEnergy + MP_KineticEnergy + MP_StrainEnergy;

        rpElement->SetValue(MP_TOTAL_ENERGY, MP_TotalEnergy);
    }


    void CalculateTotalEnergy(ModelPart& rModelPart)
    {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){

            auto element_itr = rModelPart.Elements().begin() + i;
            CalculateTotalEnergy(*(element_itr.base()));
        }
    }

} // end namespace MPMEnergyCalculationUtility

} // end namespace Kratos

#endif // KRATOS_MPM_ENERGY_CALCULATION_UTILITY


