//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andi Makarim Katili
//

// #define KRATOS_MPM_BASE_MAPPING_UTILITIES

#pragma once

// Project includes
#include "custom_utilities/mapping_utilities/mpm_base_particle_mapping_utility.hpp"
#include "utilities/atomic_utilities.h"



namespace Kratos
{
    /**
     * @class MPMPicParticleMappingUtility
     * @ingroup KratosMPM
     * @brief PIC particle mapping utility for Material Point Method
     * @details ToDo: details
     * continuation of details
     */
    class MPMPicParticleMappingUtility
        : public MPMFlipParticleMappingUtility
    {

    public:

    KRATOS_CLASS_POINTER_DEFINITION(MPMPicParticleMappingUtility);
    using IndexType = std::size_t;
    using MappingBaseType = MPMBaseParticleMappingUtility;

    ///@name Life Cycle
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * constructors, copy constructors and destructor: MANDATORY
     */

    /**
     * Constructor.
     */
    MPMPicParticleMappingUtility(ModelPart& rMPMModelPart, ModelPart& rGridModelPart, int EchoLevel)
        : MPMFlipParticleMappingUtility(rMPMModelPart, rGridModelPart, EchoLevel)
    {
    }

    ///@}
    ///@name Operations
    ///@{


    #pragma region Particle to Grid Mapping (P2G)

    // Nodal momentum and inertia is the same as FLIP

    #pragma endregion
    // end of P2G Mapping

    #pragma region Grid to Particle Mapping (G2P)

    void G2PVelocity(Element& rElement, const array_1d<double, 3>& rNewMPAcceleration) override
    {
        array_1d<double,3> new_mp_velocity{};
        MappingBaseType::EvaluateVariableOnMaterialPoint(rElement, VELOCITY, new_mp_velocity);
        rElement.SetValuesOnIntegrationPoints(MP_VELOCITY, {new_mp_velocity}, mrProcessInfo);
    }

    #pragma endregion
    // end of G2P Mapping

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

};// end of MPMBaseParticleMappingUtility
}// end of kratos namespace