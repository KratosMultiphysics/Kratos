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
     * @class MPMBaseParticleMappingUtility
     * @ingroup KratosMPM
     * @brief Base particle mapping utility for Material Point Method
     * @details ToDo: details
     * continuation of details
     */
    class MPMFlipParticleMappingUtility
        : public MPMBaseParticleMappingUtility
    {
        
    public:
        
    typedef std::size_t IndexType;
    KRATOS_CLASS_POINTER_DEFINITION(MPMFlipParticleMappingUtility);

    ///@name Life Cycle
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * constructors, copy constructors and destructor: MANDATORY
     */

    /**
     * Constructor.
     */
    MPMFlipParticleMappingUtility(ModelPart& rMPMModelPart, ModelPart& rGridModelPart, int EchoLevel)
        : MPMBaseParticleMappingUtility(rMPMModelPart, rGridModelPart, EchoLevel)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    
    #pragma region Particle to Grid Mapping (P2G)
    /**
     * Do Particle to Grid mapping for nodal mass.
     */
    void P2GMomentum(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo) override
    {
        std::vector<array_1d<double, 3 >> mp_velocity;
        std::vector<double> mp_mass;
        
        rElement.CalculateOnIntegrationPoints(MP_VELOCITY, mp_velocity, rCurrentProcessInfo);
        rElement.CalculateOnIntegrationPoints(MP_MASS, mp_mass, rCurrentProcessInfo);

        array_1d<double,3> r_nodal_momentum = mp_velocity[0] * rN_i * mp_mass[0];
        AtomicAdd(rNode.FastGetSolutionStepValue(NODAL_MOMENTUM, 0), r_nodal_momentum);
    }

    /**
     * Do Particle to Grid mapping for nodal inertia.
     */
    void P2GInertia(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo) override
    {
        std::vector<array_1d<double, 3 >> mp_acceleration;
        std::vector<double> mp_mass;

        rElement.CalculateOnIntegrationPoints(MP_ACCELERATION, mp_acceleration, rCurrentProcessInfo);
        rElement.CalculateOnIntegrationPoints(MP_MASS, mp_mass, rCurrentProcessInfo);


        array_1d<double,3> r_nodal_inertia = mp_acceleration[0] * rN_i * mp_mass[0];
        
        AtomicAdd(rNode.FastGetSolutionStepValue(NODAL_INERTIA, 0), r_nodal_inertia);
    }

    #pragma endregion
    // end of P2G Mapping
    
    #pragma region Grid to Particle Mapping (G2P)

    void G2PVelocity(Element& rElement, array_1d<double, 3>& rNewMPAcceleration, const ProcessInfo& rCurrentProcessInfo) override
    {
        // array_1d<double,3> mp_previous_velocity;
        std::vector<array_1d<double, 3 >> mp_new_velocity;
        rElement.CalculateOnIntegrationPoints(MP_VELOCITY, mp_new_velocity, rCurrentProcessInfo);
        std::vector<array_1d<double, 3 >> mp_previous_acceleration;
        rElement.CalculateOnIntegrationPoints(MP_ACCELERATION, mp_previous_acceleration, rCurrentProcessInfo);

        
        /* NOTE:
        Another way to update the MP velocity (see paper Guilkey and Weiss, 2003).
        This assume newmark (or trapezoidal, since n.gamma=0.5) rule of integration*/
        const double gamma_n = 0.5;
        const double delta_time = rCurrentProcessInfo[DELTA_TIME];
        mp_new_velocity[0] += delta_time * (mp_previous_acceleration[0] * (1.0-gamma_n) + gamma_n * rNewMPAcceleration);
        rElement.SetValuesOnIntegrationPoints(MP_VELOCITY, mp_new_velocity, rCurrentProcessInfo);
    }


    void G2PAdditionalVariables(Element& rElement, const ProcessInfo& rCurrentProcessInfo) override
    {
        MPMBaseParticleMappingUtility::G2PAdditionalVariables(rElement, rCurrentProcessInfo);
        // Flip requires no additional variables
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