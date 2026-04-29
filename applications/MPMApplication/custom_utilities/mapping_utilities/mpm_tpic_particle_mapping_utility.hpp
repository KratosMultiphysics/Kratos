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
    class MPMTpicParticleMappingUtility
        : public MPMPicParticleMappingUtility
    {

    public:

    KRATOS_CLASS_POINTER_DEFINITION(MPMTpicParticleMappingUtility);
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
    MPMTpicParticleMappingUtility(ModelPart& rMPMModelPart, ModelPart& rGridModelPart, int EchoLevel)
        : MPMPicParticleMappingUtility(rMPMModelPart, rGridModelPart, EchoLevel)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * is called to initialize the mapping scheme
     * if the particle mapping scheme  needs to perform any operation before any calculation is done
     * the quantities will be initialized and set using this method
     */
    void Initialize() override
    {

        for (Element& rElement : mrMPMModelPart.Elements())
        {
            const unsigned int dimension = rElement.GetGeometry().WorkingSpaceDimension();
            const Matrix InitiateGradient = ZeroMatrix(dimension, dimension);

            // Initialize gradients to zero (not suitable for simulation with prescribed velocities at MP)
            rElement.SetValuesOnIntegrationPoints(MP_VELOCITY_GRADIENT, {InitiateGradient}, mrProcessInfo);
            rElement.SetValuesOnIntegrationPoints(MP_ACCELERATION_GRADIENT, {InitiateGradient}, mrProcessInfo);
        }
    }

    #pragma region Particle to Grid Mapping (P2G)
    void CalculateMpToNodeVector(Element& rElement, Node& rNode, array_1d<double,3>& rMpToNodeVector)
    {
        std::vector<array_1d<double,3>> mp_coordinate{};
        rElement.CalculateOnIntegrationPoints(MP_COORD, mp_coordinate, mrProcessInfo);

        rMpToNodeVector = rNode.Coordinates() - mp_coordinate[0];
    }
    /**
     * Calculate and add TPIC affine velocity for the corresponding node
     */
    void CalculateAndAddAffineVelocity(Element& rElement, Node& rNode, std::vector<array_1d<double, 3 >> rNodalAffineVelocity)
    {
        array_1d<double,3> mp_to_node_vector = ZeroVector(3);
        this->CalculateMpToNodeVector(rElement, rNode, mp_to_node_vector);



        std::vector<Matrix> mp_velocity_gradient;
        rElement.CalculateOnIntegrationPoints(MP_VELOCITY_GRADIENT, mp_velocity_gradient, mrProcessInfo);
        // KRATOS_WATCH(rElement.Id())
        // KRATOS_WATCH(rNode.Id())
        // KRATOS_WATCH(mp_to_node_vector)
        // KRATOS_WATCH(mp_velocity_gradient[0])
        // KRATOS_WATCH(rNodalAffineVelocity[0])
        rNodalAffineVelocity[0] += prod(mp_velocity_gradient[0], mp_to_node_vector);
        // KRATOS_WATCH(rNodalAffineVelocity[0])

    }

    /**
     * Calculate and add TPIC affine acceleration for the corresponding node
     */
    void CalculateAndAddAffineAcceleration(Element& rElement, Node& rNode, std::vector<array_1d<double, 3 >> rNodalAffineAcceleration)
    {
        array_1d<double,3> mp_to_node_vector = ZeroVector(3);
        this->CalculateMpToNodeVector(rElement, rNode, mp_to_node_vector);

        std::vector<Matrix> mp_acceleration_gradient;
        rElement.CalculateOnIntegrationPoints(MP_ACCELERATION_GRADIENT, mp_acceleration_gradient, mrProcessInfo);

        rNodalAffineAcceleration[0] += prod(mp_acceleration_gradient[0], mp_to_node_vector);
    }

    /**
     * Do Particle to Grid mapping for nodal momentum.
     */
    void P2GMomentum(Element& rElement, Node& rNode, const double& rN_i) override
    {
        std::vector<array_1d<double, 3 >> node_affine_velocity;
        std::vector<double> mp_mass;

        // node_affine_velocity = mp_velocity + affine contribution
        rElement.CalculateOnIntegrationPoints(MP_VELOCITY, node_affine_velocity, mrProcessInfo);
        this->CalculateAndAddAffineVelocity(rElement, rNode, node_affine_velocity);

        rElement.CalculateOnIntegrationPoints(MP_MASS, mp_mass, mrProcessInfo);

        array_1d<double,3> r_nodal_momentum = node_affine_velocity[0] * rN_i * mp_mass[0];
        AtomicAdd(rNode.FastGetSolutionStepValue(NODAL_MOMENTUM, 0), r_nodal_momentum);
    }

    /**
     * Do Particle to Grid mapping for nodal inertia.
     */
    void P2GInertia(Element& rElement, Node& rNode, const double& rN_i) override
    {
        std::vector<array_1d<double, 3 >> node_affine_acceleration;
        std::vector<double> mp_mass;

        // node_affine_acceleration = mp_acceleration + affine contribution
        rElement.CalculateOnIntegrationPoints(MP_ACCELERATION, node_affine_acceleration, mrProcessInfo);
        this->CalculateAndAddAffineAcceleration(rElement, rNode, node_affine_acceleration);

        rElement.CalculateOnIntegrationPoints(MP_MASS, mp_mass, mrProcessInfo);

        array_1d<double,3> r_nodal_inertia = node_affine_acceleration[0] * rN_i * mp_mass[0];

        AtomicAdd(rNode.FastGetSolutionStepValue(NODAL_INERTIA, 0), r_nodal_inertia);
    }

    #pragma endregion
    // end of P2G Mapping

    #pragma region Grid to Particle Mapping (G2P)

    void CalculateMPVelocityGradient(Element& rElement, const Matrix& dN_dX)
    {
        const unsigned int dimension = rElement.GetGeometry().WorkingSpaceDimension();
        Matrix mp_velocity_gradient = ZeroMatrix(dimension, dimension);

        IndexType node_index = 0;
        for (const Node& node_i : rElement.GetGeometry())
        {
            const array_1d<double,3> velocity_i = node_i.FastGetSolutionStepValue(VELOCITY, 0);
            const Vector& dN_dX_i = row(dN_dX, node_index);

            // calculate velocity gradient for current node
            for ( IndexType dim_i = 0; dim_i < dimension; dim_i++ )
            {
                for ( IndexType dim_j = 0; dim_j < dimension; dim_j++ )
                {
                    mp_velocity_gradient(dim_i,dim_j) += velocity_i[dim_i] * dN_dX_i(dim_j);
                }
            }

            node_index++;
        }
        rElement.SetValuesOnIntegrationPoints(MP_VELOCITY_GRADIENT, {mp_velocity_gradient}, mrProcessInfo);
    }

    void CalculateMPAccelerationGradient(Element& rElement, const Matrix& dN_dX)
    {
        const unsigned int dimension = rElement.GetGeometry().WorkingSpaceDimension();
        Matrix mp_acceleration_gradient = ZeroMatrix(dimension, dimension);

        IndexType node_index = 0;
        for (const Node& node_i : rElement.GetGeometry())
        {
            const array_1d<double,3> acceleration_i = node_i.FastGetSolutionStepValue(ACCELERATION, 0);
            const Vector& dN_dX_i = row(dN_dX, node_index);

            // calculate acceleration gradient for current node
            for ( IndexType dim_i = 0; dim_i < dimension; dim_i++ )
            {
                for ( IndexType dim_j = 0; dim_j < dimension; dim_j++ )
                {
                    mp_acceleration_gradient(dim_i,dim_j) += acceleration_i[dim_i] * dN_dX_i(dim_j);
                }
            }

            node_index++;
        }

        rElement.SetValuesOnIntegrationPoints(MP_ACCELERATION_GRADIENT, {mp_acceleration_gradient}, mrProcessInfo);
    }


    void G2PAdditionalVariables(Element& rElement) override
    {
        MappingBaseType::G2PAdditionalVariables(rElement);

        Matrix Jacobian;
        rElement.GetGeometry().Jacobian(Jacobian, 0);
        Matrix InvJ;
        double detJ;
        MathUtils<double>::InvertMatrix(Jacobian, InvJ, detJ);
        Matrix r_dN_de = rElement.GetGeometry().ShapeFunctionLocalGradient(0);
        Matrix dN_dX = prod(r_dN_de, InvJ); // cartesian gradients
        // Calculate velocity and acceleration derivatives
        this->CalculateMPVelocityGradient(rElement, dN_dX);
        this->CalculateMPAccelerationGradient(rElement, dN_dX);
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