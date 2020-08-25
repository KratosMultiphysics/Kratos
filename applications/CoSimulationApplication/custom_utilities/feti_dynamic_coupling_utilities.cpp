//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Big T
//

// System includes

// External includes

// Project includes
#include "feti_dynamic_coupling_utilities.h"

namespace Kratos
{

//void FetiDynamicCouplingUtilities::FindIntersection1DGeometries2D(


	void FetiDynamicCouplingUtilities::EquilibrateDomains()
	{
        // calculate unbalanced interface free velocity
        const SizeType dim = mrOriginModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const SizeType origin_interface_dofs = dim* mrOriginModelPart.NumberOfNodes();
        Vector unbalanced_interface_free_velocity(origin_interface_dofs);

        KRATOS_ERROR_IF_NOT(origin_interface_dofs == dim * mrDestinationModelPart.NumberOfNodes())
            << "Origin and Destination have different number of interface DOFS.";

        CalculateUnbalancedInterfaceFreeVelocities(unbalanced_interface_free_velocity);


        // calculate sensitivity accelerations
        Vector origin_response_accel(mpKOrigin->size1());
        Vector destination_response_accel(mpKDestination->size1());


        // calculate condensation matrix
        // calculate lagrange mults
        // calculate correction velocities
        // apply correct

	}

	void FetiDynamicCouplingUtilities::CalculateUnbalancedInterfaceFreeVelocities(Vector& rUnbalancedVelocities)
	{
        // TODO - the nodal ordering in the model parts will generally not coincide.
        // TODO - we need to make sure we are subtracting the correct node from the other
        auto origin_interface_nodes = mrOriginModelPart.NodesArray();
        auto destination_interface_nodes = mrDestinationModelPart.NodesArray();
        const SizeType dim = mrOriginModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        // TODO this will need to be mortar mapped in the future
        Vector mapped_destination_interface_velocities(origin_interface_nodes.size() * dim);
        for (size_t i = 0; i < destination_interface_nodes.size(); ++i)
        {
            mapped_destination_interface_velocities[dim * i] =
                destination_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_X);
            mapped_destination_interface_velocities[dim * i + 1] =
                destination_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_Y);
            if (dim == 3) mapped_destination_interface_velocities[dim * i + 2] =
                destination_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_Z);
        }

        // Fill unbalanced velocities with origin velocities
        for (size_t i = 0; i < origin_interface_nodes.size(); ++i)
        {
            rUnbalancedVelocities[dim*i] = origin_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_X);
            rUnbalancedVelocities[dim*i + 1] = origin_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_Y);
            if(dim == 3) rUnbalancedVelocities[dim*i + 2] =
                origin_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_Z);
        }

        // Subtract mapped destination velocities
        KRATOS_ERROR_IF_NOT(mapped_destination_interface_velocities.size() == rUnbalancedVelocities.size())
            << "Mapped destination interface velocities and origin interface velocities must have the same size";
        for (size_t i = 0; i < mapped_destination_interface_velocities.size(); i++)
        {
            rUnbalancedVelocities[i] -= mapped_destination_interface_velocities[i];
        }
	}

} // namespace Kratos.
