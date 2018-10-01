// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

// System includes

// External includes

// Project includes
#include "nodal_displacement_response_function.h"

namespace Kratos {

NodalDisplacementResponseFunction::NodalDisplacementResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
: mrModelPart(rModelPart), mResponseSettings(ResponseSettings)
{
        // Get id of node where a displacement should be traced
        const int id_traced_node = ResponseSettings["traced_node_id"].GetInt();
        // Get pointer to traced node
        mpTracedNode = rModelPart.pGetNode(id_traced_node);

        // Get the corresponding dof to the displacement which should be traced
        // by this response function e.g. DISPLACEMENT_X, ROTATION_X,...
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();

        // Check if variable for traced dof is valid
        if( !( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mTracedDofLabel)) )
            KRATOS_ERROR << "Specified traced DOF is not available. Specified DOF: " << mTracedDofLabel << std::endl;
        else
        {
            const VariableComponentType& r_traced_dof =
                KratosComponents<VariableComponentType>::Get(mTracedDofLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(r_traced_dof) )
                << "Specified DOF is not available at traced node." << std::endl;
        }

    // Check if there are primal elements, because the primal state is required
    ProcessInfo &r_current_process_info = mrModelPart.GetProcessInfo();
    KRATOS_ERROR_IF( r_current_process_info.Has(IS_ADJOINT) && r_current_process_info[IS_ADJOINT] )
            << "NodalDisplacementResponseFunction: Can not use adjoint model part!" << std::endl;
}

double NodalDisplacementResponseFunction::CalculateValue()
{
    KRATOS_TRY;

    const VariableComponentType& r_traced_dof =
        KratosComponents<VariableComponentType>::Get(mTracedDofLabel);

    return mpTracedNode->FastGetSolutionStepValue(r_traced_dof, 0);

    KRATOS_CATCH("");
}


} // namespace Kratos.

