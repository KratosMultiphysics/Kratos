// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
// 					 Navaneeth K Narayanan
//					 Rishith Ellath Meethal
// ==============================================================================
//

// System includes

// External includes

// Project includes

// Application includes
#include "apply_chimera_process_fractional_step.h"

namespace Kratos
{


template <int TDim>
void ApplyChimeraProcessFractionalStep<TDim>::ApplyChimeraProcessFractionalStep(ModelPart &rMainModelPart, Parameters iParameters) : BaseType(rMainModelPart, iParameters)
{

}


template <int TDim>
void ApplyChimeraProcessFractionalStep<TDim>::ApplyContinuityWithElement(Geometry<Node<3>> &rGeometry,
                                                                     Node<3> &rBoundaryNode,
                                                                     Vector &rShapeFuncWeights,
                                                                     unsigned int StartId,
                                                                     std::vector<int> &ConstraintIdVector,
                                                                     typename BaseType::MasterSlaveConstraintContainerType &rMsContainer)
{
    const auto &r_clone_constraint = (LinearMasterSlaveConstraint)KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
    // Initialise the boundary nodes dofs to 0 at ever time steps
    rBoundaryNode.FastGetSolutionStepValue(VELOCITY_X, 0) = 0.0;
    rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Y, 0) = 0.0;
    if (TDim == 3)
        rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Z, 0) = 0.0;
    rBoundaryNode.FastGetSolutionStepValue(PRESSURE, 0) = 0.0;

    for (std::size_t i = 0; i < rGeometry.size(); i++)
    {
        //Interpolation of velocity
        rBoundaryNode.FastGetSolutionStepValue(VELOCITY_X, 0) += rGeometry[i].GetDof(VELOCITY_X).GetSolutionStepValue(0) * rShapeFuncWeights[i];
        rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Y, 0) += rGeometry[i].GetDof(VELOCITY_Y).GetSolutionStepValue(0) * rShapeFuncWeights[i];
        //Interpolation of pressure
        rBoundaryNode.FastGetSolutionStepValue(PRESSURE, 0) += rGeometry[i].GetDof(PRESSURE).GetSolutionStepValue(0) * rShapeFuncWeights[i];

        //Define master slave relation for velocity X and Y
        AddMasterSlaveRelationVelocity(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], VELOCITY_X, rBoundaryNode, VELOCITY_X, rShapeFuncWeights[i]);
        AddMasterSlaveRelationVelocity(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], VELOCITY_Y, rBoundaryNode, VELOCITY_Y, rShapeFuncWeights[i]);
        if (TDim == 3)
        {
            //Interpolation of velocity Z
            rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Z, 0) += rGeometry[i].GetDof(VELOCITY_Z).GetSolutionStepValue(0) * rShapeFuncWeights[i];
            AddMasterSlaveRelationVelocity(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], VELOCITY_Z, rBoundaryNode, VELOCITY_Z, rShapeFuncWeights[i]);
        }
        //Defining master slave relation for pressure
        AddMasterSlaveRelationPressure(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], PRESSURE, rBoundaryNode, PRESSURE, rShapeFuncWeights[i]);
    } // end of loop over host element nodes

    // Setting the buffer 1 same buffer 0
    rBoundaryNode.FastGetSolutionStepValue(VELOCITY_X, 1) = rBoundaryNode.FastGetSolutionStepValue(VELOCITY_X, 0);
    rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Y, 1) = rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Y, 0);
    if (TDim == 3)
        rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Z, 1) = rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Z, 0);
    rBoundaryNode.FastGetSolutionStepValue(PRESSURE, 1) = rBoundaryNode.FastGetSolutionStepValue(PRESSURE, 0);
}


} // namespace Kratos.
