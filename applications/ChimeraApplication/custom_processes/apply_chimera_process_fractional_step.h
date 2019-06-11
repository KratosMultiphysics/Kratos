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
#if !defined(KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED)
#define KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "apply_chimera_process_monolithic.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <int TDim>
class ApplyChimeraProcessFractionalStep : public ApplyChimeraProcessMonolithic<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ApplyChimeraProcessFractionalStep
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessFractionalStep);
    typedef ApplyChimeraProcessMonolithic<TDim> BaseType;
    typedef typename BaseType::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;

    ///@}
    ///@name Life Cycle
    ///@{
    ApplyChimeraProcessFractionalStep(ModelPart &rMainModelPart, Parameters iParameters)
        : BaseType(rMainModelPart, iParameters)
    {
    }


    /// Destructor.
    virtual ~ApplyChimeraProcessFractionalStep()
    {
    }
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

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
    /**
     * @brief Applies the master-slave constraint to enforce the continuity between a given geometry/element and a boundary node
     * @param rGeometry The geometry of the element
     * @param rBoundaryNode The boundary node for which the connections are to be made.
     * @param rShapeFuncWeights The shape function weights for the node in the rGeometry.
     * @param StartIndex The start Index of the constraints which are to be added.
     * @param rConstraintIdVector The vector of the constraints Ids which is accessed with StartIndex.
     * @param rMsContainer The Constraint container to which the contraints are added.
     */
    void ApplyContinuityWithElement(Geometry<Node<3>> &rGeometry,
                                    Node<3> &rBoundaryNode,
                                    Vector &rShapeFuncWeights,
                                    unsigned int StartId,
                                    std::vector<int> &ConstraintIdVector,
                                    typename BaseType::MasterSlaveConstraintContainerType &rMsContainer) override
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Applies the master-slave constraint between the given master and slave nodes with corresponding variable.
     *          In addition it sets the FS_CHIMERA_VEL_CONSTRAINT flag on the constraint which is used in the FS-Strategy for chimera
     * @param rMasterSlaveContainer The container to which the constraint to be added (useful to so OpenMP loop)
     * @param rCloneConstraint The prototype of constraint which is to be added.
     * @param ConstraintId The ID of the constraint to be added.
     * @param rMasterNode The Master node of the constraint.
     * @param rMasterVariable The variable for the master node.
     * @param rSlaveNode The Slave node of the constraint.
     * @param rSlaveVariable The variable for the slave node.
     * @param Weight The weight of the Master node.
     * @param Constant The constant of the master slave relation.
     */
    template <typename TVariableType>
    inline void AddMasterSlaveRelationVelocity(typename BaseType::MasterSlaveConstraintContainerType &rMasterSlaveContainer,
                                       const LinearMasterSlaveConstraint &rCloneConstraint,
                                       unsigned int ConstraintId,
                                       Node<3> &rMasterNode,
                                       TVariableType &rMasterVariable,
                                       Node<3> &rSlaveNode,
                                       TVariableType &rSlaveVariable,
                                       const double Weight,
                                       const double Constant = 0.0)
    {
        rSlaveNode.Set(SLAVE);
        ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = rCloneConstraint.Create(ConstraintId, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        p_new_constraint->Set(TO_ERASE);
        BaseType::mNodeIdToConstraintIdsMap[rSlaveNode.Id()].push_back(ConstraintId);
        rMasterSlaveContainer.insert(rMasterSlaveContainer.begin(), p_new_constraint);

        // TODO: Set the FS_CHIMERA_VEL_CONSTRAINT variable to true
        p_new_constraint->Set(FS_CHIMERA_VEL_CONSTRAINT);
    }

    /**
     * @brief Applies the master-slave constraint between the given master and slave nodes with corresponding variable.
     *          In addition it sets the FS_CHIMERA_PRE_CONSTRAINT flag on the constraint which is used in the FS-Strategy for chimera
     * @param rMasterSlaveContainer The container to which the constraint to be added (useful to so OpenMP loop)
     * @param rCloneConstraint The prototype of constraint which is to be added.
     * @param ConstraintId The ID of the constraint to be added.
     * @param rMasterNode The Master node of the constraint.
     * @param rMasterVariable The variable for the master node.
     * @param rSlaveNode The Slave node of the constraint.
     * @param rSlaveVariable The variable for the slave node.
     * @param Weight The weight of the Master node.
     * @param Constant The constant of the master slave relation.
     */
    template <typename TVariableType>
    inline void AddMasterSlaveRelationPressure(typename BaseType::MasterSlaveConstraintContainerType &rMasterSlaveContainer,
                                       const LinearMasterSlaveConstraint &rCloneConstraint,
                                       unsigned int ConstraintId,
                                       Node<3> &rMasterNode,
                                       TVariableType &rMasterVariable,
                                       Node<3> &rSlaveNode,
                                       TVariableType &rSlaveVariable,
                                       const double Weight,
                                       const double Constant = 0.0)
    {
        rSlaveNode.Set(SLAVE);
        ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = rCloneConstraint.Create(ConstraintId, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        p_new_constraint->Set(TO_ERASE);
        BaseType::mNodeIdToConstraintIdsMap[rSlaveNode.Id()].push_back(ConstraintId);
        rMasterSlaveContainer.insert(rMasterSlaveContainer.begin(), p_new_constraint);

        // TODO: Set the FS_CHIMERA_PRE_CONSTRAINT variable to true
        p_new_constraint->Set(FS_CHIMERA_PRE_CONSTRAINT);
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ApplyChimeraProcessFractionalStep &operator=(ApplyChimeraProcessFractionalStep const &rOther);

    ///@}




}; // Class ApplyChimeraProcessFractionalStep

} // namespace Kratos.

#endif // KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED
