//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// ==============================================================================
//
#if !defined(KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED) 
#define KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "apply_chimera_process.h"

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
class KRATOS_API(CHIMERA_APPLICATION) ApplyChimeraProcessFractionalStep : public ApplyChimera<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ApplyChimeraProcessFractionalStep
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessFractionalStep);
    typedef ApplyChimera<TDim> BaseType;
    typedef typename BaseType::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef typename BaseType::PointLocatorType PointLocatorType;
    typedef typename BaseType::PointLocatorPointerType PointLocatorPointerType;
    typedef typename BaseType::MasterSlaveContainerVectorType MasterSlaveContainerVectorType;
    typedef typename BaseType::ConstraintIdsVectorType ConstraintIdsVectorType;

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

    virtual void ExecuteFinalizeSolutionStep() override
    {
        if(BaseType::mReformulateEveryStep)
        {
            auto &vel_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_velocity_model_part");
            vel_modelpart.MasterSlaveConstraints().clear();

            auto &pre_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_pressure_model_part");
            pre_modelpart.MasterSlaveConstraints().clear();
        }

        BaseType::ExecuteFinalizeSolutionStep();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    std::string Info() const override
    {
        return "ApplyChimeraProcessFractionalStep";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "ApplyChimeraProcessFractionalStep" << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        KRATOS_INFO("ApplyChimeraProcessFractionalStep") << std::endl;
    }

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
     * @brief Applies the continuity between the boundary modelpart and the background.
     * @param rBoundaryModelPart The boundary modelpart for which the continuity is to be enforced.
     * @param rBinLocator The bin based locator formulated on the background. This is used to locate nodes on rBoundaryModelPart.
     */
    void ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorType &rBinLocator) override
    {
        //loop over nodes and find the triangle in which it falls, then do interpolation
        MasterSlaveContainerVectorType velocity_ms_container_vector;
        MasterSlaveContainerVectorType pressure_ms_container_vector;
        BaseType::ReserveMemoryForConstraintContainers(rBoundaryModelPart, velocity_ms_container_vector);
        BaseType::ReserveMemoryForConstraintContainers(rBoundaryModelPart, pressure_ms_container_vector);

        const int n_boundary_nodes = static_cast<int> (rBoundaryModelPart.Nodes().size());
        for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
        {
            ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
            Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());

            BaseType::mNodeIdToConstraintIdsMap[p_boundary_node->Id()].reserve(150);
        }

        BaseType::FormulateConstraints(rBoundaryModelPart, rBinLocator, velocity_ms_container_vector, pressure_ms_container_vector);

        BuiltinTimer mpc_add_time;
        auto &vel_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_velocity_model_part");
        BaseType::AddConstraintsToModelpart(vel_modelpart, velocity_ms_container_vector);
        VariableUtils().SetFlag(FS_CHIMERA_PRESSURE_CONSTRAINT, false, vel_modelpart.MasterSlaveConstraints());
        VariableUtils().SetFlag(FS_CHIMERA_VELOCITY_CONSTRAINT, true, vel_modelpart.MasterSlaveConstraints());
        VariableUtils().SetFlag(ACTIVE, true, vel_modelpart.MasterSlaveConstraints());

        auto &pre_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_pressure_model_part");
        BaseType::AddConstraintsToModelpart(pre_modelpart, pressure_ms_container_vector);
        VariableUtils().SetFlag(FS_CHIMERA_PRESSURE_CONSTRAINT, true, vel_modelpart.MasterSlaveConstraints());
        VariableUtils().SetFlag(FS_CHIMERA_VELOCITY_CONSTRAINT, false, vel_modelpart.MasterSlaveConstraints());
        VariableUtils().SetFlag(ACTIVE, true, pre_modelpart.MasterSlaveConstraints());

        KRATOS_INFO_IF("Adding of MPCs from containers to modelpart took : ", BaseType::mEchoLevel > 1)<< mpc_add_time.ElapsedSeconds()<< " seconds"<< std::endl;
        KRATOS_INFO_IF("Number of boundary nodes in : ", BaseType::mEchoLevel > 1) << rBoundaryModelPart.Name() << " : " << rBoundaryModelPart.NumberOfNodes() << std::endl;
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
