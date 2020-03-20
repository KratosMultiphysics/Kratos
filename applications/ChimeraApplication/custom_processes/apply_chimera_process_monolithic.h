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

#if !defined(KRATOS_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED)
#define KRATOS_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED

// System includes
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include "omp.h"

// External includes

// Project includes
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
class KRATOS_API(CHIMERA_APPLICATION) ApplyChimeraProcessMonolithic : public ApplyChimera<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessMonolithic);
    typedef ApplyChimera<TDim> BaseType;
    typedef typename BaseType::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef typename BaseType::PointLocatorType PointLocatorType;
    typedef typename BaseType::PointLocatorPointerType PointLocatorPointerType;
    typedef typename BaseType::MasterSlaveContainerVectorType MasterSlaveContainerVectorType;
    typedef typename BaseType::ConstraintIdsVectorType ConstraintIdsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor
     * @param rMainModelPart The reference to the modelpart which will be used for computations later on.
     * @param iParameters The settings parameters.
     */
    explicit ApplyChimeraProcessMonolithic(ModelPart &rMainModelPart, Parameters iParameters) : BaseType(rMainModelPart, iParameters)
    {
    }

    /// Destructor.
    virtual ~ApplyChimeraProcessMonolithic() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    std::string Info() const override
    {
        return "ApplyChimeraProcessMonolithic";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "ApplyChimeraProcessMonolithic" << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        KRATOS_INFO("ApplyChimeraProcessMonolithic") << std::endl;
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
     * @param rBinLocator The bin based locator formulated on the background. This is used to locate nodes of rBoundaryModelPart on background.
     */
    void ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorType &rBinLocator) override
    {
        auto& r_comm = BaseType::mrMainModelPart.GetCommunicator().GetDataCommunicator();
        //loop over nodes and find the triangle in which it falls, then do interpolation
        MasterSlaveContainerVectorType master_slave_container_vector;
        BaseType::ReserveMemoryForConstraintContainers(rBoundaryModelPart, master_slave_container_vector);
        std::vector<int> constraints_id_vector;
        int num_constraints_required = (TDim + 1) * (rBoundaryModelPart.Nodes().size());
        BaseType::CreateConstraintIds(constraints_id_vector, num_constraints_required);

        const int n_boundary_nodes = static_cast<int> (rBoundaryModelPart.Nodes().size());
        for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
        {
            ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
            BaseType::mNodeIdToConstraintIdsMap[i_boundary_node->Id()].reserve(150);
        }

        BaseType::FormulateConstraints(rBoundaryModelPart, rBinLocator, master_slave_container_vector, master_slave_container_vector);
        BuiltinTimer mpc_add_time;
        BaseType::AddConstraintsToModelpart(BaseType::mrMainModelPart, master_slave_container_vector);
        auto add_time = mpc_add_time.ElapsedSeconds();
        KRATOS_INFO_IF("Adding of MPCs from containers to modelpart took         : ", BaseType::mEchoLevel > 1)<< r_comm.Max(add_time, 0)<< " seconds"<< std::endl;
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
    ApplyChimeraProcessMonolithic &operator=(ApplyChimeraProcessMonolithic const &rOther);

    ///@}
}; // Class ApplyChimeraProcessMonolithic

} // namespace Kratos.

#endif //  KRATOS_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED defined 
