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

template <int TDim, class TSparseSpaceType, class TLocalSpaceType>
class KRATOS_API(CHIMERA_APPLICATION) ApplyChimeraProcessMonolithic : public ApplyChimera<TDim, TSparseSpaceType, TLocalSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessMonolithic);
    typedef ApplyChimera<TDim, TSparseSpaceType, TLocalSpaceType> BaseType;
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
     * @param pBinLocator The bin based locator formulated on the background. This is used to locate nodes of rBoundaryModelPart on background.
     */
    void ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorPointerType &pBinLocator) override
    {
        //loop over nodes and find the triangle in which it falls, then do interpolation
        MasterSlaveContainerVectorType master_slave_container_vector;
#pragma omp parallel
        {
#pragma omp single
            {
                master_slave_container_vector.resize(omp_get_num_threads());
                for (auto &container : master_slave_container_vector)
                    container.reserve(1000);
            }
        }
        std::vector<int> constraints_id_vector;

        int num_constraints_required = (TDim + 1) * (rBoundaryModelPart.Nodes().size());
        BaseType::CreateConstraintIds(constraints_id_vector, num_constraints_required);

        const int max_results = 10000;
        const unsigned int n_boundary_nodes = rBoundaryModelPart.Nodes().size();
        IndexType counter = 0;
        IndexType removed_counter = 0;
        IndexType not_found_counter = 0;

        for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
        {
            ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
            BaseType::mNodeIdToConstraintIdsMap[i_boundary_node->Id()].reserve(150);
        }

#pragma omp parallel for shared(constraints_id_vector, master_slave_container_vector, pBinLocator) reduction(+                                                             \
                                                                                                             : not_found_counter) reduction(+                              \
                                                                                                                                            : removed_counter) reduction(+ \
                                                                                                                                                                         : counter)
        for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
        {
            Vector shape_fun_weights;
            typename PointLocatorType::ResultContainerType results(max_results);
            auto &ms_container = master_slave_container_vector[omp_get_thread_num()];

            ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
            Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());
            unsigned int start_constraint_id = i_bn * (TDim + 1) * (TDim + 1);
            bool node_coupled = false;
            if ((p_boundary_node)->IsDefined(VISITED))
                node_coupled = (p_boundary_node)->Is(VISITED);

            typename PointLocatorType::ResultIteratorType result_begin = results.begin();
            Element::Pointer p_element;
            bool is_found = false;
            is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), shape_fun_weights, p_element, result_begin, max_results);

            if (node_coupled && is_found)
            {
                auto &constrainIds_for_the_node = BaseType::mNodeIdToConstraintIdsMap[p_boundary_node->Id()];
                for (auto const &constraint_id : constrainIds_for_the_node)
                {
#pragma omp critical
                    {
                        BaseType::mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
                    }
                }
                constrainIds_for_the_node.clear();
                p_boundary_node->Set(VISITED, false);
                removed_counter++;
            }

            if (is_found)
            {
                Geometry<Node<3>> &r_geom = p_element->GetGeometry();
                int init_index = 0;
                BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_X, start_constraint_id + init_index, constraints_id_vector, ms_container);
                init_index += 3;
                BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_Y, start_constraint_id + init_index, constraints_id_vector, ms_container);
                init_index += 3;
                if (TDim == 3)
                {
                    BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_Z, start_constraint_id + init_index, constraints_id_vector, ms_container);
                    init_index += 3;
                }
                BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, PRESSURE, start_constraint_id + init_index, constraints_id_vector, ms_container);
                init_index += 3;
                counter += 1;
            }
            else
            {
                not_found_counter += 1;
            }
            p_boundary_node->Set(VISITED, true);
        } // end of loop over boundary nodes

        for (auto &container : master_slave_container_vector)
            BaseType::mrMainModelPart.AddMasterSlaveConstraints(container.begin(), container.end());

        KRATOS_INFO_IF("Number of boundary nodes in : ", BaseType::mEchoLevel > 1) << rBoundaryModelPart.Name() << " is coupled " << rBoundaryModelPart.NumberOfNodes() << std::endl;
        KRATOS_INFO_IF("Number of Boundary nodes found : ", BaseType::mEchoLevel > 1) << counter << ". Number of constraints : " << counter * 9 << std::endl;
        KRATOS_INFO_IF("Number of Boundary nodes not found  : ", BaseType::mEchoLevel > 1) << not_found_counter << std::endl;
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
