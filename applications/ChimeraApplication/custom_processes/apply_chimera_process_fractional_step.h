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

template <int TDim, class TSparseSpaceType, class TLocalSpaceType>
class KRATOS_API(CHIMERA_APPLICATION) ApplyChimeraProcessFractionalStep : public ApplyChimera<TDim, TSparseSpaceType, TLocalSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ApplyChimeraProcessFractionalStep
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessFractionalStep);
    typedef ApplyChimera<TDim, TSparseSpaceType, TLocalSpaceType> BaseType;
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
     * @param pBinLocator The bin based locator formulated on the background. This is used to locate nodes on rBoundaryModelPart.
     */
    void ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorPointerType &pBinLocator) override
    {
        //loop over nodes and find the triangle in which it falls, then do interpolation
        MasterSlaveContainerVectorType velocity_ms_container_vector;
        MasterSlaveContainerVectorType pressure_ms_container_vector;
#pragma omp parallel
        {
            const IndexType num_constraints_per_thread = (rBoundaryModelPart.NumberOfNodes() * 4)/omp_get_num_threads();
#pragma omp single
            {
                velocity_ms_container_vector.resize(omp_get_num_threads());
                for (auto &container : velocity_ms_container_vector)
                    container.reserve(num_constraints_per_thread);

                pressure_ms_container_vector.resize(omp_get_num_threads());
                for (auto &container : pressure_ms_container_vector)
                    container.reserve(num_constraints_per_thread);
            }
        }
        std::vector<int> constraints_id_vector;

        int num_constraints_required = (TDim + 1) * (rBoundaryModelPart.Nodes().size());
        BaseType::CreateConstraintIds(constraints_id_vector, num_constraints_required);

        const int max_results = 10000;
        const int n_boundary_nodes = static_cast<int> (rBoundaryModelPart.Nodes().size());
        IndexType counter = 0;
        IndexType removed_counter = 0;
        IndexType not_found_counter = 0;

        for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
        {
            ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
            Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());

            BaseType::mNodeIdToConstraintIdsMap[p_boundary_node->Id()].reserve(150);
        }

        BuiltinTimer loop_over_b_nodes;

#pragma omp parallel for shared(constraints_id_vector, velocity_ms_container_vector, pressure_ms_container_vector, pBinLocator) reduction(+                                                             \
                                                                                                                                          : not_found_counter) reduction(+                              \
                                                                                                                                                                         : removed_counter) reduction(+ \
                                                                                                                                                                                                      : counter)
        for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
        {

            Vector shape_fun_weights;
            typename PointLocatorType::ResultContainerType results(max_results);
            auto &velocity_ms_container = velocity_ms_container_vector[omp_get_thread_num()]; //TODO: change this to out of loop.
            auto &pressure_ms_container = pressure_ms_container_vector[omp_get_thread_num()]; //TODO: change this to out of loop.

            ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
            Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());
            ConstraintIdsVectorType constrainIds_for_the_node;
            IndexType start_constraint_id = i_bn * (TDim + 1) * (TDim + 1);
            bool node_coupled = false;
            if ((p_boundary_node)->IsDefined(VISITED))
                node_coupled = (p_boundary_node)->Is(VISITED);

            typename PointLocatorType::ResultIteratorType result_begin = results.begin();
            Element::Pointer p_element;
            bool is_found = false;
            is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), shape_fun_weights, p_element, result_begin, max_results);

            if (node_coupled && is_found)
            {
                constrainIds_for_the_node = BaseType::mNodeIdToConstraintIdsMap[p_boundary_node->Id()];
                for (auto const &constraint_id : constrainIds_for_the_node)
                {
#pragma omp critical
                    {
                        BaseType::mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
                        removed_counter++;
                    }
                }
                constrainIds_for_the_node.clear();
                p_boundary_node->Set(VISITED, false);
            }

            if (is_found == true)
            {
                Geometry<Node<3>> &r_geom = p_element->GetGeometry();
                int init_index = 0;
                BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_X, start_constraint_id + init_index, constraints_id_vector, velocity_ms_container);
                init_index += (TDim+1);
                BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_Y, start_constraint_id + init_index, constraints_id_vector, velocity_ms_container);
                init_index += (TDim+1);
                if (TDim == 3)
                {
                    BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_Z, start_constraint_id + init_index, constraints_id_vector, velocity_ms_container);
                    init_index += (TDim+1);
                }
                BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, PRESSURE, start_constraint_id + init_index, constraints_id_vector, pressure_ms_container);
                init_index += (TDim+1);
                counter += 1;
            }
            else
            {
                not_found_counter += 1;
            }
            p_boundary_node->Set(VISITED, true);
        } // end of loop over boundary nodes

        KRATOS_INFO_IF("Loop over boundary nodes took : ", BaseType::mEchoLevel > 1)<< loop_over_b_nodes.ElapsedSeconds()<< " seconds"<< std::endl;

        BuiltinTimer mpc_add_time;
        IndexType n_vel_constraints = 0;
        for (auto &container : velocity_ms_container_vector)
        {
            const int n_constraints = static_cast<int> (container.size());
            n_vel_constraints += n_constraints;
            #pragma omp parallel for
            for (int i_con = 0; i_con < n_constraints; ++i_con)
            {
                auto i_container = container.begin() + i_con;
                i_container->Set(FS_CHIMERA_PRE_CONSTRAINT, false);
                i_container->Set(FS_CHIMERA_VEL_CONSTRAINT, true);
                i_container->Set(ACTIVE);
            }
            // auto &vel_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_velocity_model_part");
            // vel_modelpart.AddMasterSlaveConstraints(container.begin(), container.end());
        }

        auto &vel_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_velocity_model_part");
        auto& vel_constraints = vel_modelpart.MasterSlaveConstraints();
        vel_constraints.reserve(n_vel_constraints);
        auto& vel_constraints_data = vel_constraints.GetContainer();
        for (auto &container : velocity_ms_container_vector)
        {
            vel_constraints_data.insert( vel_constraints_data.end(), container.ptr_begin(), container.ptr_end() );
        }
        vel_constraints.Sort();


        IndexType n_pre_constraints = 0;
        for (auto &container : pressure_ms_container_vector)
        {
            const int n_constraints = static_cast<int> (container.size());
            n_pre_constraints += n_constraints;
            #pragma omp parallel for
            for (int i_con = 0; i_con < n_constraints; ++i_con)
            {
                auto i_container = container.begin() + i_con;
                i_container->Set(FS_CHIMERA_PRE_CONSTRAINT, true);
                i_container->Set(FS_CHIMERA_VEL_CONSTRAINT, false);
                i_container->Set(ACTIVE);
            }
            // auto &pre_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_pressure_model_part");
            // pre_modelpart.AddMasterSlaveConstraints(container.begin(), container.end());
        }

        auto &pre_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_pressure_model_part");
        auto& pre_constraints = pre_modelpart.MasterSlaveConstraints();
        pre_constraints.reserve(n_pre_constraints);
        auto& pre_constraints_data = pre_constraints.GetContainer();
        for (auto &container : pressure_ms_container_vector)
        {
            pre_constraints_data.insert( pre_constraints_data.end(), container.ptr_begin(), container.ptr_end() );
        }
        pre_constraints.Sort();


        KRATOS_INFO_IF("Adding of MPCs from containers to modelpart took : ", BaseType::mEchoLevel > 1)<< mpc_add_time.ElapsedSeconds()<< " seconds"<< std::endl;

        KRATOS_INFO_IF("Number of boundary nodes in : ", BaseType::mEchoLevel > 1) << rBoundaryModelPart.Name() << " is coupled " << rBoundaryModelPart.NumberOfNodes() << std::endl;
        KRATOS_INFO_IF("Number of Boundary nodes found : ", BaseType::mEchoLevel > 1) << counter << ". Number of constraints : " << counter * 9 << std::endl;
        KRATOS_INFO_IF("Number of Boundary nodes not found  : ", BaseType::mEchoLevel > 1) << not_found_counter << std::endl;
        KRATOS_INFO_IF("Number of constraints removed  : ", BaseType::mEchoLevel > 1) << removed_counter << std::endl;
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
