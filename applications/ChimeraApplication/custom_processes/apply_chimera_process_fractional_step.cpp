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

template <int TDim, class TDistanceCalculatorType>
void ApplyChimeraProcessFractionalStep<TDim, TDistanceCalculatorType>::ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorPointerType &pBinLocator)
{
    //loop over nodes and find the triangle in which it falls, then do interpolation
    MasterSlaveContainerVectorType velocity_ms_container_vector;
    MasterSlaveContainerVectorType pressure_ms_container_vector;
#pragma omp parallel
    {
#pragma omp single
        {
            velocity_ms_container_vector.resize(omp_get_num_threads());
            for (auto &container : velocity_ms_container_vector)
                container.reserve(1000);

            pressure_ms_container_vector.resize(omp_get_num_threads());
            for (auto &container : pressure_ms_container_vector)
                container.reserve(1000);
        }
    }
    std::vector<int> constraints_id_vector;

    int num_constraints_required = (TDim + 1) * (rBoundaryModelPart.Nodes().size());
    BaseType::CreateConstraintIds(constraints_id_vector, num_constraints_required);

    const int max_results = 10000;
    const unsigned int n_boundary_nodes = rBoundaryModelPart.Nodes().size();
    std::size_t counter = 0;
    std::size_t removed_counter = 0;
    std::size_t not_found_counter = 0;

    for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
    {
        ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
        Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());

        BaseType::mNodeIdToConstraintIdsMap[p_boundary_node->Id()].reserve(150);
    }

#pragma omp parallel for shared(constraints_id_vector, velocity_ms_container_vector, pressure_ms_container_vector, pBinLocator) reduction(+:not_found_counter) reduction(+:removed_counter) reduction(+ \
                                                                                                                                           : counter)
    for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
    {

        Vector shape_fun_weights;
        typename PointLocatorType::ResultContainerType results(max_results);
        auto &velocity_ms_container = velocity_ms_container_vector[omp_get_thread_num()]; //TODO: change this to out of loop.
        auto &pressure_ms_container = pressure_ms_container_vector[omp_get_thread_num()]; //TODO: change this to out of loop.

        ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
        Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());
        ConstraintIdsVectorType constrainIds_for_the_node;
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
            constrainIds_for_the_node = BaseType::mNodeIdToConstraintIdsMap[p_boundary_node->Id()];
            for (auto const &constraint_id : constrainIds_for_the_node)
            {
                BaseType::mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
                removed_counter++;
            }
            p_boundary_node->Set(VISITED, false);
        }

        if (is_found == true)
        {
            Geometry<Node<3>> &r_geom = p_element->GetGeometry();
            int init_index = 0;
            BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_X, start_constraint_id+init_index, constraints_id_vector, velocity_ms_container);
            init_index+=3;
            BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_Y, start_constraint_id+init_index, constraints_id_vector, velocity_ms_container);
            init_index+=3;
            if(TDim == 3){
                BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_Z, start_constraint_id+init_index, constraints_id_vector, velocity_ms_container);
                init_index+=3;
            }
            BaseType::ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, PRESSURE, start_constraint_id+init_index, constraints_id_vector, pressure_ms_container);
            init_index+=3;
            counter += 1;
        }else{
            not_found_counter+=1;
        }
        p_boundary_node->Set(VISITED, true);
    } // end of loop over boundary nodes

    for (auto &container : velocity_ms_container_vector){
        for(auto& constraint : container){
            // TODO: Set the FS_CHIMERA_VEL_CONSTRAINT variable to true
            constraint.Set(FS_CHIMERA_PRE_CONSTRAINT, false);
            constraint.Set(FS_CHIMERA_VEL_CONSTRAINT, true);
            constraint.Set(ACTIVE);
        }
        auto& vel_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_velocity_model_part");
        vel_modelpart.AddMasterSlaveConstraints(container.begin(), container.end());
    }

    for (auto &container : pressure_ms_container_vector){
        for(auto& constraint : container){
            // TODO: Set the FS_CHIMERA_PRE_CONSTRAINT variable to true
            constraint.Set(FS_CHIMERA_PRE_CONSTRAINT, true);
            constraint.Set(FS_CHIMERA_VEL_CONSTRAINT, false);
            constraint.Set(ACTIVE);
        }
        auto& pre_modelpart = BaseType::mrMainModelPart.GetSubModelPart("fs_pressure_model_part");
        pre_modelpart.AddMasterSlaveConstraints(container.begin(), container.end());
    }

    KRATOS_INFO("Number of boundary nodes in : ") << rBoundaryModelPart.Name() << " is coupled " << rBoundaryModelPart.NumberOfNodes() << std::endl;
    KRATOS_INFO("Number of Boundary nodes found : ") << counter<<". Number of constraints : "<<counter*9<< std::endl;
    KRATOS_INFO("Number of Boundary nodes not found  : ") << not_found_counter << std::endl;
}


typedef CustomCalculateSignedDistanceProcess<2> DistanceCalculator2DType;
typedef CustomCalculateSignedDistanceProcess<3> DistanceCalculator3DType;

template class ApplyChimeraProcessFractionalStep<2, DistanceCalculator2DType>;
template class ApplyChimeraProcessFractionalStep<3, DistanceCalculator3DType>;

} // namespace Kratos.
