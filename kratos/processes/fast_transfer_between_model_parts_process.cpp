//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "processes/fast_transfer_between_model_parts_process.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

FastTransferBetweenModelPartsProcess::FastTransferBetweenModelPartsProcess(
    ModelPart& rDestinationModelPart,
    ModelPart& rOriginModelPart,
    const EntityTransfered Entity,
    const Flags Flag,
    const bool ReplicateEntities
    ) : Process(),
        mrDestinationModelPart(rDestinationModelPart),
        mrOriginModelPart(rOriginModelPart),
        mEntity(Entity),
        mFlag(Flag)
{
    KRATOS_TRY

    // If the entities are replicated or transferred
    if (ReplicateEntities) {
        this->Set(MODIFIED, true);
    } else {
        this->Set(MODIFIED, false);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::operator()()
{
    Execute();
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::Execute()
{
    KRATOS_TRY;

    if (this->IsNot(MODIFIED)) {
        // In case of not flag defined we transfer all the elements
        if (mFlag == Flags()) {
            TransferWithoutFlags();
        } else {
            TransferWithFlags();
        }
    } else {
        // In case of not flag defined we transfer all the elements
        if (mFlag == Flags()) {
            ReplicateWithoutFlags();
        } else {
            ReplicateWithFlags();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::TransferWithoutFlags()
{
    const SizeType num_nodes = mrOriginModelPart.NumberOfNodes();

    if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS || mEntity == EntityTransfered::NODESANDGEOMETRIES))
        mrDestinationModelPart.AddNodes(mrOriginModelPart.NodesBegin(),mrOriginModelPart.NodesEnd());

    const SizeType num_elements = mrOriginModelPart.NumberOfElements();

    if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
        mrDestinationModelPart.AddElements(mrOriginModelPart.ElementsBegin(),mrOriginModelPart.ElementsEnd());

    const SizeType num_conditions = mrOriginModelPart.NumberOfConditions();

    if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
        mrDestinationModelPart.AddConditions(mrOriginModelPart.ConditionsBegin(),mrOriginModelPart.ConditionsEnd());

    const SizeType num_constraints = mrOriginModelPart.NumberOfMasterSlaveConstraints();

    if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS))
        mrDestinationModelPart.AddMasterSlaveConstraints(mrOriginModelPart.MasterSlaveConstraintsBegin(),mrOriginModelPart.MasterSlaveConstraintsEnd());

    const SizeType num_geometries = mrOriginModelPart.NumberOfGeometries();

    if (num_geometries != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::GEOMETRIES || mEntity == EntityTransfered::NODESANDGEOMETRIES))
        mrDestinationModelPart.AddGeometries(mrOriginModelPart.GeometriesBegin(),mrOriginModelPart.GeometriesEnd());
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::TransferWithFlags()
{
    // Auxiliar sizes
    const int num_nodes = static_cast<int>(mrOriginModelPart.NumberOfNodes());
    const int num_elements = static_cast<int>(mrOriginModelPart.NumberOfElements());
    const int num_conditions = static_cast<int>(mrOriginModelPart.NumberOfConditions());
    const int num_constraints = static_cast<int>(mrOriginModelPart.NumberOfMasterSlaveConstraints());
    const int num_geometries = static_cast<int>(mrOriginModelPart.NumberOfGeometries());

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        NodesArrayType nodes_buffer_vector;
        ElementsArrayType elements_buffer_vector;
        ConditionsArrayType conditions_buffer_vector;
        MasterSlaveConstraintArrayType constraints_buffer_vector;
        GeometriesMapType geometries_buffer_vector;

        if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS || mEntity == EntityTransfered::NODESANDGEOMETRIES)) {
            const auto it_node_begin = mrOriginModelPart.NodesBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_nodes; ++i) {
                auto it_node = it_node_begin + i;
                if (it_node->Is(mFlag)) {
                    nodes_buffer_vector.insert(nodes_buffer_vector.begin(), *(it_node.base()));
                }
            }
        }

        if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS)) {
            const auto it_elem_begin = mrOriginModelPart.ElementsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_elements; ++i) {
                auto it_elem = it_elem_begin + i;
                if (it_elem->Is(mFlag)) {
                    elements_buffer_vector.insert(elements_buffer_vector.begin(), *(it_elem.base()));
                }
            }
        }

        if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            const auto it_cond_begin = mrOriginModelPart.ConditionsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_conditions; ++i) {
                auto it_cond = it_cond_begin + i;
                if (it_cond->Is(mFlag)) {
                    conditions_buffer_vector.insert(conditions_buffer_vector.begin(), *(it_cond.base()));
                }
            }
        }

        if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS)) {
            const auto it_const_begin = mrOriginModelPart.MasterSlaveConstraintsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_constraints; ++i) {
                auto it_const = it_const_begin + i;
                if (it_const->Is(mFlag)) {
                    constraints_buffer_vector.insert(constraints_buffer_vector.begin(), *(it_const.base()));
                }
            }
        }

        if (num_geometries != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::GEOMETRIES || mEntity == EntityTransfered::NODESANDGEOMETRIES)) {
            const auto it_geom_begin = mrOriginModelPart.GeometriesBegin();
            #pragma omp for
            for(int i = 0; i < num_geometries; ++i) {
                auto it_geom = it_geom_begin;
                for (int j = 0; j < i; ++j) it_geom++;
                geometries_buffer_vector.insert(*(it_geom.base()));
            }
        }

        // We transfer
        #pragma omp critical
        {
            if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS || mEntity == EntityTransfered::NODESANDGEOMETRIES))
                mrDestinationModelPart.AddNodes(nodes_buffer_vector.begin(),nodes_buffer_vector.end());

            if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
                mrDestinationModelPart.AddElements(elements_buffer_vector.begin(),elements_buffer_vector.end());

            if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                mrDestinationModelPart.AddConditions(conditions_buffer_vector.begin(),conditions_buffer_vector.end());

            if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS))
                mrDestinationModelPart.AddMasterSlaveConstraints(constraints_buffer_vector.begin(),constraints_buffer_vector.end());

            if (num_geometries != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::GEOMETRIES || mEntity == EntityTransfered::NODESANDGEOMETRIES))
                mrDestinationModelPart.AddGeometries(geometries_buffer_vector.begin(),geometries_buffer_vector.end());
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::ReorderAllIds(ModelPart& rThisModelPart)
{
    auto& r_nodes_array = rThisModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    IndexPartition<std::size_t>(r_nodes_array.size()).for_each([&](std::size_t i) {
        (it_node_begin + i)->SetId(i + 1);
    });

    auto& r_element_array = rThisModelPart.Elements();
    const auto it_elem_begin = r_element_array.begin();
    IndexPartition<std::size_t>(r_element_array.size()).for_each([&](std::size_t i) {
        (it_elem_begin + i)->SetId(i + 1);
    });

    auto& r_condition_array = rThisModelPart.Conditions();
    const auto it_cond_begin = r_condition_array.begin();
    IndexPartition<std::size_t>(r_condition_array.size()).for_each([&](std::size_t i) {
        (it_cond_begin + i)->SetId(i + 1);
    });

    auto& r_constraint_array = rThisModelPart.MasterSlaveConstraints();
    const auto it_const_begin = r_constraint_array.begin();
    IndexPartition<std::size_t>(r_constraint_array.size()).for_each([&](std::size_t i) {
        (it_const_begin + i)->SetId(i + 1);
    });

    auto& r_geometries_array = rThisModelPart.Geometries();
    const auto it_geom_begin = r_geometries_array.begin();
    IndexPartition<std::size_t>(r_geometries_array.size()).for_each([&](std::size_t i) {
        auto it_geom = it_geom_begin;
        for (IndexType j = 0; j < i; ++j) it_geom++;
        it_geom->SetId(i + 1);
    });
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::ReplicateWithoutFlags()
{
    // Reordering ids (necessary for consistency)
    ModelPart& r_root_model_part = mrOriginModelPart.GetRootModelPart();
    ReorderAllIds(r_root_model_part);

    // Getting the auxiliar values
    const SizeType total_num_nodes = r_root_model_part.NumberOfNodes();
    const int num_nodes = static_cast<int>(mrOriginModelPart.NumberOfNodes());
    const SizeType total_num_elements = r_root_model_part.NumberOfElements();
    const int num_elements = static_cast<int>(mrOriginModelPart.NumberOfElements());
    const SizeType total_num_conditions = r_root_model_part.NumberOfConditions();
    const int num_conditions = static_cast<int>(mrOriginModelPart.NumberOfConditions());
    const SizeType total_num_constraints = r_root_model_part.NumberOfMasterSlaveConstraints();
    const int num_constraints = static_cast<int>(mrOriginModelPart.NumberOfMasterSlaveConstraints());
    const SizeType total_num_geometries = r_root_model_part.NumberOfGeometries();
    const int num_geometries = static_cast<int>(mrOriginModelPart.NumberOfGeometries());

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        NodesArrayType nodes_buffer_vector;
        ElementsArrayType elements_buffer_vector;
        ConditionsArrayType conditions_buffer_vector;
        MasterSlaveConstraintArrayType constraints_buffer_vector;
        GeometriesMapType geometries_buffer_vector;

        if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS || mEntity == EntityTransfered::NODESANDGEOMETRIES)) {
            const auto it_node_begin = mrOriginModelPart.NodesBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_nodes; ++i) {
                auto it_node = it_node_begin + i;
                NodeType::Pointer p_new_node = it_node->Clone();
                p_new_node->SetId(total_num_nodes + i + 1);
                nodes_buffer_vector.insert(nodes_buffer_vector.begin(), p_new_node);
            }
        }

        if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS)) {
            const auto it_elem_begin = mrOriginModelPart.ElementsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_elements; ++i) {
                auto it_elem = it_elem_begin + i;
                Element::Pointer p_new_elem = it_elem->Clone(total_num_elements + i + 1, it_elem->GetGeometry());
                elements_buffer_vector.insert(elements_buffer_vector.begin(), p_new_elem);
            }
        }

        if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            const auto it_cond_begin = mrOriginModelPart.ConditionsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_conditions; ++i) {
                auto it_cond = it_cond_begin + i;
                Condition::Pointer p_new_cond = it_cond->Clone(total_num_conditions + i + 1, it_cond->GetGeometry());
                conditions_buffer_vector.insert(conditions_buffer_vector.begin(), p_new_cond);
            }
        }

        if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS)) {
            const auto it_const_begin = mrOriginModelPart.MasterSlaveConstraintsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_constraints; ++i) {
                auto it_const = it_const_begin + i;
                MasterSlaveConstraint::Pointer p_new_const = it_const->Clone(total_num_constraints + i + 1);
                constraints_buffer_vector.insert(constraints_buffer_vector.begin(), p_new_const);
            }
        }

        if (num_geometries != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::GEOMETRIES || mEntity == EntityTransfered::NODESANDGEOMETRIES)) {
            const auto it_geom_begin = mrOriginModelPart.GeometriesBegin();
            #pragma omp for
            for(int i = 0; i < num_geometries; ++i) {
                auto it_geom = it_geom_begin;
                for (int j = 0; j < i; ++j) it_geom++;
                auto p_new_geom = it_geom->Create(total_num_geometries + i + 1, *it_geom);
                geometries_buffer_vector.insert(p_new_geom);
            }
        }

        // We add to the model part
        #pragma omp critical
        {
            if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS || mEntity == EntityTransfered::NODESANDGEOMETRIES))
                mrDestinationModelPart.AddNodes(nodes_buffer_vector.begin(),nodes_buffer_vector.end());

            if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
                mrDestinationModelPart.AddElements(elements_buffer_vector.begin(),elements_buffer_vector.end());

            if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                mrDestinationModelPart.AddConditions(conditions_buffer_vector.begin(),conditions_buffer_vector.end());

            if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS))
                mrDestinationModelPart.AddMasterSlaveConstraints(constraints_buffer_vector.begin(),constraints_buffer_vector.end());

            if (num_geometries != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::GEOMETRIES || mEntity == EntityTransfered::NODESANDGEOMETRIES))
                mrDestinationModelPart.AddGeometries(geometries_buffer_vector.begin(),geometries_buffer_vector.end());
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FastTransferBetweenModelPartsProcess::ReplicateWithFlags()
{
    // Reordering ids (necessary for consistency)
    ModelPart& r_root_model_part = mrOriginModelPart.GetRootModelPart();
    ReorderAllIds(r_root_model_part);

    // Getting the auxiliar values
    const SizeType total_num_nodes = r_root_model_part.NumberOfNodes();
    const int num_nodes = static_cast<int>(mrOriginModelPart.NumberOfNodes());
    const SizeType total_num_elements = r_root_model_part.NumberOfElements();
    const int num_elements = static_cast<int>(mrOriginModelPart.NumberOfElements());
    const SizeType total_num_conditions = r_root_model_part.NumberOfConditions();
    const int num_conditions = static_cast<int>(mrOriginModelPart.NumberOfConditions());
    const SizeType total_num_constraints = r_root_model_part.NumberOfMasterSlaveConstraints();
    const int num_constraints = static_cast<int>(mrOriginModelPart.NumberOfMasterSlaveConstraints());
    const SizeType total_num_geometries = r_root_model_part.NumberOfGeometries();
    const int num_geometries = static_cast<int>(mrOriginModelPart.NumberOfGeometries());

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        NodesArrayType nodes_buffer_vector;
        ElementsArrayType elements_buffer_vector;
        ConditionsArrayType conditions_buffer_vector;
        MasterSlaveConstraintArrayType constraints_buffer_vector;
        GeometriesMapType geometries_buffer_vector;

        if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS || mEntity == EntityTransfered::NODESANDGEOMETRIES)) {
            const auto it_node_begin = mrOriginModelPart.NodesBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_nodes; ++i) {
                auto it_node = it_node_begin + i;
                if (it_node->Is(mFlag)) {
                    NodeType::Pointer p_new_node = it_node->Clone();
                    p_new_node->SetId(total_num_nodes + i + 1);
                    (nodes_buffer_vector).insert(nodes_buffer_vector.begin(), p_new_node);
                }
            }
        }

        if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS)) {
            const auto it_elem_begin = mrOriginModelPart.ElementsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_elements; ++i) {
                auto it_elem = it_elem_begin + i;
                if (it_elem->Is(mFlag)) {
                    Element::Pointer p_new_elem = it_elem->Clone(total_num_elements + i + 1, it_elem->GetGeometry());
                    (elements_buffer_vector).insert(elements_buffer_vector.begin(), p_new_elem);
                }
            }
        }

        if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS)) {
            const auto it_cond_begin = mrOriginModelPart.ConditionsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_conditions; ++i) {
                auto it_cond = it_cond_begin + i;
                if (it_cond->Is(mFlag)) {
                    Condition::Pointer p_new_cond = it_cond->Clone(total_num_conditions + i + 1, it_cond->GetGeometry());
                    (conditions_buffer_vector).insert(conditions_buffer_vector.begin(), p_new_cond);
                }
            }
        }

        if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS)) {
            const auto it_const_begin = mrOriginModelPart.MasterSlaveConstraintsBegin();
            #pragma omp for schedule(guided, 512)
            for(int i = 0; i < num_constraints; ++i) {
                auto it_const = it_const_begin + i;
                if (it_const->Is(mFlag)) {
                    MasterSlaveConstraint::Pointer p_new_const = it_const->Clone(total_num_constraints + i + 1);
                    (constraints_buffer_vector).insert(constraints_buffer_vector.begin(), p_new_const);
                }
            }
        }

        if (num_geometries != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::GEOMETRIES || mEntity == EntityTransfered::NODESANDGEOMETRIES)) {
            const auto it_geom_begin = mrOriginModelPart.GeometriesBegin();
            #pragma omp for
            for(int i = 0; i < num_geometries; ++i) {
                auto it_geom = it_geom_begin;
                for (int j = 0; j < i; ++j) it_geom++;
                auto p_new_geom = it_geom->Create(total_num_geometries + i + 1, *it_geom);
                geometries_buffer_vector.insert(p_new_geom);
            }
        }

        // We add to the model part
        #pragma omp critical
        {
            if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS || mEntity == EntityTransfered::NODESANDCONDITIONS || mEntity == EntityTransfered::NODESANDGEOMETRIES))
                mrDestinationModelPart.AddNodes(nodes_buffer_vector.begin(),nodes_buffer_vector.end());

            if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
                mrDestinationModelPart.AddElements(elements_buffer_vector.begin(),elements_buffer_vector.end());

            if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS || mEntity == EntityTransfered::NODESANDCONDITIONS))
                mrDestinationModelPart.AddConditions(conditions_buffer_vector.begin(),conditions_buffer_vector.end());

            if (num_constraints != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONSTRAINTS || mEntity == EntityTransfered::NODESANDCONSTRAINTS))
                mrDestinationModelPart.AddMasterSlaveConstraints(constraints_buffer_vector.begin(),constraints_buffer_vector.end());

            if (num_geometries != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::GEOMETRIES || mEntity == EntityTransfered::NODESANDGEOMETRIES))
                mrDestinationModelPart.AddGeometries(geometries_buffer_vector.begin(),geometries_buffer_vector.end());
        }
    }
}

}
