//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_SUB_MODEL_PART_OPERATIONS_H_INCLUDED
#define KRATOS_SUB_MODEL_PART_OPERATIONS_H_INCLUDED


// External includes


// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "processes/entity_erase_process.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SubModelPartEntitiesBooleanOperationUtility
 * @brief Wrapper of boolean operations: union, intersection and difference
 * @author Miguel Maso Sotomayor
 * @ingroup KratosCore
*/
template<class TEntityType, class TContainerType>
class KRATOS_API(KRATOS_CORE) SubModelPartEntitiesBooleanOperationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SubModelPartEntitiesBooleanOperationUtility
    KRATOS_CLASS_POINTER_DEFINITION(SubModelPartEntitiesBooleanOperationUtility);

    typedef std::size_t IndexType;

    ///@}
    ///@name  Enum's
    ///@{

    enum class BooleanOperators {Union, Intersection, Difference};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    /// Destructor.

    ///@}
    ///@name Operations
    ///@{

    static void Union(ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination)
    {
        BooleanOperation(rModelPartA, rModelPartB, rDestination, BooleanOperators::Union);
    }

    static void Intersection(ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination)
    {
        BooleanOperation(rModelPartA, rModelPartB, rDestination, BooleanOperators::Intersection);
    }

    static void Difference(ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination)
    {
        BooleanOperation(rModelPartA, rModelPartB, rDestination, BooleanOperators::Difference);
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    static void AddEntities(const std::vector<IndexType>& rIds, ModelPart& rModelPart);

    static TContainerType& GetContainer(ModelPart& rModelPart);

    static std::vector<IndexType> GetContainerIds(ModelPart& rModelPart)
    {
        const TContainerType& r_container = GetContainer(rModelPart);
        std::vector<IndexType> ids_vector(r_container.size());
        IndexPartition<std::size_t>(r_container.size()).for_each([&](std::size_t i){
            ids_vector[i] = (r_container.begin()+i)->Id();
        });
        return ids_vector;
    }

    static void BooleanOperation(
        ModelPart& rModelPartA,
        ModelPart& rModelPartB,
        ModelPart& rDestination,
        BooleanOperators ThisOperator)
    {
        KRATOS_ERROR_IF(!rDestination.IsSubModelPart()) << "The destination model part must be a sub model part." << std::endl;
        const ModelPart& r_root_a = rModelPartA.GetRootModelPart();
        const ModelPart& r_root_b = rModelPartB.GetRootModelPart();
        const ModelPart& r_root_d = rDestination.GetRootModelPart();
        KRATOS_ERROR_IF(&r_root_a != &r_root_b) << "The first and second model parts must belong to the same root model part." << std::endl;
        KRATOS_ERROR_IF(&r_root_a != &r_root_d) << "The destination model part must belong to the same root model part than the first and the second." << std::endl;
        std::vector<IndexType> ids_a = GetContainerIds(rModelPartA);
        std::vector<IndexType> ids_b = GetContainerIds(rModelPartB);
        std::vector<IndexType> ids_destination;
        std::sort(ids_a.begin(), ids_a.end());
        std::sort(ids_b.begin(), ids_b.end());

        if (ThisOperator == BooleanOperators::Union) {
            std::set_union(
                ids_a.begin(), ids_a.end(),
                ids_b.begin(), ids_b.end(),
                std::back_inserter(ids_destination));
        } else if (ThisOperator == BooleanOperators::Intersection) {
            std::set_intersection(
                ids_a.begin(), ids_a.end(),
                ids_b.begin(), ids_b.end(),
                std::back_inserter(ids_destination));
        } else if (ThisOperator == BooleanOperators::Difference) {
            std::set_difference(
                ids_a.begin(), ids_a.end(),
                ids_b.begin(), ids_b.end(),
                std::back_inserter(ids_destination));
        }

        EntitiesEraseProcess<TEntityType>(rDestination, EntitiesEraseProcessFlags::ERASE_ALL_ENTITIES).Execute();
        AddEntities(ids_destination, rDestination);
    }

    ///@}

}; // Class SubModelPartEntitiesBooleanOperationUtility

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SUB_MODEL_PART_OPERATIONS_H_INCLUDED  defined
