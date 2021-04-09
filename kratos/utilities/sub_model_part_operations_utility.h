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
 * @class SubModelPartOperationsUtility
 * @brief Wrapper of set operations: union, intersection and difference
 * @author Miguel Maso Sotomayor
 * @ingroup KratosCore
*/
class KRATOS_API(KRATOS_CORE) SubModelPartOperationsUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SubModelPartOperationsUtility
    KRATOS_CLASS_POINTER_DEFINITION(SubModelPartOperationsUtility);

    typedef std::size_t IndexType;

    ///@}
    ///@name  Enum's
    ///@{

    enum class SetOperators {Union, Intersection, Difference};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    /// Destructor.

    ///@}
    ///@name Operations
    ///@{

    template<class TEntityType, class TContainerType>
    void Union(ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination)
    {
        SetOperation<TEntityType, TContainerType>(
            rModelPartA, rModelPartB, rDestination, SetOperators::Union);
    }

    template<class TEntityType, class TContainerType>
    void Intersection(ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination)
    {
        SetOperation<TEntityType, TContainerType>(
            rModelPartA, rModelPartB, rDestination, SetOperators::Intersection);
    }

    template<class TEntityType, class TContainerType>
    void Difference(ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination)
    {
        SetOperation<TEntityType, TContainerType>(
            rModelPartA, rModelPartB, rDestination, SetOperators::Difference);
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    void AddEntities(const std::vector<IndexType>& rIds, ModelPart& rModelPart);

    template<class TContainerType>
    TContainerType& GetContainer(ModelPart& rModelPart);

    template<class TContainerType>
    std::vector<IndexType> GetContainerIds(ModelPart& rModelPart)
    {
        TContainerType& r_container = GetContainer<TContainerType>(rModelPart);
        std::vector<IndexType> ids_vector;
        for (auto it = r_container.begin(); it != r_container.end(); ++it)
        {
            ids_vector.push_back(it->Id());
        }
        return ids_vector;
    }

    template<class TEntityType, class TContainerType>
    void SetOperation(
        ModelPart& rModelPartA,
        ModelPart& rModelPartB,
        ModelPart& rDestination,
        SetOperators ThisOperator)
    {
        std::vector<IndexType> ids_a = GetContainerIds<TContainerType>(rModelPartA);
        std::vector<IndexType> ids_b = GetContainerIds<TContainerType>(rModelPartB);
        std::vector<IndexType> ids_destination;
        std::sort(ids_a.begin(), ids_a.end());
        std::sort(ids_b.begin(), ids_b.end());

        if (ThisOperator == SetOperators::Union) {
            std::set_union(
                ids_a.begin(), ids_a.end(),
                ids_b.begin(), ids_b.end(),
                std::back_inserter(ids_destination));
        }
        else if (ThisOperator == SetOperators::Intersection) {
            std::set_intersection(
                ids_a.begin(), ids_a.end(),
                ids_b.begin(), ids_b.end(),
                std::back_inserter(ids_destination));
        }
        else if (ThisOperator == SetOperators::Difference) {
            std::set_difference(
                ids_a.begin(), ids_a.end(),
                ids_b.begin(), ids_b.end(),
                std::back_inserter(ids_destination));
        }

        EntitiesEraseProcess<TEntityType>(rDestination, EntitiesEraseProcessFlags::ERASE_ALL_ENTITIES).Execute();
        AddEntities<TEntityType>(ids_destination, rDestination);
    }

    ///@}

}; // Class SubModelPartOperationsUtility

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SUB_MODEL_PART_OPERATIONS_H_INCLUDED  defined
