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

    static std::vector<IndexType> GetContainerIds(ModelPart& rModelPart);

    static void BooleanOperation(
        ModelPart& rModelPartA,
        ModelPart& rModelPartB,
        ModelPart& rDestination,
        BooleanOperators ThisOperator);

    ///@}

}; // Class SubModelPartEntitiesBooleanOperationUtility

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SUB_MODEL_PART_OPERATIONS_H_INCLUDED  defined
