//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "mapper_utilities.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
template<class TSparseSpace, class TDenseSpace>
class KRATOS_API(MAPPING_APPLICATION) InterfaceVectorContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceVectorContainer
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceVectorContainer);

    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef Kratos::unique_ptr<TSystemVectorType> TSystemVectorUniquePointerType;

    enum class InterfaceEntityType {
        NODES,
        ELEMENTS,
        CONDITIONS,
        GEOMETRIES
    };

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructor.
    explicit InterfaceVectorContainer(ModelPart& rModelPart)
        : InterfaceVectorContainer(rModelPart, InterfaceEntityType::NODES)
    {}

    /// Constructor with entity type
    InterfaceVectorContainer(ModelPart& rModelPart, InterfaceEntityType EntityType)
        : mrModelPart(rModelPart)
        , mInterfaceEntityType(EntityType)
    {}

    /// Destructor.
    virtual ~InterfaceVectorContainer() = default;


    ///@}
    ///@name Operations
    ///@{

    void UpdateSystemVectorFromModelPart(const Variable<double>& rVariable,
                                         const Kratos::Flags& rMappingOptions);

    void UpdateModelPartFromSystemVector(const Variable<double>& rVariable,
                                         const Kratos::Flags& rMappingOptions);

    ///@}
    ///@name Access
    ///@{

    TSystemVectorType& GetVector()
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpInterfaceVector)
            << "The Interface-Vector was not initialized" << std::endl;
        return *mpInterfaceVector;
    }

    const TSystemVectorType& GetVector() const
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpInterfaceVector)
            << "The Interface-Vector was not initialized" << std::endl;
        return *mpInterfaceVector;
    }

    TSystemVectorUniquePointerType& pGetVector() { return mpInterfaceVector; }
    const TSystemVectorUniquePointerType& pGetVector() const { return mpInterfaceVector; }

    ModelPart& GetModelPart() { return mrModelPart; }
    const ModelPart& GetModelPart() const { return mrModelPart; }

    ///@}

private:
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;
    InterfaceEntityType mInterfaceEntityType = InterfaceEntityType::NODES;
    TSystemVectorUniquePointerType mpInterfaceVector = nullptr;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    InterfaceVectorContainer& operator=(InterfaceVectorContainer const& rOther);

    /// Copy constructor.
    InterfaceVectorContainer(InterfaceVectorContainer const& rOther);

    ///@}

}; // Class InterfaceVectorContainer

///@}

///@} addtogroup block

}  // namespace Kratos.
