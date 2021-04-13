//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_VEC_CONTAINER_H_INCLUDED )
#define  KRATOS_INTERFACE_VEC_CONTAINER_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
template<class TVector>
class KRATOS_API(MAPPING_APPLICATION) InterfaceVecContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceVecContainer
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceVecContainer);

    using VectorType = TVector;

    using VectorUniquePointerType = Kratos::unique_ptr<VectorType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit InterfaceVecContainer(ModelPart& rModelPart) : mrModelPart(rModelPart) {}

    /// Destructor.
    virtual ~InterfaceVecContainer() = default;

    /// Copy constructor.
    InterfaceVecContainer(InterfaceVecContainer const& rOther) = delete;

    /// Assignment operator.
    InterfaceVecContainer& operator=(InterfaceVecContainer const& rOther) = delete;

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

    VectorType& GetVector()
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpInterfaceVector) << "The Interface-Vector was not initialized" << std::endl;
        return *mpInterfaceVector;
    }

    const VectorType& GetVector() const
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpInterfaceVector) << "The Interface-Vector was not initialized" << std::endl;
        return *mpInterfaceVector;
    }

    VectorUniquePointerType& pGetVector() { return mpInterfaceVector; }
    const VectorUniquePointerType& pGetVector() const { return mpInterfaceVector; }

    ModelPart& GetModelPart() { return mrModelPart; }
    const ModelPart& GetModelPart() const { return mrModelPart; }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    VectorUniquePointerType mpInterfaceVector = nullptr;

}; // Class InterfaceVecContainer

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_VEC_CONTAINER_H_INCLUDED  defined
