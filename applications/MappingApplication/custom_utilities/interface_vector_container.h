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

#if !defined(KRATOS_INTERFACE_VECTOR_CONTAINER_H_INCLUDED )
#define  KRATOS_INTERFACE_VECTOR_CONTAINER_H_INCLUDED

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
class InterfaceVectorContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceVectorContainer
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceVectorContainer);

    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef Kratos::unique_ptr<TSystemVectorType> TSystemVectorUniquePointerType;

    typedef Variable<double> DoubleVariableType;
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    InterfaceVectorContainer(ModelPart& rModelPart) : mrModelPart(rModelPart) {}

    /// Destructor.
    virtual ~InterfaceVectorContainer() = default;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void UpdateSystemVectorFromModelPart(const DoubleVariableType& rVariable,
                                         const Kratos::Flags& rMappingOptions);

    void UpdateSystemVectorFromModelPart(const ComponentVariableType& rVariable,
                                         const Kratos::Flags& rMappingOptions);

    void UpdateModelPartFromSystemVector(const DoubleVariableType& rVariable,
                                         const Kratos::Flags& rMappingOptions);

    void UpdateModelPartFromSystemVector(const ComponentVariableType& rVariable,
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

    TSystemVectorUniquePointerType& pGetVector()
    {
        return mpInterfaceVector;
    }

    ModelPart& GetModelPart()
    {
        return mrModelPart;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    TSystemVectorUniquePointerType mpInterfaceVector = nullptr;

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


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

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_VECTOR_CONTAINER_H_INCLUDED  defined
