//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

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
    typedef Variable<array_1d<double, 3>> Array3VariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    InterfaceVectorContainer(ModelPart& rModelPart);

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

    void UpdateSystemVectorFromModelPart(const Array3VariableType& rVariable,
                                         const Kratos::Flags& rMappingOptions);

    void UpdateModelPartFromSystemVector(const DoubleVariableType& rVariable,
                                         const Kratos::Flags& rMappingOptions);

    void UpdateModelPartFromSystemVector(const ComponentVariableType& rVariable,
                                         const Kratos::Flags& rMappingOptions);

    void UpdateModelPartFromSystemVector(const Array3VariableType& rVariable,
                                         const Kratos::Flags& rMappingOptions);


    // {
    //     // Here we construct a function pointer to not have the if all the time inside the loop
    //     const auto fill_fct = MapperUtilities::GetFillFunction<TVarType>(rMappingOptions);

    //     const int num_local_nodes = mrModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    //     const auto nodes_begin = mrModelPart.GetCommunicator().LocalMesh().NodesBegin();

    //     #pragma omp parallel for
    //     for (int i=0; i<num_local_nodes; i++) {
    //         fill_fct(*(nodes_begin + i), rVariable, *mpInterfaceVector[i]);
    //     }
    // }

    // {
    //     const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    //     // Here we construct a function pointer to not have the if all the time inside the loop
    //     const auto update_fct = std::bind(MapperUtilities::GetUpdateFunction<TVarType>(rMappingOptions),
    //                                         std::placeholders::_1,
    //                                         std::placeholders::_2,
    //                                         std::placeholders::_3,
    //                                         factor);
    //     const int num_local_nodes = mrModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    //     const auto nodes_begin = mrModelPart.GetCommunicator().LocalMesh().NodesBegin();

    //     #pragma omp parallel for
    //     for (int i=0; i<num_local_nodes; i++) {
    //         update_fct(*(nodes_begin + i), rVariable, *mpInterfaceVector[i]);
    //     }
    // }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

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

    ModelPart& mrModelPart
    TSystemVectorUniquePointerType mpInterfaceVector = nullptr;

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
    InterfaceVectorContainer& operator=(InterfaceVectorContainer const& rOther);

    /// Copy constructor.
    InterfaceVectorContainer(InterfaceVectorContainer const& rOther);


    ///@}

    }; // Class InterfaceVectorContainer

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                InterfaceVectorContainer& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const InterfaceVectorContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_VECTOR_CONTAINER_H_INCLUDED  defined
