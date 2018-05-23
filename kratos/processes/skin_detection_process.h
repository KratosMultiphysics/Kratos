//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SKIN_DETECTION_PROCESS_H_INCLUDED )
#define  KRATOS_SKIN_DETECTION_PROCESS_H_INCLUDED

// System includes

// External includes
#include <unordered_set>
#include <unordered_map>

// Project includes
#include "processes/process.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // General geometry type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SkinDetectionProcess
 * @ingroup KratosCore
 * @brief An algorithm that looks for neighbour elements in a mesh and creates a submodelpart containing the skin of the disconnected elements (interface elements)
 * @details For that pourpose if builds an unordered map of the surrounding elements and nodes and performs different checks.
 * @tparam TDim The dimension where the problem is computed
 * @author Vicente Mataix Ferrandiz
*/
template<SizeType TDim>
class KRATOS_API(KRATOS_CORE) SkinDetectionProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SkinDetectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(SkinDetectionProcess);

    // Containers definition
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    typedef ModelPart::ElementsContainerType        ElementsArrayType;

    // Containers iterators definition
    typedef NodesArrayType::iterator                NodesIterarorType;
    typedef ConditionsArrayType::iterator      ConditionsIteratorType;
    typedef ElementsArrayType::iterator          ElementsIteratorType;

    // Weak pointers vectors types
    typedef WeakPointerVector<NodeType> NodePointerVector;
    typedef WeakPointerVector<Element> ElementPointerVector;

    /// Definition of the vector indexes considered
    typedef std::vector<IndexType> VectorIndexType;

    /// Definition of the hasher considered
    typedef VectorIndexHasher<VectorIndexType> VectorIndexHasherType;

    /// Definition of the key comparor considered
    typedef VectorIndexComparor<VectorIndexType> VectorIndexComparorType;

    /// Define the set considered for elemento pointers
    typedef std::unordered_set<VectorIndexType, VectorIndexHasherType, VectorIndexComparorType > HashSetVectorIntType;

    /// Define the HashMapVectorIntElementPointerType iterator type
    typedef HashSetVectorIntType::iterator HashSetVectorIntTypeIteratorType;

    /// Define the map considered for elemento pointers
    typedef std::unordered_map<VectorIndexType, VectorIndexType, VectorIndexHasherType, VectorIndexComparorType > HashMapVectorIntType;

    /// Define the HashMapVectorIntElementPointerType iterator type
    typedef HashMapVectorIntType::iterator HashMapVectorIntTypeIteratorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part where the search of neighbours is performed
     * @param ThisParameters The parameters of configuration
     */
    SkinDetectionProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    virtual ~SkinDetectionProcess() {}

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method executes the algorithm that looks for neighbour nodes and elements in a  mesh of prismatic elements
     */
    void Execute() override;

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
    virtual std::string Info() const override
    {
        return "SkinDetectionProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SkinDetectionProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
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

    ModelPart& mrModelPart;     /// The main model part
    Parameters mThisParameters; /// The parameters (can be used for general pourposes)

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method computes the potential size for neighbours
     * @param itElem The element iterator where to check the size and so on
     * @return The reserve size for the neighbour elements vector
     * @todo Check that EdgesNumber() and FacesNumber() are properly implemeted on the geometries of interest
     */
    SizeType ComputePotentialNeighboursSize(ElementsIteratorType itElem);

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
    SkinDetectionProcess& operator=(SkinDetectionProcess const& rOther);

    ///@}

}; // Class SkinDetectionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<SizeType TDim>
inline std::istream& operator >> (std::istream& rIStream,
                                  SkinDetectionProcess<TDim>& rThis);

/// output stream function
template<SizeType TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SkinDetectionProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SKIN_DETECTION_PROCESS_H_INCLUDED  defined
