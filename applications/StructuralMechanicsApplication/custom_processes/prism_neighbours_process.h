// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PRISM_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_PRISM_NEIGHBOURS_PROCESS_H_INCLUDED

// System includes

// External includes
#include <unordered_map>

// Project includes
#include "processes/process.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"

namespace Kratos
{

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

/**
 * @class PrismNeighboursProcess
 * @ingroup StructuralMechanicsApplication
 * @brief An algorithm that looks for neighbour nodes and elements in a mesh of prismatic elements
 * @details For that pourpose if builds an unordered map of the surrounding elements and nodes and performs different checks
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PrismNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PrismNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(PrismNeighboursProcess);

    // General geometry type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;

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

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Definition of the vector indexes considered
    typedef vector<IndexType> VectorIndexType;

    /// Definition of the hasher considered
    typedef VectorIndexHasher<VectorIndexType> VectorIndexHasherType;

    /// Definition of the key comparor considered
    typedef VectorIndexComparor<VectorIndexType> VectorIndexComparorType;

    /// Define the map considered for indexes
    typedef std::unordered_map<VectorIndexType, IndexType, VectorIndexHasherType, VectorIndexComparorType > HashMapVectorIntIntType;

    /// Define the HashMapVectorIntIntType iterator type
    typedef HashMapVectorIntIntType::iterator HashMapVectorIntIntIteratorType;

    /// Define the map considered for elemento pointers
    typedef std::unordered_map<VectorIndexType, Element::Pointer, VectorIndexHasherType, VectorIndexComparorType > HashMapVectorIntElementPointerType;

    /// Define the HashMapVectorIntElementPointerType iterator type
    typedef HashMapVectorIntElementPointerType::iterator HashMapVectorIntElementPointerIteratorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part where the search of neighbours is performed
     * @param ComputeOnNodes If true it will compute neighbours on nodes besides the elements. False by default to save memory
     */
    PrismNeighboursProcess(
        ModelPart& rModelPart,
        const bool ComputeOnNodes = false
        ) : mrModelPart(rModelPart),
            mComputeOnNodes(ComputeOnNodes)
    {
    }

    /// Destructor.
    virtual ~PrismNeighboursProcess() {}

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

    /**
     * @brief This function is designed for being called at the end of the computations right after reading the model and the groups
     */
    void ExecuteFinalize() override;

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
        return "PrismNeighboursProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PrismNeighboursProcess";
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
    const bool mComputeOnNodes; /// If true it will compute neighbours on nodes besides the elements. False by default to save memory

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method should be called in case that the current list of neighbour must be drop
     */
    void ClearNeighbours();

    /**
     * @brief This method add a unique weak pointer for any class
     * @param rPointerVector The vector containing the pointers of the defined class
     * @param Candidate The potential candidate  to add to the vector of pointers
     * @tparam TDataType The class type of the pointer
     * @todo Move this method to a common class (it is reused in the prism neighbour)
     */
    template< class TDataType >
    void  AddUniqueWeakPointer(
        WeakPointerVector< TDataType >& rPointerVector,
        const typename TDataType::WeakPointer Candidate
        )
    {
        typename WeakPointerVector< TDataType >::iterator beginit = rPointerVector.begin();
        typename WeakPointerVector< TDataType >::iterator endit   = rPointerVector.end();
        while ( beginit != endit && beginit->Id() != (Candidate.lock())->Id()) {
            beginit++;
        }
        if( beginit == endit ) {
            rPointerVector.push_back(Candidate);
        }

    }

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
    PrismNeighboursProcess& operator=(PrismNeighboursProcess const& rOther);

    ///@}

}; // Class PrismNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PrismNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PrismNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PRISM_NEIGHBOURS_PROCESS_H_INCLUDED  defined
