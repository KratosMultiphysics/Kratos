// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ASSIGN_NODAL_ELEMENTS_TO_NODES_PROCESS)
#define KRATOS_ASSIGN_NODAL_ELEMENTS_TO_NODES_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
    /// The size definition
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
/** 
 * @class AssignNodalElementsToNodes
 * @ingroup StructuralMechanicsApplication
 * @brief This
 * @details It
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AssignNodalElementsToNodes
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignNodalElementsToNodes
    KRATOS_CLASS_POINTER_DEFINITION(AssignNodalElementsToNodes);
    
    /// The index definition
    typedef std::size_t                                     IndexType;

    /// Geometric type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;

    /// The definition of the containers
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    typedef ModelPart::ElementsContainerType        ElementsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rThisModelPart The model part to compute
     * @param ThisParameters The parameters of configuration
     */
    AssignNodalElementsToNodes(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~AssignNodalElementsToNodes() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
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
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitialize() override;
    
    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

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
    std::string Info() const override
    {
        return "AssignNodalElementsToNodes";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignNodalElementsToNodes";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
    
    ModelPart& mrThisModelPart;   /// The main model part
    Parameters mThisParameters;   /// The parameters (can be used for general pourposes)
    bool mConstantValues = true;  /// If the values to assign are constant or functions

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief After we have transfer the information from the previous modelpart we initilize the elements
     * @param rModelPart The model part contining the elements
     */

    void InitializeElements(ModelPart& rModelPart);

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
    AssignNodalElementsToNodes& operator=(AssignNodalElementsToNodes const& rOther) = delete;

    /// Copy constructor.
    //AssignNodalElementsToNodes(AssignNodalElementsToNodes const& rOther);


    ///@}

}; // Class AssignNodalElementsToNodes

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   AssignNodalElementsToNodes& rThis);
// 
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const AssignNodalElementsToNodes& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
// 
//     return rOStream;
// }

}
#endif /* KRATOS_ASSIGN_NODAL_ELEMENTS_TO_NODES_PROCESS defined */
