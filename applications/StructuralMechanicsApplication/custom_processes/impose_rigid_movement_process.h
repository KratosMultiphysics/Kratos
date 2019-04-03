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

#if !defined(KRATOS_IMPOSE_RIGID_MOVEMENT_PROCESS)
#define KRATOS_IMPOSE_RIGID_MOVEMENT_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
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
    
/** 
 * @class ImposeRigidMovementProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This method assign linear kinematic constrains to a certain submodelpart
 * @details It assigns to the first node of the submodelpart by default or to a certain node if provided
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ImposeRigidMovementProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImposeRigidMovementProcess
    KRATOS_CLASS_POINTER_DEFINITION(ImposeRigidMovementProcess);
    
    /// General type definitions
    typedef Node<3>                                                      NodeType;

    /// General containers type definitions
    typedef ModelPart::MasterSlaveConstraintContainerType ConstraintContainerType;

    /// Component definition
    typedef VectorComponentAdaptor< array_1d< double, 3 > >         ComponentType;

    /// Definitions of the integers
    typedef std::size_t                                                 IndexType;
    typedef std::size_t                                                  SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rThisModelPart The model part to compute
     * @param ThisParameters The parameters of configuration
     */
    ImposeRigidMovementProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~ImposeRigidMovementProcess() override
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
     * @brief This function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override;
    
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
        return "ImposeRigidMovementProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ImposeRigidMovementProcess";
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
    
    ModelPart& mrThisModelPart; /// The model part to compute
    Parameters mThisParameters; /// The parameters (can be used for general pourposes)

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
    ImposeRigidMovementProcess& operator=(ImposeRigidMovementProcess const& rOther) = delete;

    /// Copy constructor.
    //ImposeRigidMovementProcess(ImposeRigidMovementProcess const& rOther);


    ///@}

}; // Class ImposeRigidMovementProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ImposeRigidMovementProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ImposeRigidMovementProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}
#endif /* KRATOS_IMPOSE_RIGID_MOVEMENT_PROCESS defined */
