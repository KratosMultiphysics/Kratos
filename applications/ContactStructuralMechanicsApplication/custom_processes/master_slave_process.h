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

#if !defined(KRATOS_MASTER_SLAVE_PROCESS)
#define KRATOS_MASTER_SLAVE_PROCESS

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
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
 * @class MasterSlaveProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process assigns as master/slave the conditions
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) MasterSlaveProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MasterSlaveProcess
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveProcess);

    /// Index type definition
    typedef std::size_t IndexType;

    /// General type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MasterSlaveProcess( ModelPart& rThisModelPart)
        :mrThisModelPart(rThisModelPart)
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~MasterSlaveProcess() override
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
    std::string Info() const override
    {
        return "MasterSlaveProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MasterSlaveProcess";
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

    ModelPart& mrThisModelPart; /// The model part where to mcompute the master/slave flags

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
    MasterSlaveProcess& operator=(MasterSlaveProcess const& rOther) = delete;

    /// Copy constructor.
    //MasterSlaveProcess(MasterSlaveProcess const& rOther);


    ///@}

}; // Class MasterSlaveProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MasterSlaveProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MasterSlaveProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_MASTER_SLAVE_PROCESS defined */
