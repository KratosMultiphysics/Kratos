//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_MASTER_SLAVE_PROCESS_H_INCLUDED )
#define  KRATOS_MASTER_SLAVE_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"

#include "processes/process.h"

#include "includes/dof.h"

#include "constraints/linear_master_slave_constraint.h"

#include "geometries/coupling_geometry.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class MasterSlaveProcess
 * @ingroup IgaApplication
 * @brief This class outputs the location of the quadrature points within the local space of the containing geometry. */
class KRATOS_API(IGA_APPLICATION) MasterSlaveProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MasterSlaveProcess
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveProcess);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    MasterSlaveProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~MasterSlaveProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    /// Called once before the solution loop and is writing the quadrature domain.
    void ExecuteBeforeSolutionLoop() override;

    const Parameters GetDefaultParameters() const override;

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

private:
    ///@name Member Variables
    ///@{

    /// Model part and different settings
    Model& mrModel;             /// The main model part
    Parameters mThisParameters; /// The parameters (can be used for general pourposes)
    ModelPart* mpModelPart;

    ///@}

}; // Class MasterSlaveProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MasterSlaveProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MasterSlaveProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_MASTER_SLAVE_PROCESS_H_INCLUDED  defined
