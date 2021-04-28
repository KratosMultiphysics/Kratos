//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_MPM_PROCESS_H_INCLUDED )
#define  KRATOS_MPM_PROCESS_H_INCLUDED

// System includes

// External includes

// Kratos includes
#include "containers/model.h"
#include "processes/process.h"

#include "utilities/quadrature_points_utility.h"
#include "utilities/variable_utils.h"

// Iga includes
#include "iga_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class OutputQuadratureDomainProcess
 * @ingroup IgaApplication
 * @brief This class outputs the location of the quadrature points within the local space of the containing geometry. */
class KRATOS_API(IGA_APPLICATION) MpmProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MpmProcess
    KRATOS_CLASS_POINTER_DEFINITION(MpmProcess);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    MpmProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~MpmProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    /// Called once before the solution loop and is writing the quadrature domain.
    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "MpmProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MpmProcess";
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

    ///@}

    void ResetNodalVariables(ModelPart& rModelPart);

}; // Class MpmProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
    MpmProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
    const MpmProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_MPM_PROCESS_H_INCLUDED  defined
