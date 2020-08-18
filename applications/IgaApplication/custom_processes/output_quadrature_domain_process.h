//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_OUTPUT_QUADRATURE_DOMAIN_PROCESS_H_INCLUDED )
#define  KRATOS_OUTPUT_QUADRATURE_DOMAIN_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"

#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class OutputQuadratureDomainProcess
 * @ingroup IgaApplication
 * @brief This class outputs the location of the quadrature points within the local space of the containing geometry. */
class KRATOS_API(IGA_APPLICATION) OutputQuadratureDomainProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of OutputQuadratureDomainProcess
    KRATOS_CLASS_POINTER_DEFINITION(OutputQuadratureDomainProcess);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    OutputQuadratureDomainProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~OutputQuadratureDomainProcess() = default;

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
        return "OutputQuadratureDomainProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "OutputQuadratureDomainProcess";
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

}; // Class OutputQuadratureDomainProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  OutputQuadratureDomainProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const OutputQuadratureDomainProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_OUTPUT_QUADRATURE_DOMAIN_PROCESS_H_INCLUDED  defined
