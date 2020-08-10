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

/**
 * @class OutputQuadratureDomainProcess
 * @ingroup KratosCore
 * @brief This class is used in order to check results using a json file containing the solution a given model part with a certain frequency
 * @details This stores the dababase in a class denominated ResultDatabase which considers Table to store the information, therefore being able to interpolate results
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(KRATOS_CORE) OutputQuadratureDomainProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of OutputQuadratureDomainProcess
    KRATOS_CLASS_POINTER_DEFINITION(OutputQuadratureDomainProcess);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    /// The node type definiton
    typedef Node<3> NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    OutputQuadratureDomainProcess(
        Model& rModel,
        Parameters ThisParameters = Parameters(R"({})"));

    /// Destructor.
    virtual ~OutputQuadratureDomainProcess() {}

    ///@}
    ///@name Operations
    ///@{

    /// Called only once befor the solution loop
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
    Model& mrModel;         /// The main model part
    Parameters mThisParameters;     /// The parameters (can be used for general pourposes)

    ///@}

}; // Class OutputQuadratureDomainProcess

///@}

///@name Type Definitions
///@{


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
