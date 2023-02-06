//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "processes/process.h"

// Application includes

namespace Kratos {
///@addtogroup OptimizationApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ExtractBoundaryNodesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ExtractBoundaryNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(ExtractBoundaryNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    ExtractBoundaryNodesProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~ExtractBoundaryNodesProcess() override = default;

    /// Assignment operator.
    ExtractBoundaryNodesProcess& operator=(ExtractBoundaryNodesProcess const& rOther) = delete;

    /// Copy constructor.
    ExtractBoundaryNodesProcess(ExtractBoundaryNodesProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void Execute() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;
    std::string mModelPartName;
    std::string mBoundaryNodesSubModelPartName;

    int mEchoLevel;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ExtractBoundaryNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.
