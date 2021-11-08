//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED)
#define KRATOS_RANS_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"

// Application includes
#include "rans_formulation_process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Clips given scalar variable to a range
 *
 * This process clips a given scalar variable to a range in all nodes in the model part.
 *
 */

class KRATOS_API(RANS_APPLICATION) RansClipScalarVariableProcess
: public RansFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansClipScalarVariableProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansClipScalarVariableProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansClipScalarVariableProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansClipScalarVariableProcess() override = default;

    /// Assignment operator.
    RansClipScalarVariableProcess& operator=(RansClipScalarVariableProcess const& rOther) = delete;

    /// Copy constructor.
    RansClipScalarVariableProcess(RansClipScalarVariableProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void Execute() override;

    void ExecuteAfterCouplingSolveStep() override
    {
        Execute();
    }

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    std::string mModelPartName;
    std::string mVariableName;
    int mEchoLevel;

    double mMinValue;
    double mMaxValue;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansClipScalarVariableProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED defined
