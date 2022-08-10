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

#if !defined(KRATOS_MAX_SINGULAR_VALUE_DECOMPOSITION_PROCESS_H_INCLUDED)
#define KRATOS_MAX_SINGULAR_VALUE_DECOMPOSITION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "includes/process_info.h"

// Application includes

namespace Kratos
{
///@addtogroup LinearSolversApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Clips given scalar variable to a range
 *
 * This process clips a given scalar variable to a range in all nodes in the model part.
 *
 */

class KRATOS_API(LINEAR_SOLVERS_APPLICATION) MaxSingularValueDecompositionProcess
: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MaxSingularValueDecompositionProcess
    KRATOS_CLASS_POINTER_DEFINITION(MaxSingularValueDecompositionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    MaxSingularValueDecompositionProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~MaxSingularValueDecompositionProcess() override = default;

    /// Assignment operator.
    MaxSingularValueDecompositionProcess& operator=(MaxSingularValueDecompositionProcess const& rOther) = delete;

    /// Copy constructor.
    MaxSingularValueDecompositionProcess(MaxSingularValueDecompositionProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

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
    std::string mInputVariableName;
    std::string mOutputVariableName;
    std::string mContainerType;
    int mEchoLevel;

    ///@}
    ///@name Private operations
    ///@{

    template<class TContainterType>
    void CalculateAndStoreMaxSingularValues(
        TContainterType& rContainer,
        const ProcessInfo& rProcessInfo) const;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const MaxSingularValueDecompositionProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_MAX_SINGULAR_VALUE_DECOMPOSITION_PROCESS_H_INCLUDED defined
