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

#if !defined(KRATOS_RANS_SMOOTH_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED)
#define KRATOS_RANS_SMOOTH_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"

// Application includes
#include "rans_point_execution_formulation_process.h"

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

class KRATOS_API(RANS_APPLICATION) RansSmoothClipScalarVariableProcess
: public RansPointExecutionFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansSmoothClipScalarVariableProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansSmoothClipScalarVariableProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansSmoothClipScalarVariableProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansSmoothClipScalarVariableProcess() override = default;

    /// Assignment operator.
    RansSmoothClipScalarVariableProcess& operator=(RansSmoothClipScalarVariableProcess const& rOther) = delete;

    /// Copy constructor.
    RansSmoothClipScalarVariableProcess(RansSmoothClipScalarVariableProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

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
    int mNumberOfSweeps;
    bool mAlwaysFindNeighbourNodes;
    bool mIsNeighbourNodesInitialized;

    double mMinValue;
    double mMaxValue;
    double mInverseDistanceWeightingPowerParameter;

    ///@}
    ///@name Private operations
    ///@{

    void FindNeighbourNodes();

    void ExecuteOperation() override;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansSmoothClipScalarVariableProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_SMOOTH_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED defined
