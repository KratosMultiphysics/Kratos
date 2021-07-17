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

#if !defined(KRATOS_RANS_COMPUTE_PREVIOUS_STEP_BOSSAK_SECOND_DERIVATIVESPROCESS_H_INCLUDED)
#define KRATOS_RANS_COMPUTE_PREVIOUS_STEP_BOSSAK_SECOND_DERIVATIVES_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "containers/model.h"
#include "containers/variable.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application incldues

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Apply a specific flag for nodes and conditions
 *
 * This process apply a given flag to nodes of the modelpart.
 * Then, if preferred, applies to all conditions in the given model, which has
 * the given flag applied to all the nodes in the specific condition.
 *
 */

class KRATOS_API(RANS_APPLICATION) RansComputePreviousStepBossakSecondDerivativesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    using CopyVariableDataListItem = std::tuple<const std::string, const bool, const std::string, const bool>;

    /// Pointer definition of RansComputePreviousStepBossakSecondDerivativesProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansComputePreviousStepBossakSecondDerivativesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansComputePreviousStepBossakSecondDerivativesProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansComputePreviousStepBossakSecondDerivativesProcess() override = default;

    /// Assignment operator.
    RansComputePreviousStepBossakSecondDerivativesProcess& operator=(RansComputePreviousStepBossakSecondDerivativesProcess const& rOther) = delete;

    /// Copy constructor.
    RansComputePreviousStepBossakSecondDerivativesProcess(RansComputePreviousStepBossakSecondDerivativesProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void Execute() override;

    void ExecuteFinalizeSolutionStep() override;

    void ExecuteFinalize() override;

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
    ///@name Private Classes
    ///@{

    enum ExecutionPoint
    {
        INITIALIZE = 0,
        INITIALIZE_SOLUTION_STEP = 1,
        EXECUTE = 2,
        FINALIZE_SOLUTION_STEP = 3,
        FINALIZE = 4
    };

    template<class TVariableDataType>
    class BossakPreviousStepSecondDerivativeCalculator
    {
    public:
        ///@name Life Cycle
        ///@{

        explicit BossakPreviousStepSecondDerivativeCalculator(
            const double BossakAlpha,
            const Variable<TVariableDataType>& rCurrentSecondDerivative,
            const Variable<TVariableDataType>& rCurrentRelaxedSecondDerivative,
            const bool IsCurrentRelaxedSecondDerivativeHistorical);

        ///@}
        ///@name Public Operations
        ///@{

        void Check(const ModelPart& rModelPart) const;

        void Calculate(NodeType& rNode) const;

        std::string Info() const;

        ///@}

    private:
        ///@name Private Member Variables
        ///@{

        const double mBossakAlpha;
        const double mOneMinusBossakAlpha;
        const Variable<TVariableDataType>& mrSecondDerivative;
        const Variable<TVariableDataType>& mrRelaxedSecondDerivative;
        const bool mIsRelaxedSecondDerivativeHistorical;

        const TVariableDataType& (BossakPreviousStepSecondDerivativeCalculator::*mRelaxedSecondDerivativeGetterMethod)(const NodeType& rNode) const;

        const TVariableDataType& GetRelaxedSecondDerivativeFromHistorical(const NodeType& rNode) const;

        const TVariableDataType& GetRelaxedSecondDerivativeFromNonHistorical(const NodeType& rNode) const;

        ///@}
    };

    ///@}
    ///@name Private Member Variables
    ///@{

    Model& mrModel;
    int mEchoLevel;

    std::string mModelPartName;

    std::vector<ExecutionPoint> mExecutionPointsList;

    std::vector<BossakPreviousStepSecondDerivativeCalculator<double>> mDoubleVariableDataList;
    std::vector<BossakPreviousStepSecondDerivativeCalculator<array_1d<double, 3>>> mArray3DVariableDataList;

    ///@}
    ///@name Private Operations
    ///@{

    void ExecuteCalculation();

    void UpdateExecutionPointsList(const std::vector<std::string>& rCopyExecutionPointsList);

    ///@}

}; // Class RansComputePreviousStepBossakSecondDerivativesProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansComputePreviousStepBossakSecondDerivativesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_COMPUTE_PREVIOUS_STEP_BOSSAK_SECOND_DERIVATIVES_H_INCLUDED defined
