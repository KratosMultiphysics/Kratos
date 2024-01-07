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

#if !defined(KRATOS_RANS_BOSSAK_PREVIOUS_STEP_BOSSAK_VARIABLE_DERIVATIVES_PROCESS_H_INCLUDED)
#define KRATOS_RANS_BOSSAK_PREVIOUS_STEP_BOSSAK_VARIABLE_DERIVATIVES_PROCESS_H_INCLUDED

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

class KRATOS_API(RANS_APPLICATION) RansInitializeBossakPreviousStepVariableDerivatives : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansInitializeBossakPreviousStepVariableDerivatives
    KRATOS_CLASS_POINTER_DEFINITION(RansInitializeBossakPreviousStepVariableDerivatives);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansInitializeBossakPreviousStepVariableDerivatives(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansInitializeBossakPreviousStepVariableDerivatives() override = default;

    /// Assignment operator.
    RansInitializeBossakPreviousStepVariableDerivatives& operator=(RansInitializeBossakPreviousStepVariableDerivatives const& rOther) = delete;

    /// Copy constructor.
    RansInitializeBossakPreviousStepVariableDerivatives(RansInitializeBossakPreviousStepVariableDerivatives const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitializeSolutionStep() override;

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

    struct BossakConstants
    {
        double Alpha;
        double Gamma;
        double C0;
        double C2;
        double C3;
        double C4;
    };

    template<class TVariableDataType>
    class BossakPreviousStepDerivativeCalculator
    {
    public:
        ///@name Life Cycle
        ///@{

        explicit BossakPreviousStepDerivativeCalculator(
            const Variable<TVariableDataType>& rFirstDerivativeVariable,
            const Variable<TVariableDataType>& rSecondDerivativeVariable,
            const Variable<TVariableDataType>& rRelaxedSecondDerivativeVariable,
            const bool IsRelaxedSecondDerivativeVariableHistorical);

        ///@}
        ///@name Public Operations
        ///@{

        void Check(const ModelPart& rModelPart) const;

        void Calculate(
            NodeType& rNode,
            const BossakConstants& rBossakConstants) const;

        std::string Info() const;

        ///@}

    private:
        ///@name Private Member Variables
        ///@{

        const Variable<TVariableDataType>& mrFirstDerivativeVariable;
        const Variable<TVariableDataType>& mrSecondDerivativeVariable;
        const Variable<TVariableDataType>& mrRelaxedSecondDerivativeVariable;

        const bool mIsRelaxedSecondDerivativeVariableHistorical;

        const TVariableDataType& (BossakPreviousStepDerivativeCalculator::*mRelaxedSecondDerivativeGetterMethod)(const NodeType& rNode) const;

        const TVariableDataType& GetRelaxedSecondDerivativeFromHistorical(const NodeType& rNode) const;

        const TVariableDataType& GetRelaxedSecondDerivativeFromNonHistorical(const NodeType& rNode) const;

        ///@}
    };

    ///@}
    ///@name Private Member Variables
    ///@{

    Model& mrModel;
    std::string mModelPartName;

    int mEchoLevel;

    double mBossakAlpha;
    bool mIsInitialized = false;

    std::vector<BossakPreviousStepDerivativeCalculator<double>> mDoubleVariableDataList;
    std::vector<BossakPreviousStepDerivativeCalculator<array_1d<double, 3>>> mArray3DVariableDataList;

    ///@}
    ///@name Private Operations
    ///@{

    void ExecuteCalculation();

    template<class TVariableDataType>
    bool CheckVariableTypes(
        const std::string& rFirstDerivativeVariableName,
        const std::string& rSecondDerivativeVariableName,
        const std::string& rRelaxedSecondDerivativeVariableName) const
    {
        KRATOS_TRY

        if (KratosComponents<Variable<TVariableDataType>>::Has(rFirstDerivativeVariableName)) {
            if (KratosComponents<Variable<TVariableDataType>>::Has(rSecondDerivativeVariableName)) {
                if (KratosComponents<Variable<TVariableDataType>>::Has(rRelaxedSecondDerivativeVariableName)) {
                    return true;
                } else {
                    KRATOS_ERROR
                        << "First derivative variable " << rFirstDerivativeVariableName
                        << " and relaxed second derivative variable "
                        << rRelaxedSecondDerivativeVariableName << " type mismatch.\n";
                }
            } else {
                KRATOS_ERROR << "First derivative variable " << rFirstDerivativeVariableName
                            << " and second derivative variable "
                            << rSecondDerivativeVariableName << " type mismatch.\n";
            }
        }

        return false;

        KRATOS_CATCH("");
    }

    ///@}

}; // Class RansInitializeBossakPreviousStepVariableDerivatives

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansInitializeBossakPreviousStepVariableDerivatives& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_COMPUTE_PREVIOUS_STEP_BOSSAK_SECOND_DERIVATIVES_H_INCLUDED defined
