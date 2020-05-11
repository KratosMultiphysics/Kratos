//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "factories/linear_solver_factory.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "processes/process.h"
#include "utilities/variable_utils.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class KRATOS_API(RANS_APPLICATION) RansWallDistanceCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;
    using BuilderSolverPointerType =
        typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer;

    /// Pointer definition of RansWallDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansWallDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansWallDistanceCalculationProcess(Model& rModel, Parameters rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"               : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "max_iterations"                : 10,
            "echo_level"                    : 0,
            "wall_flag_variable_name"       : "STRUCTURE",
            "wall_flag_variable_value"      : true,
            "re_calculate_at_each_time_step": false,
            "linear_solver_settings" : {
                "solver_type"     : "amgcl"
            }
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        this->CreateLinearSolver();

        mMaxIterations = mrParameters["max_iterations"].GetInt();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mWallFlagVariableName = mrParameters["wall_flag_variable_name"].GetString();
        mWallFlagVariableValue = mrParameters["wall_flag_variable_value"].GetBool();
        mRecalculateAtEachTimeStep =
            mrParameters["re_calculate_at_each_time_step"].GetBool();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansWallDistanceCalculationProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        KRATOS_TRY

        RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

        const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, DISTANCE);

        return 0.0;

        KRATOS_CATCH("");
    }

    void SetBuilderAndSolver(BuilderSolverPointerType pBuilderAndSolver)
    {
        mpBuilderAndSolver = pBuilderAndSolver;
    }

    typename TLinearSolver::Pointer GetLinearSolver()
    {
        return mpLinearSolver;
    }

    void ExecuteInitialize() override
    {
        CalculateWallDistances();
    }

    void ExecuteInitializeSolutionStep() override
    {
        if (mRecalculateAtEachTimeStep)
        {
            CalculateWallDistances();
        }
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return std::string("RansWallDistanceCalculationProcess");
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Model& mrModel;
    Parameters mrParameters;
    std::string mModelPartName;

    typename TLinearSolver::Pointer mpLinearSolver;
    BuilderSolverPointerType mpBuilderAndSolver;

    int mMaxIterations;
    int mEchoLevel;
    std::string mWallFlagVariableName;
    bool mWallFlagVariableValue;
    bool mRecalculateAtEachTimeStep;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CreateLinearSolver()
    {
        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(
            mrParameters["linear_solver_settings"]);
    }

    void ExecuteVariationalDistanceCalculationProcess();

    void CalculateWallDistances()
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const Flags& r_wall_flag = KratosComponents<Flags>::Get(mWallFlagVariableName);

        VariableUtils variable_utilities;

        variable_utilities.SetVariable(DISTANCE, 1.0, r_model_part.Nodes());
        variable_utilities.SetVariable(DISTANCE, 0.0, r_model_part.Nodes(),
                                       r_wall_flag, mWallFlagVariableValue);

        ExecuteVariationalDistanceCalculationProcess();

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Wall distances calculated in " << mModelPartName << ".\n";

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RansWallDistanceCalculationProcess& operator=(RansWallDistanceCalculationProcess const& rOther);

    /// Copy constructor.
    RansWallDistanceCalculationProcess(RansWallDistanceCalculationProcess const& rOther);

    ///@}

}; // Class RansWallDistanceCalculationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED defined
