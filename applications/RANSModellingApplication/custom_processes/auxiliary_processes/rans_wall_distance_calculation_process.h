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
#include "custom_utilities/rans_variable_utils.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "factories/linear_solver_factory.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "processes/variational_distance_calculation_process.h"
#include "rans_modelling_application_variables.h"

namespace Kratos
{
///@addtogroup RANSModellingApplication
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
class RansWallDistanceCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansWallDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansWallDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansWallDistanceCalculationProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "max_iterations"           : 10,
            "echo_level"               : 0,
            "wall_flag_variable_name"  : "STRUCTURE",
            "wall_flag_variable_value" : true,
            "linear_solver_settings" : {
                "solver_type"     : "amgcl"
            }
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(
            mrParameters["linear_solver_settings"]);

        mMaxIterations = mrParameters["max_iterations"].GetInt();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mWallFlagVariableName = mrParameters["wall_flag_variable_name"].GetString();
        mWallFlagVariableValue = mrParameters["wall_flag_variable_value"].GetBool();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansWallDistanceCalculationProcess() override
    {
        // delete mpDistanceCalculator;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        KRATOS_TRY

        const ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
        }

        return 0;

        KRATOS_CATCH("");
    }

    void ExecuteInitialize() override
    {
        CalculateWallDistances();
    }

    void Execute() override
    {
        // CalculateWallDistances();
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
    Parameters& mrParameters;
    std::string mModelPartName;

    typename TLinearSolver::Pointer mpLinearSolver;

    int mMaxIterations;
    int mEchoLevel;
    std::string mWallFlagVariableName;
    bool mWallFlagVariableValue;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateWallDistances()
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const Flags& r_wall_flag = KratosComponents<Flags>::Get(mWallFlagVariableName);

        RansVariableUtils().SetScalarVar(DISTANCE, 1.0, r_model_part.Nodes());
        RansVariableUtils().SetScalarVarForFlag(
            DISTANCE, 0.0, r_model_part.Nodes(), r_wall_flag, mWallFlagVariableValue);

        const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

        if (domain_size == 2)
        {
            VariationalDistanceCalculationProcess<2, TSparseSpace, TDenseSpace, TLinearSolver> distance_calculation_process(
                r_model_part, mpLinearSolver, mMaxIterations);
            distance_calculation_process.Execute();
        }
        else if (domain_size == 3)
        {
            VariationalDistanceCalculationProcess<3, TSparseSpace, TDenseSpace, TLinearSolver> distance_calculation_process(
                r_model_part, mpLinearSolver, mMaxIterations);
            distance_calculation_process.Execute();
        }
        else
        {
            KRATOS_ERROR << "Unknown domain size = " << domain_size;
        }

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
