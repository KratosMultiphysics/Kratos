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

#if !defined(KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/communicator.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "linear_solvers/linear_solver.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/process.h"
#include "processes/variational_distance_calculation_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class KRATOS_API(RANS_APPLICATION) RansWallDistanceCalculationBaseProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using SparseSpaceType = TSparseSpace;
    using DenseSpaceType = TDenseSpace;
    using LinearSolverType = TLinearSolver;
    using BuilderSolverPointerType =
        typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer;

    using VariationalDistanceCalculationProcessType2D =
        VariationalDistanceCalculationProcess<2, SparseSpaceType, DenseSpaceType, LinearSolverType>;
    using VariationalDistanceCalculationProcessType3D =
        VariationalDistanceCalculationProcess<3, SparseSpaceType, DenseSpaceType, LinearSolverType>;

    /// Pointer definition of RansWallDistanceCalculationBaseProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansWallDistanceCalculationBaseProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansWallDistanceCalculationBaseProcess(
        Model& rModel,
        Parameters rParameters)
    : mrModel(rModel)
    {
        KRATOS_TRY

        rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        mMaxIterations = rParameters["max_iterations"].GetInt();
        mEchoLevel = rParameters["echo_level"].GetInt();
        mModelPartName = rParameters["model_part_name"].GetString();
        mWallFlagVariableName = rParameters["wall_flag_variable_name"].GetString();
        mWallFlagVariableValue = rParameters["wall_flag_variable_value"].GetBool();
        mRecalculateAtEachTimeStep =
            rParameters["re_calculate_at_each_time_step"].GetBool();
        mCorrectDistancesUsingNeighbors =
            rParameters["correct_distances_using_neighbors"].GetBool();
        mLinearSolverParameters = rParameters["linear_solver_settings"];

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansWallDistanceCalculationBaseProcess() override = default;

    /// Assignment operator.
    RansWallDistanceCalculationBaseProcess& operator=(RansWallDistanceCalculationBaseProcess const& rOther) = delete;

    /// Copy constructor.
    RansWallDistanceCalculationBaseProcess(RansWallDistanceCalculationBaseProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        KRATOS_TRY

        const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

        KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(DISTANCE))
            << "DISTANCE is not found in nodal solution step variables list of "
            << mModelPartName << ".";

        return 0.0;

        KRATOS_CATCH("");
    }

    void ExecuteInitialize() override
    {
        KRATOS_TRY

        this->CreateLinearSolver();
        this->CreateBuilderAndSolver();

        auto& r_model_part = mrModel.GetModelPart(mModelPartName);
        const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

        if (domain_size == 2) {
            mpVariationalDistanceCalculationProcess2D =
                Kratos::make_shared<VariationalDistanceCalculationProcessType2D>(
                    r_model_part, mpLinearSolver, mpBuilderAndSolver, mMaxIterations);
        } else if (domain_size == 3) {
            mpVariationalDistanceCalculationProcess3D =
                Kratos::make_shared<VariationalDistanceCalculationProcessType3D>(
                    r_model_part, mpLinearSolver, mpBuilderAndSolver, mMaxIterations);
        } else {
            KRATOS_ERROR << "Unknown domain size = " << domain_size;
        }

        CalculateWallDistances();

        KRATOS_CATCH("");
    }

    void ExecuteInitializeSolutionStep() override
    {
        if (mRecalculateAtEachTimeStep) {
            CalculateWallDistances();
        }
    }

    const Parameters GetDefaultParameters() const override
    {

        const auto default_parameters = Parameters(R"(
        {
            "model_part_name"                  : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "max_iterations"                   : 10,
            "echo_level"                       : 0,
            "wall_flag_variable_name"          : "STRUCTURE",
            "wall_flag_variable_value"         : true,
            "re_calculate_at_each_time_step"   : false,
            "correct_distances_using_neighbors": true,
            "linear_solver_settings" : {
                "solver_type"     : "amgcl"
            }
        })");

        return default_parameters;
    }

    ///@}
    ///@name Input and output
    ///@{

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

protected:
    ///@name Protected member Variables
    ///@{

    Model& mrModel;

    Parameters mLinearSolverParameters;
    std::string mModelPartName;
    typename TLinearSolver::Pointer mpLinearSolver;
    BuilderSolverPointerType mpBuilderAndSolver;

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Create a Linear Solver object
     *
     * This method should create linear solver in respective SparseSpace
     * to be used in Variational distance calculation process
     *
     */
    virtual void CreateLinearSolver() = 0;

    /**
     * @brief Create a Builder And Solver object
     *
     * This method should create builder and solver in respective SparseSpace
     * to be used in Variational distance calculation process.
     *
     */
    virtual void CreateBuilderAndSolver() = 0;

    ///@}

private:
    ///@name Member Variables
    ///@{

    int mMaxIterations;
    int mEchoLevel;
    std::string mWallFlagVariableName;
    bool mWallFlagVariableValue;
    bool mRecalculateAtEachTimeStep;
    bool mCorrectDistancesUsingNeighbors;

    typename VariationalDistanceCalculationProcessType2D::Pointer mpVariationalDistanceCalculationProcess2D;
    typename VariationalDistanceCalculationProcessType3D::Pointer mpVariationalDistanceCalculationProcess3D;

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateWallDistances()
    {
        KRATOS_TRY

        auto& r_model_part = mrModel.GetModelPart(mModelPartName);

        const Flags& r_wall_flag = KratosComponents<Flags>::Get(mWallFlagVariableName);

        VariableUtils variable_utilities;

        variable_utilities.SetVariable(DISTANCE, 1.0, r_model_part.Nodes());
        variable_utilities.SetVariable(DISTANCE, 0.0, r_model_part.Nodes(),
                                       r_wall_flag, mWallFlagVariableValue);

        const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

        if (domain_size == 2) {
            mpVariationalDistanceCalculationProcess2D->Execute();
        } else if (domain_size == 3) {
            mpVariationalDistanceCalculationProcess3D->Execute();
        } else {
            KRATOS_ERROR << "Unknown domain size = " << domain_size;
        }

        if (mCorrectDistancesUsingNeighbors) {
            CorrectWallDistances();
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Wall distances calculated in " << mModelPartName << ".\n";

        KRATOS_CATCH("");
    }

    /**
     * @brief Corrects wall distances
     *
     * VariationalDistanceCalculationProcess gives incorrect negative wall distances
     * when mesh is not refined enough. This method iterates through NEIGHBOUR nodes for a node with negative
     * DISTANCE and corrects it by calculated positive DISTANCE average from the NEIGHBOURs
     *
     * Works in MPI as well.
     *
     */
    void CorrectWallDistances()
    {
        KRATOS_TRY

        auto& r_model_part = mrModel.GetModelPart(mModelPartName);
        auto& r_communicator = r_model_part.GetCommunicator();
        auto& r_data_communicator = r_communicator.GetDataCommunicator();

        auto& r_nodes = r_communicator.LocalMesh().Nodes();

        VariableUtils().SetNonHistoricalVariableToZero(DISTANCE, r_model_part.Nodes());

        FindGlobalNodalNeighboursProcess find_nodal_neighbours_process(
            r_data_communicator, r_model_part);
        find_nodal_neighbours_process.Execute();

        class GlobalPointerAdder
        {
        public:
            typedef GlobalPointersVector<ModelPart::NodeType> value_type;
            value_type gp_vector;

            value_type GetValue()
            {
                gp_vector.Unique();
                return gp_vector;
            }

            void LocalReduce(const value_type& rGPVector)
            {
                for (auto& r_gp : rGPVector.GetContainer()) {
                    this->gp_vector.push_back(r_gp);
                }
            }
            void ThreadSafeReduce(GlobalPointerAdder& rOther)
            {
#pragma omp critical
                {
                    for (auto& r_gp : rOther.gp_vector.GetContainer()) {
                        this->gp_vector.push_back(r_gp);
                    }
                }
            }
        };

        auto all_global_pointers =
            BlockPartition<ModelPart::NodesContainerType>(r_nodes).for_each<GlobalPointerAdder>(
                [](ModelPart::NodeType& rNode) {
                    return rNode.GetValue(NEIGHBOUR_NODES);
                });

        GlobalPointerCommunicator<ModelPart::NodeType> pointer_comm(
            r_data_communicator, all_global_pointers);

        auto distance_proxy =
            pointer_comm.Apply([](const GlobalPointer<ModelPart::NodeType>& gp) -> double {
                return gp->FastGetSolutionStepValue(DISTANCE);
            });


        class CustomReducer
        {
        public:
            typedef int value_type;
            int number_of_modified_nodes;

            value_type GetValue()
            {
                return number_of_modified_nodes;
            }

            void LocalReduce(int ModifiedNodesCount)
            {
                this->number_of_modified_nodes += ModifiedNodesCount;
            }
            void ThreadSafeReduce(CustomReducer& rOther)
            {
#pragma omp critical
                {
                    this->number_of_modified_nodes += rOther.number_of_modified_nodes;
                }
            }
        };

        int number_of_modified_nodes =
            BlockPartition<ModelPart::NodesContainerType>(r_nodes).for_each<CustomReducer>(
                [&](ModelPart::NodeType& rNode) -> int {
                    if (rNode.FastGetSolutionStepValue(DISTANCE) < 0.0) {
                        const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES);

                        int count = 0;
                        double average_value = 0.0;
                        for (int j_node = 0;
                             j_node < static_cast<int>(r_neighbours.size()); ++j_node) {
                            const double j_distance =
                                distance_proxy.Get(r_neighbours(j_node));
                            if (j_distance > 0.0) {
                                average_value += j_distance;
                                count++;
                            }
                        }

                        if (count > 0) {
                            rNode.SetValue(
                                DISTANCE, average_value / static_cast<double>(count));
                        } else {
                            KRATOS_ERROR << "Node " << rNode.Id() << " at "
                                         << rNode.Coordinates() << " didn't find any neighbour("
                                         << r_neighbours.size() << ") with positive DISTANCE. Please recheck "
                                         << mModelPartName << ".\n";
                        }

                        return 1;
                    }
                    return 0;
                });

        BlockPartition<ModelPart::NodesContainerType>(r_nodes).for_each(
            [&](ModelPart::NodeType& rNode) {
                double& r_distance = rNode.FastGetSolutionStepValue(DISTANCE);
                const double avg_distance = rNode.GetValue(DISTANCE);
                if (r_distance < 0.0) {
                    r_distance = avg_distance;
                }
            });

        r_communicator.SynchronizeVariable(DISTANCE);
        number_of_modified_nodes =
            r_communicator.GetDataCommunicator().SumAll(number_of_modified_nodes);

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0 && number_of_modified_nodes > 0)
            << "Corrected " << number_of_modified_nodes
            << " nodal wall distances in " << mModelPartName << ".\n";

        KRATOS_CATCH("");
    }

    ///@}

}; // Class RansWallDistanceCalculationBaseProcess

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class KRATOS_API(RANS_APPLICATION) RansWallDistanceCalculationProcess
    : public RansWallDistanceCalculationBaseProcess<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    using SparseSpaceType = TSparseSpace;
    using DenseSpaceType = TDenseSpace;
    using LinearSolverType = TLinearSolver;
    using BuilderSolverPointerType =
        typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer;
    using BaseType =
        RansWallDistanceCalculationBaseProcess<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Pointer definition of RansWallDistanceCalculationBaseProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansWallDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansWallDistanceCalculationProcess(Model& rModel, Parameters rParameters)
        : BaseType(rModel, rParameters)
    {
    }

    /// Destructor.
    ~RansWallDistanceCalculationProcess() override = default;

    /// Assignment operator.
    RansWallDistanceCalculationProcess& operator=(RansWallDistanceCalculationProcess const& rOther) = delete;

    /// Copy constructor.
    RansWallDistanceCalculationProcess(RansWallDistanceCalculationProcess const& rOther) = delete;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}

private:
    ///@name Private Operations
    ///@{

    void CreateLinearSolver() override;

    void CreateBuilderAndSolver() override;

    ///@}
};

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
