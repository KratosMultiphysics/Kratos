//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <numeric>
#include <variant>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "containers/edge_based_data_structure.h"
#include "includes/define.h"
#include "processes/process.h"
#include "solving_strategies/strategies/butcher_tableau.h"
#include "utilities/atomic_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{
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

/// @brief Process to solve a pure convection problem using a Flux Corrected Transport scheme and an edge-based data structure
/// @tparam TDim Number of dimensions
template<unsigned int TDim>
class FluxCorrectedTransportConvectionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    using VariantTableauType = std::variant<ButcherTableauForwardEuler, ButcherTableauMidPointMethod, ButcherTableauRK3TVD, ButcherTableauRK4>;

    ///@}
    ///@name Pointer Definitions
    ///@{

    /// Pointer definition of FluxCorrectedTransportConvectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(FluxCorrectedTransportConvectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    FluxCorrectedTransportConvectionProcess(
        Model& rModel,
        Parameters ThisParameters)
        : Process()
    {
        // Validate the common settings as well as the element formulation specific ones
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Reference to the model part to which the convection is to be performed
        mpModelPart = &(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()));

        // Set the edge-based data structure pointer
        mpEdgeDataStructure = Kratos::make_unique<EdgeBasedDataStructure<TDim>>();

        // Assign all the required member variables
        mEchoLevel = ThisParameters["echo_level"].GetInt();
        mTimeSchemeName = ThisParameters["time_scheme"].GetString();
        mMaxAllowedCFL = ThisParameters["max_CFL"].GetDouble();
        mMaxAllowedDt = ThisParameters["max_delta_time"].GetDouble();
        mDiffusionConstant = ThisParameters["diffusion_constant"].GetDouble();
        mpConvectedVar = &KratosComponents<Variable<double>>::Get(ThisParameters["convected_variable_name"].GetString());
        mpConvectionVar = &KratosComponents<Variable<array_1d<double, 3>>>::Get(ThisParameters["convection_variable_name"].GetString());

        // Check provided values
        KRATOS_ERROR_IF(mDiffusionConstant < 1.0e-12 || mDiffusionConstant > 1.0)
            << "Provided 'diffusion_constant' " << mDiffusionConstant << " is not valid. Value must be between 0 and 1." << std::endl;
    }

    /// Copy constructor.
    FluxCorrectedTransportConvectionProcess(FluxCorrectedTransportConvectionProcess const& rOther) = delete;

    /// Destructor.
    ~FluxCorrectedTransportConvectionProcess() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    FluxCorrectedTransportConvectionProcess &operator=(FluxCorrectedTransportConvectionProcess const &rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override
    {
        // Check that there are DOFs for the convected variable
        // These are required to check the fixity of the convected variable
        block_for_each(mpModelPart->Nodes(), [this](auto& rNode){
            KRATOS_ERROR_IF_NOT(rNode.HasDofFor(*mpConvectedVar))
                << "No DOF for convected variable '" << mpConvectedVar->Name() << "' in node " << rNode.Id() << std::endl;
        });

        // Fill the edge-based data structure
        mpEdgeDataStructure->CalculateEdgeDataStructure(*mpModelPart);
        KRATOS_INFO_IF("FluxCorrectedTransportConvectionProcess", mEchoLevel > 0) << "Edge-based data structure completed." << std::endl;

        // Auxiliary flag to avoid resetting the edge data structure at each Execute call
        mPerformInitialize = false;

        // Allocate auxiliary arrays
        // Note that we give the size such that the index matches the id of the corresponding node
        // Also note that this implies having a "fake" 0-id entry in the first position of the arrays
        mAuxSize = mpModelPart->NumberOfNodes() + 1;
        mSolution.resize(mAuxSize);
        mSolutionOld.resize(mAuxSize);
        mLowOrderUpdate.resize(mAuxSize);
        mHighOrderUpdate.resize(mAuxSize);
        mConvectionValues.resize(mAuxSize);
    }

    void Execute() override
    {
        KRATOS_TRY;

        // Check if the edge data structure is already set
        if (mPerformInitialize) {
            this->ExecuteInitialize();
        }

        // Evaluate current time (n) and target time (n+1)
        // Note that it is assumed that the CloneTimeStep of the model part has been already done
        const auto& r_process_info = mpModelPart->GetProcessInfo();
        const auto &r_prev_process_info = r_process_info.GetPreviousTimeStepInfo();
        double prev_time = r_prev_process_info[TIME];
        const double target_time = r_process_info[TIME];

        // Initialize data vectors
        const SizeType n_nodes = mpModelPart->NumberOfNodes();
        IndexPartition<IndexType>(n_nodes).for_each([this](IndexType iNode){
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();
            mSolutionOld[i_node_id] = it_node->GetSolutionStepValue(*mpConvectedVar,1); //solution at n
            mConvectionValues[i_node_id] = it_node->GetSolutionStepValue(*mpConvectionVar,1); //convective velocity at n
        });

        // Get maximum time increment according to global CFL
        const double target_dt = target_time - prev_time; // Time to be covered by the substepping loop
        const double cfl_dt = this->CalculateSubStepDeltaTime(); // Maximum time increment according to edge mesh CFL
        const double max_dt = std::min(cfl_dt, mMaxAllowedDt); // Check if the user defined delta time is more restrictive than the CFL one
        const SizeType n_steps = std::ceil(target_dt / max_dt); // Required number of substeps (rounded up)
        const double dt = target_dt / n_steps; // Time increment to be used in the substep loop
        KRATOS_INFO_IF("FluxCorrectedTransportConvectionProcess", mEchoLevel > 0) << u8"Solving FCT convection with \u0394t = " << dt << u8" (max allowed \u0394t = " <<  max_dt << ")" << std::endl;

        // Substepping time loop
        for (IndexType step = 1; step <= n_steps; ++step) {
            // Solve current substep
            KRATOS_INFO_IF("FluxCorrectedTransportConvectionProcess", mEchoLevel > 1) << "Substep = " << step << " - Time = " <<  prev_time << u8" - \u0394t = " << dt << std::endl;
            this->SolveSubStep(dt);

            // Advance in time
            prev_time += dt;
        }

        // Set final solution in the model part database
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
            auto it_node = mpModelPart->NodesBegin() + iNode;
            if (!it_node->IsFixed(*mpConvectedVar)) {
                it_node->FastGetSolutionStepValue(*mpConvectedVar) = mSolution[it_node->Id()];
            }
        });

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        mSolution.clear();
        mSolutionOld.clear();
        mLowOrderUpdate.clear();
        mHighOrderUpdate.clear();
        mConvectionValues.clear();
        mpEdgeDataStructure->Clear();
    }

    const Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"({
            "echo_level" : 0,
            "model_part_name" : "",
            "convected_variable_name" : "DISTANCE",
            "convection_variable_name" : "VELOCITY",
            "max_CFL" : 1.0,
            "max_delta_time" : 1.0,
            "diffusion_constant" : 1.0,
            "time_scheme" : "RK3-TVD"
        })");

        return default_parameters;
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
        return "FluxCorrectedTransportConvectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FluxCorrectedTransportConvectionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart; /// Pointer to the background mesh model part

    std::string mTimeSchemeName; /// String containing the time integration scheme name

    SizeType mAuxSize; /// Size of the internal vector containers that equals the number of nodes plus one

    SizeType mEchoLevel; /// Level of information that is output

    double mMaxAllowedDt; /// Maximum allowed time step for the substepping

    double mMaxAllowedCFL; /// Maximum allowed CFL number for the substepping

    double mDiffusionConstant; /// Diffusion constant for the artificial diffusion in the low order scheme

    bool mPerformInitialize = true; /// Flag to indicate if the ExecuteInitialize is required to be done

    const Variable<double>* mpConvectedVar = nullptr; /// Pointer to the convected variable (e.g. DISTANCE)

    const Variable<array_1d<double,3>>* mpConvectionVar = nullptr; /// Pointer to the convection variable (e.g. VELOCITY)

    std::vector<double> mSolution; /// Auxiliary vector containing the current solution

    std::vector<double> mSolutionOld; /// Auxiliary vector containing the previous step solution

    std::vector<double> mLowOrderUpdate; /// Auxiliary vector containing the low order solution update

    std::vector<double> mHighOrderUpdate; /// Auxiliary vector containing the high order solution update

    std::vector<array_1d<double,3>> mConvectionValues; /// Auxiliary vector to store the convection variable values (e.g. VELOCITY)

    typename EdgeBasedDataStructure<TDim>::UniquePointer mpEdgeDataStructure; /// Pointer to the edge-based data structure

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    //TODO: Do it with the edge projected velocity
    double CalculateSubStepDeltaTime()
    {
        // Get edge data structure containers
        const auto& r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto& r_col_indices = mpEdgeDataStructure->GetColIndices();

        // Calculate the CFL number at each edge to get the minimum allowed time step
        // Note that we don't loop the last entry of the row container as it is NNZ
        SizeType n_nodes = r_row_indices.size() - 1;
        const double min_dt = IndexPartition<IndexType>(n_nodes).for_each<MinReduction<double>>([&](IndexType iNode){
            double dt_ij = std::numeric_limits<double>::max();
            const auto it_row = r_row_indices.begin() + iNode;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;
            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // i-node nodal data
                const double i_vel_norm = norm_2(mConvectionValues[iNode]);

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double j_vel_norm = norm_2(mConvectionValues[j_node_id]);

                    // Get ij-edge length from CSR data structure
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iNode, j_node_id);
                    const double ij_length = r_ij_edge_data.GetLength();

                    // Calculate the max time step from the velocities at both edge ends and get the minimum
                    dt_ij = CalculateEdgeLocalDeltaTime(ij_length, i_vel_norm, j_vel_norm);
                }
            }

            return dt_ij;
        });

        // Synchronize among MPI nodes
        mpModelPart->GetCommunicator().GetDataCommunicator().MinAll(min_dt);

        // Return maximum allowable time step
        return min_dt;
    }

    double CalculateEdgeLocalDeltaTime(
        const double Length,
        const double NormVelI,
        const double NormVelJ)
    {
        const double dx_CFL = Length * mMaxAllowedCFL;
        const double dt_i = NormVelI > 1.0e-12 ? dx_CFL / NormVelI : mMaxAllowedDt;
        const double dt_j = NormVelJ > 1.0e-12 ? dx_CFL / NormVelJ : mMaxAllowedDt;
        return std::min(dt_i, dt_j);
    }

    void SolveSubStep(const double DeltaTime)
    {
        // Calculate the low order Runge-Kutta update
        const double low_order_max_iteration = 1;
        const double low_order_diff_constant = mDiffusionConstant;

        // Calculate the high order Runge-Kutta update
        const double high_order_max_iteration = 3;
        const double high_order_diff_constant = 0.0;

        // Get Butcher tableau variant with the selected Runge-Kutta scheme
        const auto butcher_tableau_variant = GetButcherTableauVariant();

        // Allocate auxiliary Runge-Kutta arrays
        const SizeType rk_steps = std::visit([](const auto& rTableau){return rTableau.SubstepCount();}, butcher_tableau_variant);
        std::vector<std::vector<double>> rk_residuals_lo(mAuxSize, std::vector<double>(rk_steps, 0.0));
        std::vector<std::vector<double>> rk_residuals_ho(mAuxSize, std::vector<double>(rk_steps, 0.0));

        // Initialize intermediate substep solution vector
        std::vector<double> rk_u_lo = mSolutionOld;
        std::vector<double> rk_u_ho = mSolutionOld;

        // Perform the Runge-Kutta intermediate steps
        for (IndexType step = 1; step <= rk_steps; ++step) {
            // Solve current Runge-Kutta step
            const IndexType residual_column = step - 1;
            const double rk_step_theta = std::visit([&step](const auto &rTableau) {return rTableau.GetIntegrationTheta(step);}, butcher_tableau_variant); // Runge-Kutta step integration theta

            // Calculate low order residual
            CalculateResidual(rk_residuals_lo, rk_u_lo, residual_column, low_order_diff_constant, rk_step_theta * DeltaTime);

            // Calculate high order residual
            CalculateResidual(rk_residuals_ho, rk_u_ho, residual_column, high_order_diff_constant, rk_step_theta * DeltaTime);

            if (step != rk_steps) {
                // Do the explicit lumped mass matrix solve to obtain the intermediate solutions
                const auto [rk_mat_row_begin, rk_mat_row_end] = std::visit([&step](const auto& rTableau) {return rTableau.GetMatrixRow(step);}, butcher_tableau_variant); // Runge-Kutta matrix row
                ExplicitSolveSolution(rk_u_lo, rk_residuals_lo, DeltaTime, rk_mat_row_begin, rk_mat_row_end, low_order_max_iteration);
                ExplicitSolveSolution(rk_u_ho, rk_residuals_ho, DeltaTime, rk_mat_row_begin, rk_mat_row_end, high_order_max_iteration);

                // Set the auxiliary arrays required to calculate the antidiffusive fluxes
                IndexPartition<IndexType>(mAuxSize).for_each([&](IndexType i){
                    mSolution[i] = rk_u_lo[i];
                    mHighOrderUpdate[i] = rk_u_ho[i] - mSolutionOld[i];
                });

                // Current substep time advance coefficient to calculate the limiting
                double time_coeff = 0.0;
                for (auto it = rk_mat_row_begin; it != rk_mat_row_end; ++it) {
                    time_coeff += *it;
                }

                // Calculate the antidiffusive edge contributions
                CalculateAntidiffusiveEdgeContributions(time_coeff * DeltaTime);

                // Evaluate limiter
                EvaluateLimiter(time_coeff * DeltaTime);

                // Set low and high order solutions to the limited one
                rk_u_lo = mSolution;
                rk_u_ho = mSolution;
            }
        }

        // Do the final solution update
        const auto rk_weights = std::visit([](const auto& rTableau) -> auto {return rTableau.GetWeights();}, butcher_tableau_variant);
        ExplicitSolveSolution(mSolution, rk_residuals_lo, DeltaTime, rk_weights.begin(), rk_weights.end(), low_order_max_iteration);
        ExplicitSolveUpdate(mHighOrderUpdate, rk_residuals_ho, DeltaTime, rk_weights.begin(), rk_weights.end(), high_order_max_iteration);

        // Calculate the antidiffusive edge contributions
        CalculateAntidiffusiveEdgeContributions(DeltaTime);

        // Evaluate limiter
        EvaluateLimiter(DeltaTime);

        // Do the current step solution update
        IndexPartition<IndexType>(mAuxSize).for_each([this](IndexType i){
            mSolutionOld[i] = mSolution[i];
        });
    }

    auto GetButcherTableauVariant()
    {
        if (mTimeSchemeName == "forward_euler") {
            return VariantTableauType(ButcherTableauForwardEuler());
        } else if (mTimeSchemeName == "midpoint_method") {
            return VariantTableauType(ButcherTableauMidPointMethod());
        } else if (mTimeSchemeName == "RK3-TVD") {
            return VariantTableauType(ButcherTableauRK3TVD());
        } else if (mTimeSchemeName == "RK4") {
            return VariantTableauType(ButcherTableauRK4());
        } else {
            KRATOS_ERROR << "Butcher tableau cannot be found for time scheme '" << mTimeSchemeName << "'. Available options are 'forward_euler','midpoint_method','RK3-TVD' and 'RK4'." << std::endl;
        }
    }

    template<class TSolveCoeffsIterType>
    void ExplicitSolveSolution(
        std::vector<double>& rSolution,
        const std::vector<std::vector<double>>& rResiduals,
        const double DeltaTime,
        const TSolveCoeffsIterType& rSolveCoeffBegin,
        const TSolveCoeffsIterType& rSolveCoeffEnd,
        const IndexType MaxIteration)
    {
        auto solution_updater = [&](std::vector<double>& rSolution, const IndexType NodeId, const double UpdateValue)
        { rSolution[NodeId] = mSolutionOld[NodeId] + UpdateValue; };
        ExplicitSolve(rSolution, rResiduals, DeltaTime, rSolveCoeffBegin, rSolveCoeffEnd, solution_updater, MaxIteration);
    }

    template<class TSolveCoeffsIterType>
    void ExplicitSolveUpdate(
        std::vector<double>& rUpdate,
        const std::vector<std::vector<double>>& rResiduals,
        const double DeltaTime,
        const TSolveCoeffsIterType& rSolveCoeffBegin,
        const TSolveCoeffsIterType& rSolveCoeffEnd,
        const IndexType MaxIteration)
    {
        auto solution_updater = [&](std::vector<double> &rUpdate, const IndexType NodeId, const double UpdateValue)
        { rUpdate[NodeId] = UpdateValue; };
        ExplicitSolve(rUpdate, rResiduals, DeltaTime, rSolveCoeffBegin, rSolveCoeffEnd, solution_updater, MaxIteration);
    }

    template<class TSolveCoeffsIterType>
    void ExplicitSolve(
        std::vector<double>& rSolution,
        const std::vector<std::vector<double>>& rResiduals,
        const double DeltaTime,
        const TSolveCoeffsIterType& rSolveCoeffBegin,
        const TSolveCoeffsIterType& rSolveCoeffEnd,
        const std::function<void(std::vector<double>&, const IndexType, const double)>& rUpdater,
        const IndexType MaxIteration)
    {
        if (MaxIteration == 1) {
            // Standard lumped mass matrix solve
            IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
                // Get nodal data
                const auto it_node = mpModelPart->NodesBegin() + iNode;
                const IndexType i_node_id = it_node->Id();

                // Do the explicit lumped mass matrix solve
                const double M_l = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(i_node_id);
                const bool is_fixed = it_node->IsFixed(*mpConvectedVar);
                const double update = is_fixed ? 0.0 : (DeltaTime / M_l) * std::inner_product(rSolveCoeffBegin, rSolveCoeffEnd, rResiduals[i_node_id].begin(), 0.0);
                rUpdater(rSolution, i_node_id, update);
            });
        } else {
            // Get edge data structure containers
            const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
            const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
            SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

            // First explicit lumped mass matrix solve
            // Note that in the case of iterative explicit solve this initializes the update and Galerkin residual vectors
            std::vector<double> aux_update(mAuxSize);
            std::vector<double> galerkin_residual(mAuxSize);
            IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
                // Get nodal data
                const auto it_node = mpModelPart->NodesBegin() + iNode;
                const IndexType i_node_id = it_node->Id();

                // Do the explicit lumped mass matrix solve
                if (it_node->IsFixed(*mpConvectedVar)) {
                    aux_update[i_node_id] = 0.0;
                } else {
                    const double M_l = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(i_node_id);
                    const double i_node_res = std::inner_product(rSolveCoeffBegin, rSolveCoeffEnd, rResiduals[i_node_id].begin(), 0.0);
                    galerkin_residual[i_node_id] = i_node_res;
                    aux_update[i_node_id] = (DeltaTime / M_l) * i_node_res;
                }
            });

            // Do the consistent mass matrix solve (iteratively with lumped mass matrix)
            std::vector<double> it_residual(mAuxSize);
            for (IndexType i = 1; i < MaxIteration; ++i) {
                // Initialize current iteration residual to the Galerkin one
                IndexPartition<IndexType>(mAuxSize).for_each([&](IndexType i){
                    it_residual[i] = galerkin_residual[i];
                });

                // Add the lumping correction
                IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
                    // Get current row (node) storage data
                    const auto it_row = r_row_indices.begin() + iRow;
                    const IndexType i_col_index = *it_row;
                    const SizeType n_cols = *(it_row+1) - i_col_index;

                    // Check that there are CSR columns (i.e. that current node involves an edge)
                    if (n_cols != 0) {
                        // i-node nodal data
                        const double delta_u_h_i = aux_update[iRow];

                        // j-node nodal loop (i.e. loop ij-edges)
                        const auto i_col_begin = r_col_indices.begin() + i_col_index;
                        for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                            // j-node nodal data
                            IndexType j_node_id = *(i_col_begin + j_node);
                            const double delta_u_h_j = aux_update[j_node_id];

                            // Add previous iteration correction
                            const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                            const double Mc_ij = r_ij_edge_data.GetOffDiagonalConsistentMass();
                            const double res_edge_i = Mc_ij * (delta_u_h_i - delta_u_h_j) / DeltaTime;
                            const double res_edge_j = Mc_ij * (delta_u_h_j - delta_u_h_i) / DeltaTime;

                            // Atomic additions
                            AtomicAdd(it_residual[iRow], res_edge_i);
                            AtomicAdd(it_residual[j_node_id], res_edge_j);
                        }
                    }
                });

                // Solve current iteration residual
                IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
                    // Get nodal data
                    const auto it_node = mpModelPart->NodesBegin() + iNode;
                    const IndexType i_node_id = it_node->Id();

                    // Do the explicit lumped mass matrix solve
                    if (it_node->IsFixed(*mpConvectedVar)) {
                        aux_update[i_node_id] = 0.0;
                    } else {
                        const double M_l = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(i_node_id);
                        aux_update[i_node_id] = (DeltaTime / M_l) * it_residual[i_node_id];
                    }
                });
            }

            // Add the iterative update to the solution vector
            IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
                const auto it_node = mpModelPart->NodesBegin() + iNode;
                if (!it_node->IsFixed(*mpConvectedVar)) {
                    const IndexType i_node_id = it_node->Id();
                    rUpdater(rSolution, i_node_id, aux_update[i_node_id]);
                }
            });
        }
    }

    void CalculateResidual(
        std::vector<std::vector<double>>& rResiduals,
        const std::vector<double>& rSolutionVector,
        const IndexType ResidualColumn,
        const double DiffusionConstant,
        const double DeltaTime)
    {
        // Get edge data structure containers
        const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Aux TLS container
        struct AuxTLS
        {
            array_1d<double, TDim> F_i;
            array_1d<double, TDim> F_j;
            array_1d<double, TDim> d_ij;
            array_1d<double, TDim> b_ij;
            array_1d<double, TDim> F_ij_num;
            // array_1d<double, 3> vel_ij_half;
        };

        // Off-diagonal (edge) contributions assembly
        IndexPartition<IndexType>(aux_n_rows).for_each(AuxTLS(), [&](IndexType iRow, AuxTLS& rTLS){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // Get TLS data
                auto& F_i = rTLS.F_i;
                auto& F_j = rTLS.F_j;
                auto& d_ij = rTLS.d_ij;
                auto& b_ij = rTLS.b_ij;
                auto& F_ij_num = rTLS.F_ij_num;
                // auto& vel_ij_half = rTLS.vel_ij_half;

                // i-node nodal data
                const double u_i = rSolutionVector[iRow];
                const auto& r_i_vel = mConvectionValues[iRow];

                // i-node convective flux calculation
                for (IndexType d = 0; d < TDim; ++d) {
                    F_i[d] = u_i * r_i_vel[d];
                }

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double u_j = rSolutionVector[j_node_id];
                    const auto& r_j_vel = mConvectionValues[j_node_id];

                    // j-node convective flux calculation
                    for (IndexType d = 0; d < TDim; ++d) {
                        F_j[d] = u_j * r_j_vel[d];
                    }

                    // Get ij-edge operators data from CSR data structure
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const auto& r_Ni_DNj = r_ij_edge_data.GetOffDiagonalConvective();
                    const auto& r_DNi_Nj = r_ij_edge_data.GetOffDiagonalConvectiveTranspose();
                    d_ij = 0.5 * (r_Ni_DNj - r_DNi_Nj);
                    double D_ij = 0.0;
                    for (IndexType d = 0; d < TDim; ++d) {
                        D_ij += std::pow(d_ij[d],2);
                    }

                    // // Calculate numerical flux "upwind" contribution at the edge midpoint
                    // // Note that this numerical flux corresponds to the Lax-Wendroff scheme
                    // const double u_ij_half = 0.5 * (u_i + u_j) - 0.5 * DeltaTime * inner_prod(-d_ij, F_i - F_j) / D_ij;
                    // noalias(vel_ij_half) = 0.5 * (r_i_vel + r_j_vel);
                    // for (IndexType d = 0; d < TDim; ++d) {
                    //     F_ij_num[d] = 2.0 * (vel_ij_half[d] * u_ij_half);
                    // }

                    // Standard flux
                    F_ij_num = F_i + F_j;

                    // Calculate convection volume residual contributions
                    double res_edge_i = inner_prod(-d_ij, F_ij_num);
                    double res_edge_j = inner_prod(d_ij, F_ij_num);

                    // If current edge belogs to a boundary, add the corresponding convection boundary integrals
                    if (r_ij_edge_data.IsBoundary()) {
                        // Get ij-edge boundary operators from CSR data structure
                        const auto &r_N_N_normal = r_ij_edge_data.GetOffDiagonalConvectiveBoundary();
                        b_ij = 0.5 * r_N_N_normal;

                        // Add boundary contribution to the residual
                        const double res_edge_bd = inner_prod(b_ij, F_ij_num);
                        res_edge_i -= res_edge_bd;
                        res_edge_j += res_edge_bd;
                    }

                    // // Laplacian contribution for the testing
                    // const double c_edge = r_ij_edge_data.GetOffDiagonalLaplacian();
                    // const double res_edge_i = 1.0 * c_edge * (u_i - u_j);
                    // const double res_edge_j = 1.0 * c_edge * (u_j - u_i);

                    // Add low order scheme diffusion
                    if (DiffusionConstant > 0.0) {
                        const double local_dt = CalculateEdgeLocalDeltaTime(r_ij_edge_data.GetLength(), norm_2(r_i_vel), norm_2(r_j_vel));
                        const double c_tau = mDiffusionConstant / local_dt;
                        const double Mc_ij = r_ij_edge_data.GetOffDiagonalConsistentMass();
                        const double ij_low_order_diff = c_tau * Mc_ij;
                        res_edge_i += ij_low_order_diff * (u_j - u_i);
                        res_edge_j += ij_low_order_diff * (u_i - u_j);
                    }

                    // Atomic additions
                    AtomicAdd(rResiduals[iRow][ResidualColumn], res_edge_i);
                    AtomicAdd(rResiduals[j_node_id][ResidualColumn], res_edge_j);
                }
            }
        });

        // Add the diagonal boundary term coming from the convective flux split
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each(array_1d<double,TDim>(), [&](IndexType iNode, array_1d<double,TDim>& rTLS){
            // Get nodal data
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();

            // Add the diagonal contribution to the residual
            rTLS = mpEdgeDataStructure->GetBoundaryMassMatrixDiagonal(i_node_id);
            double aux_res = 0.0;
            const double u_i = rSolutionVector[i_node_id];
            const auto &r_i_vel = mConvectionValues[i_node_id];
            for (IndexType d = 0; d < TDim; ++d) {
                aux_res += u_i * r_i_vel[d] * rTLS[d];
            }
            rResiduals[i_node_id][ResidualColumn] -= aux_res;
        });
    }

    void CalculateAntidiffusiveEdgeContributions(const double DeltaTime)
    {
        // Check input time increment
        KRATOS_ERROR_IF(DeltaTime < 1.0e-12)
            << "Antidiffusive edge contributions cannot be computed as provided time increment " << DeltaTime << "is close to zero." << std::endl;

        // Get edge data structure containers
        const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Add the diffusion to the Taylor-Galerkin already computed residual
        IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // i-node nodal data
                const double u_i = mSolution[iRow]; // i-node low order solution (prediction)
                const double du_h_i = mHighOrderUpdate[iRow]; // i-node high order solution update
                const auto& r_i_vel = mConvectionValues[iRow];

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double u_j = mSolution[j_node_id]; // j-node low order solution (prediction)
                    const double du_h_j = mHighOrderUpdate[j_node_id]; // j-node high order solution update
                    const auto& r_j_vel = mConvectionValues[j_node_id];

                    // Calculate and store the antidiffusive flux edge contribution
                    auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const double local_dt = CalculateEdgeLocalDeltaTime(r_ij_edge_data.GetLength(), norm_2(r_i_vel), norm_2(r_j_vel));
                    const double c_tau = mDiffusionConstant * DeltaTime / local_dt;
                    const double Mc_ij = r_ij_edge_data.GetOffDiagonalConsistentMass();
                    const double AEC_ij = Mc_ij * (c_tau * (u_i - u_j) + (du_h_i - du_h_j)) / DeltaTime;

                    // Prelimiting step (avoids steepening)
                    // This step removes the antidiffusive fluxes that flatten the solution profile instead of steepening it
                    if (AEC_ij * (u_j - u_i) > 0.0) {
                        r_ij_edge_data.SetAntidiffusiveEdgeContribution(0.0);
                    } else {
                        r_ij_edge_data.SetAntidiffusiveEdgeContribution(AEC_ij);
                    }
                }
            }
        });
    }

    void EvaluateLimiter(const double DeltaTime)
    {
        // Allocate and initialize auxiliary limiter arrays
        std::vector<double> P_min(mAuxSize, 0.0);
        std::vector<double> P_max(mAuxSize, 0.0);
        std::vector<double> Q_min(mAuxSize, 0.0);
        std::vector<double> Q_max(mAuxSize, 0.0);

        // Get edge data structure containers
        const auto &r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto &r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Calculate the summation of positive/negative antidiffusive fluxes and bound fluxes
        IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // i-node nodal data
                const double u_i = mSolution[iRow]; //u_l
                // const double u_i = mSolutionOld[iRow]; //u_n
                const double M_l_i = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(iRow);

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double u_j = mSolution[j_node_id]; //u_l
                    // const double u_j = mSolutionOld[j_node_id]; //u_n
                    const double M_l_j = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(j_node_id);

                    // Get the antidiffusive flux edge contribution
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const double AEC_ij = r_ij_edge_data.GetAntidiffusiveEdgeContribution();

                    // Assemble the positive/negative antidiffusive fluxes
                    AtomicAdd(P_min[iRow], std::min(0.0, AEC_ij));
                    AtomicAdd(P_max[iRow], std::max(0.0, AEC_ij));
                    AtomicAdd(P_min[j_node_id], std::min(0.0, -AEC_ij));
                    AtomicAdd(P_max[j_node_id], std::max(0.0, -AEC_ij));

                    // Assemble the positive/negative antidiffusive bounds
                    const double aux_val_i = M_l_i * (u_j - u_i) / DeltaTime;
                    const double aux_val_j = M_l_j * (u_i - u_j) / DeltaTime;
                    Q_min[iRow] = std::min(Q_min[iRow], aux_val_i);
                    Q_max[iRow] = std::max(Q_max[iRow], aux_val_i);
                    Q_min[j_node_id] = std::min(Q_min[j_node_id], aux_val_j);
                    Q_max[j_node_id] = std::max(Q_max[j_node_id], aux_val_j);
                }
            }
        });

        // Calculate the nodal correction factors
        const double zero_tol = 1.0e-12;
        std::vector<double> R_min(mAuxSize);
        std::vector<double> R_max(mAuxSize);
        IndexPartition<IndexType>(mAuxSize).for_each([&](IndexType i){
            // R_{i}^{+} nodal correction factor calculation
            // Note that we check if the summation of positive antidiffusive fluxes is zero
            const double P_max_i = std::abs(P_max[i]) < zero_tol ? zero_tol : P_max[i];
            R_max[i] = std::min(1.0, Q_max[i] / P_max_i);

            // R_{i}^{-} nodal correction factor calculation
            // Note that we check if the summation of negative antidiffusive fluxes is zero
            const double P_min_i = std::abs(P_min[i]) < zero_tol ? -zero_tol : P_min[i];
            R_min[i] = std::min(1.0, Q_min[i] / P_min_i);
        });

        // Assembly of the final antidiffusive contribution
        std::vector<double> antidiff_assembly(mAuxSize, 0.0);
        IndexPartition<IndexType>(aux_n_rows).for_each([&](IndexType iRow){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iRow;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);

                    // Get the antidiffusive flux edge contribution
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const double AEC_ij = r_ij_edge_data.GetAntidiffusiveEdgeContribution();

                    // Set as correction factor the minimum in the edge
                    const double alpha_ij = AEC_ij > 0.0 ? std::min(R_max[iRow], R_min[j_node_id]) : std::min(R_min[iRow], R_max[j_node_id]);
                    KRATOS_ERROR_IF(alpha_ij < 0.0 || alpha_ij > 1.0) << "Correction factor is " << alpha_ij << " for edge " << iRow << "-" << j_node_id << std::endl;

                    // Assemble current edge antidiffusive contribution
                    const double aux = alpha_ij * AEC_ij;
                    AtomicAdd(antidiff_assembly[iRow], aux);
                    AtomicAdd(antidiff_assembly[j_node_id], -aux);
                }
            }
        });

        // Solve and add final antidiffusive contributions
        // Note that in here it is assumed that the low order solution is already added
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
            // Get nodal data
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();

            // Solve antidiffusive contribution
            const double M_l = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(i_node_id);
            const double solve_coeff = DeltaTime / M_l;
            mSolution[i_node_id] += solve_coeff * antidiff_assembly[i_node_id];
        });
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

    ///@}
}; // Class FluxCorrectedTransportConvectionProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// Input stream function
template<unsigned int TDim>
inline std::istream& operator >> (
    std::istream& rIStream,
    FluxCorrectedTransportConvectionProcess<TDim>& rThis);

/// Output stream function
template< unsigned int TDim>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const FluxCorrectedTransportConvectionProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
