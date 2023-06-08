//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "containers/edge_based_data_structure.h"
#include "includes/define.h"
#include "includes/gid_io.h"
#include "processes/process.h"
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

/// Short class definition.
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and convects a level set distance
* on the top of it
*/
template<unsigned int TDim>
class FluxCorrectedTransportConvectionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using SizeType = std::size_t;

    using IndexType = std::size_t;

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
        mMaxAllowedCFL = ThisParameters["max_CFL"].GetDouble();
        mMaxAllowedDt = ThisParameters["max_delta_time"].GetDouble();
        mpConvectedVar = &KratosComponents<Variable<double>>::Get(ThisParameters["convected_variable_name"].GetString());
        mpConvectionVar = &KratosComponents<Variable<array_1d<double, 3>>>::Get(ThisParameters["convection_variable_name"].GetString());
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
        // Fill the edge-based data structure
        mpEdgeDataStructure->CalculateEdgeDataStructure(*mpModelPart);

        // Auxiliary flag to avoid resetting the edge data structure at each Execute call
        mPerformInitialize = false;

        // Allocate auxiliary arrays
        // Note that we give the size such that the index matches the id of the corresponding node
        // Also note that this implies having a "fake" 0-id entry in the first position of the arrays
        mAuxSize = mpModelPart->NumberOfNodes() + 1;
        mResidual.resize(mAuxSize);
        mUpdateLowOrder.resize(mAuxSize);
        mUpdateHighOrder.resize(mAuxSize);
        mSolution.resize(mAuxSize);
        mConvectionValues.resize(mAuxSize);
        mConvectionValuesOld.resize(mAuxSize);
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
        const double time = r_process_info[TIME];
        double prev_time = r_prev_process_info[TIME];

        // Update data vectors
        const SizeType n_nodes = mpModelPart->NumberOfNodes();
        IndexPartition<IndexType>(n_nodes).for_each([this](IndexType iNode){
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();
            mSolution[i_node_id] = it_node->GetSolutionStepValue(*mpConvectedVar,1); //solution at n
            mConvectionValues[i_node_id] = it_node->GetSolutionStepValue(*mpConvectionVar,1); //convective velocity at n
            mConvectionValuesOld[i_node_id] = it_node->GetSolutionStepValue(*mpConvectionVar,2); //convective velocity at n-1
        });

        GidIO<> gid_io_convection("/home/rzorrilla/Desktop/FluxCorrectedTransportConvectionProcess2D", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        gid_io_convection.InitializeMesh(0);
        gid_io_convection.WriteMesh(mpModelPart->GetMesh());
        gid_io_convection.FinalizeMesh();
        gid_io_convection.InitializeResults(0, mpModelPart->GetMesh());
        gid_io_convection.WriteNodalResults(DISTANCE, mpModelPart->Nodes(), 0, 0);
        gid_io_convection.WriteNodalResults(VELOCITY, mpModelPart->Nodes(), 0, 0);

        // Substepping loop according to max CFL
        IndexType step = 1;
        while (time - prev_time > 1.0e-12) {
            // Evaluate current substep maximum allowable delta time
            const double dt = this->CalculateSubStepDeltaTime(prev_time, time);

            // Solve current substep
            KRATOS_INFO_IF("FluxCorrectedTransportConvectionProcess", mEchoLevel > 0) << "Substep " << step << " - \u0394t " << dt << std::endl;
            this->SolveSubStep(dt);

            // TODO:Remove after debugging
            IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
                auto it_node = mpModelPart->NodesBegin() + iNode;
                it_node->FastGetSolutionStepValue(*mpConvectedVar) = mSolution[it_node->Id()];
            });
            gid_io_convection.WriteNodalResults(DISTANCE, mpModelPart->Nodes(), step, 0);
            gid_io_convection.WriteNodalResults(VELOCITY, mpModelPart->Nodes(), step, 0);

            // Advance in time
            step++;
            prev_time += dt;
        }

        gid_io_convection.FinalizeResults();

        // Set final solution in the model part database
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
            auto it_node = mpModelPart->NodesBegin() + iNode;
            it_node->FastGetSolutionStepValue(*mpConvectedVar) = mSolution[it_node->Id()];
        });

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        mResidual.clear();
        mUpdateLowOrder.clear();
        mUpdateHighOrder.clear();
        mSolution.clear();
        mConvectionValues.clear();
        mConvectionValuesOld.clear();
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
            "max_delta_time" : 1.0
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

    SizeType mAuxSize; /// Size of the internal vector containers that equals the number of nodes plus one

    SizeType mEchoLevel; /// Level of information that is output

    double mMaxAllowedDt; /// Maximum allowed time step for the substepping

    double mMaxAllowedCFL; /// Maximum allowed CFL number for the substepping

    bool mPerformInitialize = true; /// Flag to indicate if the ExecuteInitialize is required to be done

    const Variable<double>* mpConvectedVar = nullptr; /// Pointer to the convected variable (e.g. DISTANCE)

    const Variable<array_1d<double,3>>* mpConvectionVar = nullptr; /// Pointer to the convection variable (e.g. VELOCITY)

    std::vector<double> mResidual; /// Auxiliary vector to assemble the edge residuals

    std::vector<double> mSolution; /// Auxiliary vector containing the solution

    std::vector<double> mUpdateLowOrder;

    std::vector<double> mUpdateHighOrder;

    std::vector<array_1d<double,3>> mConvectionValues; /// Auxiliary vector to store the convection variable values (e.g. VELOCITY) 

    std::vector<array_1d<double,3>> mConvectionValuesOld; /// Auxiliary vector to store the previous step convection variable values (e.g. VELOCITY)

    typename EdgeBasedDataStructure<TDim>::UniquePointer mpEdgeDataStructure; /// Pointer to the edge-based data structure

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    double CalculateSubStepDeltaTime(
        const double PreviousTime,
        const double TargetTime)
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
                    // const double dx_CFL = ij_length*mMaxAllowedCFL;
                    // const double dt_i = i_vel_norm > 1.0e-12 ? dx_CFL/i_vel_norm : mMaxAllowedDt;
                    // const double dt_j = j_vel_norm > 1.0e-12 ? dx_CFL/j_vel_norm : mMaxAllowedDt;
                    // dt_ij = std::min(dt_i, dt_j);
                    dt_ij = CalculateEdgeLocalDeltaTime(ij_length, i_vel_norm, j_vel_norm);
                }
            }

            return dt_ij;
        });

        // Check that current time step doesn't exceed the target time
        // Use maximum allowable time step if it doesn't exceed the target time
        // Otherwise the difference between current time and target one is used (potentially in last substep)
        double dt = PreviousTime + min_dt > TargetTime ? TargetTime - PreviousTime : min_dt;

        // Synchronize among MPI nodes
        mpModelPart->GetCommunicator().GetDataCommunicator().MinAll(dt);

        // Return maximum allowable time step
        return dt;
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
        // Get edge data structure containers
        const auto& r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto& r_col_indices = mpEdgeDataStructure->GetColIndices();
        SizeType aux_n_rows = r_row_indices.size() - 1; // Note that the last entry of the row container is the NNZ

        // Initialize containers for the edge assembly
        IndexPartition<IndexType>(mAuxSize).for_each([this](IndexType i){
            mResidual[i] = 0.0;
            mUpdateLowOrder[i] = 0.0;
            mUpdateHighOrder[i] = 0.0;
        });

        // Aux TLS container
        struct AuxTLS
        {
            array_1d<double,TDim> F_i;
            array_1d<double,TDim> F_j;
            array_1d<double,TDim> d_ij;
            array_1d<double,TDim> b_ij;
            array_1d<double,TDim> b_node;
        };

        // Edge contributions assembly
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
                auto& b_node = rTLS.b_node;

                // i-node nodal data
                const double u_i = mSolution[iRow];
                const auto& r_i_vel = mConvectionValues[iRow]; //TODO: Use substep velocity

                // // i-node mass matrices diagonal contributions
                // const double M_c_i = mpEdgeDataStructure->GetMassMatrixDiagonal(iRow);
                // const double M_l_i = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(iRow);

                // i-node convective flux calculation
                for (IndexType d = 0; d < TDim; ++d) {
                    F_i[d] = u_i * r_i_vel[d];
                }

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    const double u_j = mSolution[j_node_id];
                    const auto& r_j_vel = mConvectionValues[j_node_id]; //TODO: Use substep velocity

                    // // j-node mass matrices diagonal contributions
                    // const double M_c_j = mpEdgeDataStructure->GetMassMatrixDiagonal(j_node_id);
                    // const double M_l_j = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(j_node_id);

                    // // Laplacian problem
                    // // TODO: Remove once we check that the explicit update is correct
                    // const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    // const double c_edge = r_ij_edge_data.GetOffDiagonalLaplacian();
                    // const double res_edge_i = 0.1 * c_edge * (u_i - u_j);
                    // const double res_edge_j = 0.1 * c_edge * (u_j - u_i);

                    // j-node convective flux calculation
                    for (IndexType d = 0; d < TDim; ++d) {
                        F_j[d] = u_j * r_j_vel[d];
                    }

                    // Get ij-edge operators data from CSR data structure
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iRow, j_node_id);
                    const auto& r_Ni_DNj = r_ij_edge_data.GetOffDiagonalConvective();
                    const auto& r_DNi_Nj = r_ij_edge_data.GetOffDiagonalConvectiveTranspose();

                    // // Calculate fluxes along edges
                    // d_ij = 0.5 * (r_DNi_Nj - r_Ni_DNj);
                    // double f_i = 0.0;
                    // double f_j = 0.0;
                    // double D_ij = 0.0;
                    // for (IndexType d = 0; d < TDim; ++d) {
                    //     D_ij += std::pow(d_ij[d],2);
                    //     f_i += d_ij[d] * F_i[d];
                    //     f_j -= d_ij[d] * F_j[d];
                    // }
                    // D_ij = std::sqrt(D_ij);
                    // f_i /= D_ij;
                    // f_j /= D_ij;

                    // // Calculate numerical flux from the fluxes along edges
                    // // Note that this numerical flux corresponds to a Lax-Wendroff scheme
                    // const double F_ij = (f_i + f_j) - DeltaTime * (f_i - f_j) / D_ij;

                    // // Calculate convection volume residual
                    // double res_edge_i = -D_ij * F_ij;
                    // double res_edge_j = D_ij * F_ij;

                    // Standard convection to test
                    d_ij = 0.5 * (r_Ni_DNj - r_DNi_Nj);
                    double res_edge_i = inner_prod(-d_ij, F_i + F_j);
                    double res_edge_j = inner_prod(d_ij, F_i + F_j);

                    // If current edge belogs to a boundary, add the corresponding convection boundary integrals
                    if (r_ij_edge_data.IsBoundary()) {
                        // Get ij-edge boundary operators from CSR data structure
                        const auto &r_N_N_normal = r_ij_edge_data.GetConvectiveBoundary();
                        b_node = r_N_N_normal;
                        b_ij = 0.5 * r_N_N_normal;

                        // Add boundary contribution to the residual
                        const double res_edge_bd = inner_prod(b_ij, F_i + F_j);
                        res_edge_i += -res_edge_bd - inner_prod(b_node, F_i);
                        res_edge_j += +res_edge_bd + inner_prod(b_node, F_j);
                    }

                    // // Add Laplacian term
                    // // TODO: Remove once we check that the explicit update is correct
                    // const double c_edge = r_ij_edge_data.GetOffDiagonalLaplacian();
                    // res_edge_i += 0.01 * c_edge * (u_i - u_j);
                    // res_edge_j += 0.01 * c_edge * (u_j - u_i);

                    // // Add low order scheme diffusion (Kuzmin)
                    // double k_ij = 0.0;
                    // double k_ji = 0.0;
                    // for (IndexType d = 0; d < TDim; ++d) {
                    //     k_ij += r_Ni_DNj[d] * r_j_vel[d];
                    //     k_ji += r_DNi_Nj[d] * r_i_vel[d];
                    // }
                    // const double d_ij_diff = std::max({-k_ij, 0.0, -k_ji});
                    // res_edge_i += d_ij_diff * (u_j - u_i);
                    // res_edge_j += d_ij_diff * (u_i - u_j);

                    // Add low order scheme diffusion (Rainald)
                    const double local_dt = CalculateEdgeLocalDeltaTime(r_ij_edge_data.GetLength(), norm_2(r_i_vel), norm_2(r_j_vel));
                    const double c_tau = DeltaTime / local_dt;
                    const double Mc_i_j = r_ij_edge_data.GetOffDiagonalConsistentMass();
                    const double ij_low_order_diff =  c_tau * Mc_i_j;
                    res_edge_i += ij_low_order_diff * (u_j - u_i);
                    res_edge_j += ij_low_order_diff * (u_i - u_j);

                    // Calculate antidiffusive fluxes (high order contribution)
                    //TODO: Implement this!

                    // Atomic additions
                    AtomicAdd(mResidual[iRow], res_edge_i);
                    AtomicAdd(mResidual[j_node_id], res_edge_j);
                }
            }
        });

        // Do the explicit lumped mass matrix solve and apply the correction to current solution
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([&](IndexType iNode){
            const auto it_node = mpModelPart->NodesBegin() + iNode;
            const IndexType i_node_id = it_node->Id();
            const double M_l = mpEdgeDataStructure->GetLumpedMassMatrixDiagonal(i_node_id);
            const double solve_coeff = DeltaTime / M_l;
            mSolution[i_node_id] += solve_coeff * mResidual[i_node_id];
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
