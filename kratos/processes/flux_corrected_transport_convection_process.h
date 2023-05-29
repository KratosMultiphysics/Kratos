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
        const SizeType n_nodes = mpModelPart->NumberOfNodes();
        mConvectedValues.resize(n_nodes);
        mConvectionValues.resize(n_nodes);
        mConvectedValuesOld.resize(n_nodes);
        mConvectionValuesOld.resize(n_nodes);
        IndexPartition<IndexType>(n_nodes).for_each([this](IndexType iNode) {
            mConvectedValues[iNode] = 0.0;
            mConvectedValuesOld[iNode] = 0.0;
            mConvectionValues[iNode] = ZeroVector(3);
            mConvectionValuesOld[iNode] = ZeroVector(3);
        });
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
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
            const auto& r_node = mpModelPart->GetNode(iNode+1);
            mConvectedValuesOld[iNode] = r_node.GetSolutionStepValue(*mpConvectedVar,1); //solution at n
            mConvectionValues[iNode] = r_node.GetSolutionStepValue(*mpConvectionVar,1); //convective velocity at n
            mConvectionValuesOld[iNode] = r_node.GetSolutionStepValue(*mpConvectionVar,2); //convective velocity at n-1
        });

        // Substepping loop according to max CFL
        IndexType step = 1;
        while (time - prev_time > 1.0e-12) {
            // Evaluate current substep maximum allowable delta time
            const double dt = this->CalculateSubStepDeltaTime(prev_time, time);

            // Solve current substep
            KRATOS_INFO_IF("FluxCorrectedTransportConvectionProcess", mEchoLevel > 0) << "Substep " << step << " - \u0394t " << dt << std::endl;
            this->SolveSubStep(dt);

            // Update solution
            IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
                mConvectedValuesOld[iNode] = mConvectedValues[iNode];
            });

            // Advance in time
            step++;
            prev_time += dt;
        }

        // Set final solution in the model part database
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
            auto it_node = mpModelPart->NodesBegin() + iNode;
            it_node->FastGetSolutionStepValue(*mpConvectedVar) = mConvectedValues[iNode];
        });

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        mConvectedValues.clear();
        mConvectionValues.clear();
        mConvectedValuesOld.clear();
        mConvectionValuesOld.clear();
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

    ModelPart* mpModelPart;

	SizeType mEchoLevel;

    double mMaxAllowedDt;

    double mMaxAllowedCFL;

    bool mPerformInitialize = true;

    const Variable<double>* mpConvectedVar = nullptr;

    const Variable<array_1d<double,3>>* mpConvectionVar = nullptr;

    std::vector<double> mConvectedValues;

    std::vector<double> mConvectedValuesOld;

    std::vector<array_1d<double,3>> mConvectionValues;

    std::vector<array_1d<double,3>> mConvectionValuesOld;

    typename EdgeBasedDataStructure<TDim>::UniquePointer mpEdgeDataStructure;

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
                    const double j_vel_norm = norm_2(mConvectionValues[j_node_id-1]);

                    // Get ij-edge length from CSR data structure
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iNode, j_node_id);
                    const double ij_length = r_ij_edge_data.GetLength();

                    // Calculate the max time step from the velocities at both edge ends and get the minimum
                    const double dx_CFL = ij_length*mMaxAllowedCFL;
                    const double dt_i = i_vel_norm > 1.0e-12 ? dx_CFL/i_vel_norm : mMaxAllowedDt;
                    const double dt_j = j_vel_norm > 1.0e-12 ? dx_CFL/j_vel_norm : mMaxAllowedDt;
                    dt_ij = std::min(dt_i, dt_j);
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

    void SolveSubStep(const double DeltaTime)
    {
        // Get edge data structure containers
        const auto& r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto& r_col_indices = mpEdgeDataStructure->GetColIndices();

        // Aux TLS container
        struct AuxTLS
        {
            array_1d<double,TDim> F_i;
            array_1d<double,TDim> F_j;
        };

        // Edge contributions assembly
        // Note that we don't loop the last entry of the row container as it is NNZ
        SizeType n_nodes = r_row_indices.size() - 1;
        IndexPartition<IndexType>(n_nodes).for_each(AuxTLS(), [&](IndexType iNode, AuxTLS& rTLS){
            // Get current row (node) storage data
            const auto it_row = r_row_indices.begin() + iNode;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;

            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // Get TLS data
                auto& F_i = rTLS.F_i;
                auto& F_j = rTLS.F_j;

                // i-node nodal data
                double& r_u_i = mConvectedValuesOld[iNode];
                const auto& r_i_vel = mConvectionValues[iNode];

                // i-node convective flux calculation
                for (IndexType d = 0; d < TDim; ++d) {
                    F_i[d] = r_u_i * r_i_vel[d];
                }

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    double& r_u_j = mConvectedValuesOld[j_node_id-1];
                    const auto& r_j_vel = mConvectionValues[j_node_id-1];

                    // j-node convective flux calculation
                    for (IndexType d = 0; d < TDim; ++d) {
                        F_j[d] = r_u_j * r_j_vel[d];
                    }

                    // Get ij-edge operators data from CSR data structure
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iNode, j_node_id);
                    const auto& r_d_ij = r_ij_edge_data.GetFirstDerivatives();

                    // Calculate fluxes along edges
                    double f_i = 0.0;
                    double f_j = 0.0;
                    double D_ij = 0.0;
                    for (IndexType d = 0; d < TDim; ++d) {
                        D_ij += std::pow(r_d_ij[d],2);
                        f_i += r_d_ij[d] * F_i[d];
                        f_j -= r_d_ij[d] * F_j[d];
                    }
                    D_ij = std::sqrt(D_ij);
                    f_i /= D_ij;
                    f_j /= D_ij;

                    // Calculate residual from numerical flux
                    double F_ij_d;
                    double r_i = 0.0;
                    double r_j = 0.0;
                    const double u_ij = 0.5*(r_u_i + r_u_j) - 0.5*DeltaTime*(f_i-f_j)/D_ij; // u_{ij}^{n+0.5}
                    for (IndexType d = 0; d < TDim; ++d) {
                        //TODO: We're using the old values. Probably we should use an estimate of the substep midpoint ones
                        F_ij_d = (r_i_vel[d] + r_j_vel[d]) * u_ij; // F(u_{ij}^{n+0.5})
                        r_i -= 0.5*r_d_ij[d]*F_ij_d;
                        r_j += 0.5*r_d_ij[d]*F_ij_d;
                    }

                    // Calculate nodal explicit low order contributions
                    const double diff_coeff = 1.0;
                    const double M_c = r_ij_edge_data.GetMassCoefficient();
                    const double M_l = r_ij_edge_data.GetLumpedMassCoefficient();
                    double delta_u_i_low = (DeltaTime / M_l) * (r_i + diff_coeff * (M_c * r_u_i + M_c * r_u_j - M_l * r_u_i));
                    double delta_u_j_low = (DeltaTime / M_l) * (r_j + diff_coeff * (M_c * r_u_i + M_c * r_u_j - M_l * r_u_j));


                    // Calculate antidiffusive fluxes (high order contribution)
                    //TODO: Implement this!

                    // Atomic additions
                    //FIXME: Remove after testing
                    double delta_u_i = delta_u_i_low;
                    double delta_u_j = delta_u_j_low;
                    AtomicAdd(mConvectedValues[iNode], delta_u_i);
                    AtomicAdd(mConvectedValues[j_node_id - 1], delta_u_j);
                }
            }
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
