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
        mMaxAllowedCFL = ThisParameters["max_delta_time"].GetDouble();
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
        const auto& r_row_indices = mpEdgeDataStructure->GetRowIndices();
        SizeType n_nodes = r_row_indices.size() - 1; // Note that last node is the NNZ
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
            const auto& r_node = mpModelPart->GetNode(iNode);
            mConvectedValuesOld[iNode] = r_node.GetSolutionStepValue(*mpConvectedVar,1); //solution at n
            mConvectionValues[iNode] = r_node.GetSolutionStepValue(*mpConvectionVar,1); //convective velocity at n
            mConvectionValuesOld[iNode] = r_node.GetSolutionStepValue(*mpConvectionVar,2); //convective velocity at n-1
        });

        // Substepping loop according to max CFL
        while (time - prev_time < 1.0e-12) {
            // Evaluate current substep maximum allowable delta time
            const double dt = this->CalculateSubStepDeltaTime(prev_time, time);

            // Solve current substep
            this->SolveSubStep(dt);

            // Advance in time
            prev_time += dt;

            // Update solution
            IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
                mConvectedValuesOld[iNode] = mConvectedValues[iNode];
            });
        }

        // Set final solution in the model part database
        IndexPartition<IndexType>(mpModelPart->NumberOfNodes()).for_each([this](IndexType iNode){
            (mpModelPart->GetNode(iNode)).FastGetSolutionStepValue(*mpConvectedVar) = mConvectedValues[iNode];
        });


        // // Set convection problem data
        // auto& r_conv_process_info = mpDistanceModelPart->GetProcessInfo();
        // const double previous_delta_time = r_conv_process_info.GetValue(DELTA_TIME);
        // const double dt =  previous_delta_time / static_cast<double>(n_substep);
        // r_conv_process_info.SetValue(DELTA_TIME, dt);
        // r_conv_process_info.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(*mpLevelSetVar);
        // r_conv_process_info.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetGradientVariable(*mpLevelSetGradientVar);
        // r_conv_process_info.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetConvectionVariable(*mpConvectVar);

        // // Save current level set value and current and previous step velocity values
        // // If the nodal stabilization tau is to be used, it is also computed in here
        // IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each([&](int i_node){
        //     const auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
        //     mVelocity[i_node] = it_node->FastGetSolutionStepValue(*mpConvectVar);
        //     mVelocityOld[i_node] = it_node->FastGetSolutionStepValue(*mpConvectVar,1);
        //     mOldDistance[i_node] = it_node->FastGetSolutionStepValue(*mpLevelSetVar,1);
        // });

        // for(unsigned int step = 1; step <= n_substep; ++step){
        //     KRATOS_INFO_IF("FluxCorrectedTransportConvectionProcess", mLevelSetConvectionSettings["echo_level"].GetInt() > 0) <<
        //         "Doing step "<< step << " of " << n_substep << std::endl;

        //     // Compute shape functions of old and new step
        //     const double Nold = 1.0 - static_cast<double>(step) / static_cast<double>(n_substep);
        //     const double Nnew = 1.0 - Nold;

        //     const double Nold_before = 1.0 - static_cast<double>(step-1) / static_cast<double>(n_substep);
        //     const double Nnew_before = 1.0 - Nold_before;

        //     // Emulate clone time step by copying the new distance onto the old one
        //     IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each(
        //     [&](int i_node){
        //         auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

        //         const array_1d<double,3>& r_v = mVelocity[i_node];
        //         const array_1d<double,3>& r_v_old = mVelocityOld[i_node];

        //         noalias(it_node->FastGetSolutionStepValue(*mpConvectVar)) = Nold * r_v_old + Nnew * r_v;
        //         noalias(it_node->FastGetSolutionStepValue(*mpConvectVar, 1)) = Nold_before * r_v_old + Nnew_before * r_v;
        //         it_node->FastGetSolutionStepValue(*mpLevelSetVar, 1) = it_node->FastGetSolutionStepValue(*mpLevelSetVar);
        //     }
        //     );

        //     mpSolvingStrategy->InitializeSolutionStep();
        //     mpSolvingStrategy->Predict();
        //     mpSolvingStrategy->SolveSolutionStep(); // forward convection to reach phi_n+1
        //     mpSolvingStrategy->FinalizeSolutionStep();

        // }

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
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
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
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
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
        const double max_dt = IndexPartition<IndexType>(n_nodes).for_each<MinReduction<double>>([&](IndexType iNode){
            KRATOS_WATCH(iNode)
            double dt_ij = std::numeric_limits<double>::max();
            const auto it_row = r_row_indices.begin() + iNode;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;
            KRATOS_WATCH(n_cols)
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
                    const double dx_CFL = ij_length*mMaxAllowedCFL;
                    const double dt_i = i_vel_norm > 1.0e-12 ? dx_CFL/i_vel_norm : mMaxAllowedDt;
                    const double dt_j = j_vel_norm > 1.0e-12 ? dx_CFL/j_vel_norm : mMaxAllowedDt;
                    dt_ij = std::min(dt_i, dt_j);

                    std::cout << "Dt in edge " << iNode << "-" << j_node_id << " : " << dt_ij << std::endl;
                }
            }

            return dt_ij;
        });

        // Check that current time step doesn't exceed the target time
        // Return maximum time step if it doesn't exceed the target time
        // Otherwise the difference between current time and target one is used (potentially in last substep)
        return PreviousTime + max_dt > TargetTime ? TargetTime - PreviousTime : max_dt;
    }

    void SolveSubStep(const double DeltaTime)
    {
        // Get edge data structure containers
        const auto& r_row_indices = mpEdgeDataStructure->GetRowIndices();
        const auto& r_col_indices = mpEdgeDataStructure->GetColIndices();
        const auto& r_values = mpEdgeDataStructure->GetValues();

        // Edge contributions assembly
        // Note that we don't loop the last entry of the row container as it is NNZ
        SizeType n_nodes = r_row_indices.size() - 1;
        IndexPartition<IndexType>(n_nodes).for_each([&](IndexType iNode){
            KRATOS_WATCH(iNode)
            const auto it_row = r_row_indices.begin() + iNode;
            const IndexType i_col_index = *it_row;
            const SizeType n_cols = *(it_row+1) - i_col_index;
            KRATOS_WATCH(n_cols)
            // Check that there are CSR columns (i.e. that current node involves an edge)
            if (n_cols != 0) {
                // i-node nodal data
                auto& r_i_node = mpModelPart->GetNode(iNode);
                auto& r_i_val = r_i_node.FastGetSolutionStepValue(*mpConvectedVar);
                auto& r_i_vel = r_i_node.FastGetSolutionStepValue(*mpConvectionVar);

                // j-node nodal loop (i.e. loop ij-edges)
                const auto i_col_begin = r_col_indices.begin() + i_col_index;
                for (IndexType j_node = 0; j_node < n_cols; ++j_node) {
                    // j-node nodal data
                    IndexType j_node_id = *(i_col_begin + j_node);
                    auto& r_j_node = mpModelPart->GetNode(j_node_id);
                    auto& r_j_val = r_j_node.FastGetSolutionStepValue(*mpConvectedVar);

                    // Get ij-edge operators data from CSR data structure
                    const auto& r_ij_edge_data = mpEdgeDataStructure->GetEdgeData(iNode, j_node_id);

                    // Atomic additions
                    double aux = 0.0;
                    AtomicAdd(r_i_val, aux);
                    AtomicAdd(r_j_val, -aux);

                    std::cout << "Assembling edge " << iNode << "-" << j_node_id << std::endl;
                }
            }
        });
    }

    unsigned int EvaluateNumberOfSubsteps()
    {
        // // First of all compute the maximum local CFL number
        // const double dt = mpDistanceModelPart->GetProcessInfo()[DELTA_TIME];
        // double max_cfl_found = block_for_each<MaxReduction<double>>(mpDistanceModelPart->Elements(), [&](Element& rElement){
        //     double vol;
        //     array_1d<double, TDim+1 > N;
        //     BoundedMatrix<double, TDim+1, TDim > DN_DX;
        //     auto& r_geom = rElement.GetGeometry();
        //     GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, vol);

        //     // Compute h
        //     double h=0.0;
        //     for(unsigned int i=0; i<TDim+1; i++){
        //         double h_inv = 0.0;
        //         for(unsigned int k=0; k<TDim; k++){
        //             h_inv += DN_DX(i,k)*DN_DX(i,k);
        //         }
        //         h += 1.0/h_inv;
        //     }
        //     h = sqrt(h)/static_cast<double>(TDim+1);

        //     // Get average velocity at the nodes
        //     array_1d<double, 3 > vgauss = ZeroVector(3);
        //     for(unsigned int i=0; i<TDim+1; i++){
        //         vgauss += N[i]* r_geom[i].FastGetSolutionStepValue(*mpConvectVar);
        //     }

        //     double cfl_local = norm_2(vgauss) / h;
        //     return cfl_local;
        // });
        // max_cfl_found *= dt;

        // // Synchronize maximum CFL between processes
        // max_cfl_found = mpDistanceModelPart->GetCommunicator().GetDataCommunicator().MaxAll(max_cfl_found);

        // unsigned int n_steps = static_cast<unsigned int>(max_cfl_found / mMaxAllowedCFL);
        // if(n_steps < 1){
        //     n_steps = 1;
        // }

		// // Now we compare with the maximum set
		// if (mMaxSubsteps > 0 && mMaxSubsteps < n_steps){
        //     n_steps = mMaxSubsteps;
        // }

        // return n_steps;
        return 0;
    }

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


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
    FluxCorrectedTransportConvectionProcess& operator=(FluxCorrectedTransportConvectionProcess const& rOther);

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
