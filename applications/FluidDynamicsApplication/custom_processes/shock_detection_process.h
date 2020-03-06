//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#ifndef KRATOS_SHOCK_DETECTION_PROCESS
#define KRATOS_SHOCK_DETECTION_PROCESS

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "processes/compute_nodal_gradient_process.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Main class for shock detection
/** This class implements some utilities for the detection of sharp discontinuitites (shocks) in the FE solution
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockDetectionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShockDetectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShockDetectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ShockDetectionProcess() = default;

    /// Constructor with default shock sensor variable
    ShockDetectionProcess(ModelPart& rModelPart)
    : Process()
    , mrModelPart(rModelPart)
    , mrShockSensorVariable(SHOCK_SENSOR)
    {}

    /// Constructor with custom shock sensor variable
    ShockDetectionProcess(
        ModelPart& rModelPart,
        Variable<double>& rShockSensorVariable)
    : Process()
    , mrModelPart(rModelPart)
    , mrShockSensorVariable(rShockSensorVariable)
    {}

    /// Destructor.
    virtual ~ShockDetectionProcess() = default;

    /// Assignment operator.
    ShockDetectionProcess &operator=(ShockDetectionProcess const &rOther) = delete;

    /// Copy constructor.
    ShockDetectionProcess(ShockDetectionProcess const &rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override
    {
        // Calculate the nodal area
        CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>(
            mrModelPart,
            mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE)).Execute();

        // Calculate the nodal neighbours for the edge-based shock capturing
        const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
        FindGlobalNodalNeighboursProcess(r_data_communicator, mrModelPart).Execute();
    }

    void EdgeBasedShockDetection(
        const Variable<double>& rShockVariable,
        const Variable<array_1d<double, 3>>& rShockGradientVariable)
    {
        // If required recompute the NODAL_AREA
        if (mUpdateNodalAreaAtEachStep) {
            CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>(
            mrModelPart,
            mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE)).Execute();
        }

        // If required recompute the NODAL_NEIGHBOURS
        if (mUpdateNodalNeighboursAtEachStep) {
            const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
            FindGlobalNodalNeighboursProcess(r_data_communicator, mrModelPart).Execute();
        }

        // Specialize the edge based shock detection
        EdgeBasedShockDetectionSpecialization<>(rShockVariable, rShockGradientVariable);
    }

    template<class TShockVariableType, class TShockGradientVariableType>
    void EdgeBasedShockDetectionSpecialization(
        const TShockVariableType& rShockVariable,
        const TShockGradientVariableType& rShockGradientVariable)
    {
        // Calculate the shock variable nodal gradients
        ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(
            mrModelPart,
            rShockVariable,
            rShockGradientVariable).Execute();

//         // TODO: THIS HAS TO BE UPDATED TO CONSIDER VECTOR VARIABLES
//         // Calculate the shock variable edge based nodal gradients
// #pragma omp parallel for
//         for (int i_node = 0; i_node < static_cast<int>(mrModelPart.NumberOfNodes()); ++i_node) {
//             auto it_node = mrModelPart.NodesBegin() + i_node;
//             const auto& r_var_i = it_node->FastGetSolutionStepValue(rShockVariable);
//             auto& r_grad_var_i = it_node->GetValue(rShockGradientVariable);

//             // Loop the neighbours to compute the nodal gradients
//             auto& r_neighbours = it_node->GetValue(NEIGHBOUR_NODES);
//             KRATOS_DEBUG_ERROR_IF(r_neighbours.size() == 0) << "Node " << i_node << " has no neighbours." << std::endl;

//             double sum_norm_l_ji = 0.0;
//             for (auto& r_neigh : r_neighbours ) {
//                 const auto& r_var_j = r_neigh.FastGetSolutionStepValue(rShockVariable);
//                 const auto l_ji = r_neigh.Coordinates() - it_node->Coordinates();
//                 const double norm_l_ji = norm_2(l_ji);
//                 const auto unit_l_ji = l_ji / norm_l_ji;
//                 sum_norm_l_ji += norm_l_ji;
//                 r_grad_var_i += (r_var_j - r_var_i) * unit_l_ji;
//             }

//             r_grad_var_i /= sum_norm_l_ji;
//         }

        // Perform the shock detection
#pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mrModelPart.NumberOfNodes()); ++i_node) {
            auto it_node = mrModelPart.NodesBegin() + i_node;
            double& r_shock_sens = it_node->GetValue(mrShockSensorVariable);
            const auto& r_var_i = it_node->FastGetSolutionStepValue(rShockVariable);
            const auto& r_grad_var_i = it_node->GetValue(rShockGradientVariable);

            // Loop the neighbours to compute the shock sensor
            r_shock_sens = 0.0;
            const double zero_tol = 1.0e-8;
            auto& r_neighbours = it_node->GetValue(NEIGHBOUR_NODES);
            KRATOS_DEBUG_ERROR_IF(r_neighbours.size() == 0) << "Node " << i_node << " has no neighbours." << std::endl;
            for (auto& r_neigh : r_neighbours ) {
                // Get the neighbour values
                const double& r_var_j = r_neigh.FastGetSolutionStepValue(rShockVariable);
                const auto& r_grad_var_j = r_neigh.GetValue(rShockGradientVariable);
                const auto l_ji = r_neigh.Coordinates() - it_node->Coordinates();

                // Calculate the density sensor auxiliary values
                // TODO THIS HAS TO BE PARTICULARIZED FOR ARRAY TYPES
                const auto aux_1 = r_var_j - r_var_i;
                const auto aux_2 = 0.5 * inner_prod(l_ji, r_grad_var_i + r_grad_var_j);
                const auto num = aux_1 - aux_2;
                const auto den = std::abs(aux_1) + std::abs(aux_2);
                double beta_ij = 0.0;

                if (it_node->Id() == 31626) {
                    KRATOS_WATCH(r_neigh.Id())
                    KRATOS_WATCH(num)
                    KRATOS_WATCH(den)
                }

                // Check if the solution is not constant (den close to 0.0)
                if (std::abs(den) > zero_tol) {
                    // Check if the density sensor is between the physical bounds (0.0 < var_sens < 1.0)
                    const double aux_beta_ij = std::abs(num / den);
                    if (aux_beta_ij < 1.0) {
                        beta_ij = aux_beta_ij;
                    } else {
                        beta_ij = 1.0;
                    }
                } else {
                    beta_ij = 0.0;
                }
                if (it_node->Id() == 31626) {
                    KRATOS_WATCH(beta_ij)
                }

                // Check against the current shock density sensor value and keep the largest one
                if (r_shock_sens < beta_ij) {
                    r_shock_sens = beta_ij;
                }
            }
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ShockDetectionProcess";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "ShockDetectionProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {}

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// Updates the NODAL_AREA at each time step (required in case the mesh deforms)
    const bool mUpdateNodalAreaAtEachStep = false;

    /// Updates the NODAL_NEIGHBOURS at each time step (required in case topology changes)
    const bool mUpdateNodalNeighboursAtEachStep = false;

    /// Reference to the model part in where the shock detection is to be performed
    ModelPart& mrModelPart;

    /// Reference to the shock sensor variable
    const Variable<double>& mrShockSensorVariable;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override{}

    void load(Serializer& rSerializer) override {}

    ///@}
}; // Class ShockDetectionProcess

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(
    std::istream &rIStream,
    ShockDetectionProcess &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const ShockDetectionProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_SHOCK_DETECTION_PROCESS  defined
