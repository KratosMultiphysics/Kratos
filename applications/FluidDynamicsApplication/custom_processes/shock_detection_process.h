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
#include "processes/calculate_nodal_area_process.h"
#include "processes/compute_nodal_gradient_process.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "containers/global_pointers_vector.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Auxiliary neighbour data class
 * Auxiliary class to retrieve the data from the neighbours when communicating the pointers
 * @tparam TShockVariableType Shock variable type
 * @tparam TShockGradientVariableType Shock gradient variable type
 */
template<class TShockVariableType, class TShockGradientVariableType>
class NeighbourData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NeighbourData
    KRATOS_CLASS_POINTER_DEFINITION(NeighbourData);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Neighbour Data object
     * Default neighbour data container constructor
     * Required to compile the class
     */
    NeighbourData() = default;

    /**
     * @brief Construct a new Neighbour Data object
     * Constructes a new neighbour data container instance
     * @param rShockVariableValue Neighbour shock variable value
     * @param rShockGradientVariableValue Neighbour shock variable gradient value
     * @param rCoordinates Neighbour node coordinates
     */
    NeighbourData(
        const typename TShockVariableType::Type& rShockVariableValue,
        const typename TShockGradientVariableType::Type& rShockGradientVariableValue,
        const array_1d<double, 3>& rCoordinates)
    {
        mCoordinates = rCoordinates;
        mShockVariableValue = rShockVariableValue;
        mShockGradientVariableValue = rShockGradientVariableValue;
    }

    ///@}
    ///@name Member Variables
    ///@{

    array_1d<double, 3> mCoordinates;
    typename TShockVariableType::Type mShockVariableValue;
    typename TShockGradientVariableType::Type mShockGradientVariableValue;

    ///@}
private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("mCoordinates",mCoordinates);
        rSerializer.save("mShockVariableValue",mShockVariableValue);
        rSerializer.save("mShockGradientVariableValue",mShockGradientVariableValue);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("mCoordinates",mCoordinates);
        rSerializer.load("mShockVariableValue",mShockVariableValue);
        rSerializer.load("mShockGradientVariableValue",mShockGradientVariableValue);
    }

    ///@}
};

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

    /// Variable component type
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > VariableComponentType;

    /// Node pointer type
    typedef typename Node<3>::Pointer NodePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ShockDetectionProcess() = default;

    /// Constructor with default shock sensor variable for double shock variable
    ShockDetectionProcess(
        ModelPart& rModelPart,
        const Variable<double>& rShockDoubleVariable,
        const Variable<array_1d<double,3>>& rShockGradientVariable,
        const bool UpdateNodalAreaAtEachStep = false,
        const bool UpdateNodalNeighboursAtEachStep = false)
    : Process()
    , mrModelPart(rModelPart)
    , mUpdateNodalAreaAtEachStep(UpdateNodalAreaAtEachStep)
    , mUpdateNodalNeighboursAtEachStep(UpdateNodalNeighboursAtEachStep)
    , mShockVariableIsDouble(true)
    , mpShockDoubleVariable(&rShockDoubleVariable)
    , mpShockGradientVariable(&rShockGradientVariable)
    , mpShockSensorVariable(&SHOCK_SENSOR)
    {}

    /// Constructor with default shock sensor variable for component shock variable
    ShockDetectionProcess(
        ModelPart& rModelPart,
        const VariableComponentType& rShockComponentVariable,
        const Variable<array_1d<double,3>>& rShockGradientVariable,
        const bool UpdateNodalAreaAtEachStep = false,
        const bool UpdateNodalNeighboursAtEachStep = false)
    : Process()
    , mrModelPart(rModelPart)
    , mUpdateNodalAreaAtEachStep(UpdateNodalAreaAtEachStep)
    , mUpdateNodalNeighboursAtEachStep(UpdateNodalNeighboursAtEachStep)
    , mShockVariableIsDouble(false)
    , mpShockComponentVariable(&rShockComponentVariable)
    , mpShockGradientVariable(&rShockGradientVariable)
    , mpShockSensorVariable(&SHOCK_SENSOR)
    {}

    /// Constructor with custom shock sensor variable for double shock variable
    ShockDetectionProcess(
        ModelPart& rModelPart,
        const Variable<double>& rShockDoubleVariable,
        const Variable<array_1d<double,3>>& rShockGradientVariable,
        const Variable<double>& rShockSensorVariable,
        const bool UpdateNodalAreaAtEachStep = false,
        const bool UpdateNodalNeighboursAtEachStep = false)
    : Process()
    , mrModelPart(rModelPart)
    , mUpdateNodalAreaAtEachStep(UpdateNodalAreaAtEachStep)
    , mUpdateNodalNeighboursAtEachStep(UpdateNodalNeighboursAtEachStep)
    , mShockVariableIsDouble(true)
    , mpShockDoubleVariable(&rShockDoubleVariable)
    , mpShockGradientVariable(&rShockGradientVariable)
    , mpShockSensorVariable(&rShockSensorVariable)
    {}

    /// Constructor with custom shock sensor variable for component shock variable
    ShockDetectionProcess(
        ModelPart& rModelPart,
        const VariableComponentType& rShockComponentVariable,
        const Variable<array_1d<double,3>>& rShockGradientVariable,
        const Variable<double>& rShockSensorVariable,
        const bool UpdateNodalAreaAtEachStep = false,
        const bool UpdateNodalNeighboursAtEachStep = false)
    : Process()
    , mrModelPart(rModelPart)
    , mUpdateNodalAreaAtEachStep(UpdateNodalAreaAtEachStep)
    , mUpdateNodalNeighboursAtEachStep(UpdateNodalNeighboursAtEachStep)
    , mShockVariableIsDouble(false)
    , mpShockComponentVariable(&rShockComponentVariable)
    , mpShockGradientVariable(&rShockGradientVariable)
    , mpShockSensorVariable(&rShockSensorVariable)
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

    /**
     * @brief Initializes the values for the shock detection
     * This method initializes the nodal mass, that is required for the nodal gradients
     * calculation, and the nodal neighbours.
     * It has to be executed once (in case there is no mesh deformation nor topology changes)
     */
    void ExecuteInitialize() override;

    /**
     * @brief Calculates the edge based shock detection
     * This method performs the edge based shock detection
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method performs all the operations
     * This method perform all the operations that are required for the shock detection
     */
    void Execute() override;

    /**
     * @brief Perform edge based shock detection
     * This method performs the edge based shock detection
     * @param rShockVariable Double variable to perform the shock detection
     * @param rShockGradientVariable Vector variable to calculate the shock variable gradients
     */
    void EdgeBasedShockDetection(
        const Variable<double>& rShockVariable,
        const Variable<array_1d<double, 3>>& rShockGradientVariable);

    /**
     * @brief Perform edge based shock detection
     * This method performs the edge based shock detection
     * @param rShockVariable Component variable to perform the shock detection
     * @param rShockGradientVariable Vector variable to calculate the shock variable gradients
     */
    void EdgeBasedShockDetection(
        const VariableComponentType& rShockVariable,
        const Variable<array_1d<double, 3>>& rShockGradientVariable);

    /**
     * @brief Template specialization of the edge based shock detection function
     * Auxiliary method to specialize the variable types for the edge based shock detection
     * @tparam TShockVariableType Shock variable type
     * @tparam TShockGradientVariableType Shock gradient variable type
     * @param rShockVariable Component variable to perform the shock detection
     * @param rShockGradientVariable Vector variable to calculate the shock variable gradients
     */
    template<class TShockVariableType, class TShockGradientVariableType>
    void EdgeBasedShockDetectionSpecialization(
        const TShockVariableType& rShockVariable,
        const TShockGradientVariableType& rShockGradientVariable)
    {
        // If required recompute the NODAL_AREA
        // This is required for the nodal gradients calculation
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

        // Calculate the shock variable nodal gradients
        ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(
            mrModelPart,
            rShockVariable,
            rShockGradientVariable).Execute();

        auto& r_comm = mrModelPart.GetCommunicator();
        auto& r_data_comm = r_comm.GetDataCommunicator();

        // Create the global pointers list
        GlobalPointersVector<Node<3>> global_pointers_list;
        if (r_comm.IsDistributed()) {
            for (auto &r_node : r_comm.LocalMesh().Nodes()) {
                auto& r_gp_to_neighbours = r_node.GetValue(NEIGHBOUR_NODES).GetContainer();
                for (auto &r_gp : r_gp_to_neighbours) {
                    global_pointers_list.push_back(r_gp);
                }
            }
            global_pointers_list.Unique();
        }

        // Now create the pointer communicator and shock values retrieve proxy
        GlobalPointerCommunicator<Node<3>> pointer_communicator(r_data_comm, global_pointers_list);
        auto shock_variables_proxy = pointer_communicator.Apply([&](const GlobalPointer<Node<3>>& rpNode)
        {
            NeighbourData<TShockVariableType, TShockGradientVariableType> neighbour_data(
                rpNode->FastGetSolutionStepValue(rShockVariable),
                rpNode->GetValue(rShockGradientVariable),
                rpNode->Coordinates());
            return neighbour_data;
        });

        // Perform the shock detection
#pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(r_comm.LocalMesh().NumberOfNodes()); ++i_node) {
            auto it_node = r_comm.LocalMesh().NodesBegin() + i_node;
            double& r_shock_sens = it_node->GetValue(*mpShockSensorVariable);
            const auto& r_var_i = it_node->FastGetSolutionStepValue(rShockVariable);
            const auto& r_grad_var_i = it_node->GetValue(rShockGradientVariable);

            // Loop the neighbours to compute the shock sensor
            r_shock_sens = 0.0;
            const double zero_tol = 1.0e-8;
            auto& r_neighbours = it_node->GetValue(NEIGHBOUR_NODES);
            KRATOS_DEBUG_ERROR_IF(r_neighbours.size() == 0) << "Node " << i_node << " has no neighbours." << std::endl;
            for (auto& r_neigh : r_neighbours ) {
                // Get the neighbour values
                const auto values_j = shock_variables_proxy.Get(&r_neigh);
                const double& r_var_j = values_j.mShockVariableValue;
                const auto& r_grad_var_j = values_j.mShockGradientVariableValue;
                const auto l_ji = values_j.mCoordinates - it_node->Coordinates();

                // Calculate the shock sensor auxiliary values
                const auto aux_1 = r_var_j - r_var_i;
                const auto aux_2 = 0.5 * inner_prod(l_ji, r_grad_var_i + r_grad_var_j);
                const auto num = aux_1 - aux_2;
                const auto den = std::abs(aux_1) + std::abs(aux_2);

                // Check if the solution is not constant (den close to 0.0)
                double beta_ij = 0.0;
                if (std::abs(den) > zero_tol) {
                    // Compute and bound the shock sensor
                    const double aux_beta_ij = std::abs(num / den);
                    beta_ij = aux_beta_ij < 1.0 ? aux_beta_ij : 1.0;
                }

                // Check against the current value of shock sensor and keep the largest one
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
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override;

    ///@}
private:
    ///@name Member Variables
    ///@{

    /// Reference to the model part in where the shock detection is to be performed
    ModelPart& mrModelPart;

    /// Updates the NODAL_AREA at each time step (required in case the mesh deforms)
    const bool mUpdateNodalAreaAtEachStep = false;

    /// Updates the NODAL_NEIGHBOURS at each time step (required in case topology changes)
    const bool mUpdateNodalNeighboursAtEachStep = false;

    /// Flag to indicate if the nodal area has been already computed
    bool mNodalAreaAlreadyComputed = false;

    /// Flag to indicate if the nodal neighbours have been already computed
    bool mNodalNeighboursAlreadyComputed = false;

    /// Flag to indicate if the shock variable type is double or component one
    const bool mShockVariableIsDouble;

    /// Pointer to the shock detection double variable
    const Variable<double>* mpShockDoubleVariable = nullptr;

    /// Pointer to the shock detection component variable
    const VariableComponentType* mpShockComponentVariable = nullptr;

    /// Name of the shock detection gradient variable
    const Variable<array_1d<double,3>>* mpShockGradientVariable = nullptr;

    /// Name of the shock sensor variable
    const Variable<double>* mpShockSensorVariable = nullptr;

    ///@}
    ///@name Serialization
    ///@{


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
