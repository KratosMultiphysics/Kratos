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

#if !defined(KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_SENSITIVITIES_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_SENSITIVITIES_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "rans_modelling_application_variables.h"

#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_check_utilities.h"

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

/**
 * @brief Calculates partial derivatives turbulent kinematic viscosity
 *
 * This process uses following formula to calculate turbulent kinematic viscosity (i.e. $\nu_t$) partial derivatives.
 * $\nu_t$ calculation is according to RansNutKEpsilonHighReSensitivitiesProcess class
 *
 * \[
 *      \frac{\partial\nu_t}{\partial u_x} = 0 \\
 *      \frac{\partial\nu_t}{\partial u_y} = 0 \\
 *      \frac{\partial\nu_t}{\partial u_z} = 0 \\
 *      \frac{\partial\nu_t}{\partial P} = 0 \\
 *      \frac{\partial\nu_t}{\partial k} = 2C_\mu\frac{k}{\epsilon} \\
 *      \frac{\partial\nu_t}{\partial \epsilon} = -C_\mu \frac{k^2}{\epsilon^2} \\
 * \]
 *
 * $k$ is the turbulent kinetic energy, $\epsilon$ is the turbulent energy dissipation rate
 * $u_x$, $u_y$, $u_z$, $P$ are velocity componenets and pressure.
 *
 * RANS_NUT_SCALAR_PARTIAL_DERIVATIVES vector variable is filled with partial derivatives of each node in following order
 *      Index 0: partial derivative of \nu_t w.r.t. $k$
 *      Index 1: partial derivative of \nu_t w.r.t. $\epsilon$
 *
 * This variable is only stored in nodal data value container (not in the historical data value container). The velocity
 * and pressure derivatives are not stored since n the $k-\epsilon$ high $Re$ implementation, they are always zero.
 *
 * @see RansNutKEpsilonHighReCalculationProcess
 */

class RansNutKEpsilonHighReSensitivitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansNutKEpsilonHighReSensitivitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutKEpsilonHighReSensitivitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutKEpsilonHighReSensitivitiesProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mCmu = mrParameters["c_mu"].GetDouble();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansNutKEpsilonHighReSensitivitiesProcess() override
    {
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

        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);

        RansCheckUtilities rans_check_utilities;

        rans_check_utilities.CheckIfModelPartExists(mrModel, mModelPartName);

        const ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();

        rans_check_utilities.CheckIfVariableExistsInNodesContainer(
            r_nodes, TURBULENT_KINETIC_ENERGY);
        rans_check_utilities.CheckIfVariableExistsInNodesContainer(
            r_nodes, TURBULENT_ENERGY_DISSIPATION_RATE);

        return 0;

        KRATOS_CATCH("");
    }

    void ExecuteInitializeSolutionStep() override
    {
        Execute();
    }

    void Execute() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        ModelPart::NodesContainerType& r_nodes = r_model_part.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            Vector nut_partial_derivatives(2);

            const double& tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double& epsilon =
                r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

            nut_partial_derivatives[0] = 2.0 * mCmu * tke / epsilon;
            nut_partial_derivatives[1] = -1.0 * mCmu * std::pow(tke / epsilon, 2);

            r_node.SetValue(RANS_NUT_SCALAR_PARTIAL_DERIVATIVES, nut_partial_derivatives);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << "Calculated k-epsilon high Re nu_t sensitivities for nodes in"
            << mModelPartName << "\n";

        KRATOS_CATCH("");
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
        return std::string("RansNutKEpsilonHighReSensitivitiesProcess");
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

    int mEchoLevel;

    double mCmu;

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
    RansNutKEpsilonHighReSensitivitiesProcess& operator=(
        RansNutKEpsilonHighReSensitivitiesProcess const& rOther);

    /// Copy constructor.
    RansNutKEpsilonHighReSensitivitiesProcess(RansNutKEpsilonHighReSensitivitiesProcess const& rOther);

    ///@}

}; // Class RansNutKEpsilonHighReSensitivitiesProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansNutKEpsilonHighReSensitivitiesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_SENSITIVITIES_PROCESS_H_INCLUDED defined
