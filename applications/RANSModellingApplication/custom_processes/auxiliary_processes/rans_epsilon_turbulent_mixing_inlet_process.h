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

#if !defined(KRATOS_RANS_EPSILON_TURBULENT_MIXING_LENGTH_INLET_PROCESS_H_INCLUDED)
#define KRATOS_RANS_EPSILON_TURBULENT_MIXING_LENGTH_INLET_PROCESS_H_INCLUDED

// System includes
#include <cmath>
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "rans_modelling_application_variables.h"

#include "custom_utilities/rans_check_utilities.h"
#include "custom_utilities/rans_variable_utils.h"

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
 * @brief Sets epsilon value best on turbulent mixing length
 *
 * This process sets epsilon values based on the following formula
 *
 * \[
 *
 *  \epsilon = \frac{C_\mu^{0.75} max\left\lbrace k, 0.0\right\rbrace^{1.5}}{L}
 *
 * \]
 *
 * In here $k$ is turbulent kinetic energy, $\epsilon$ is turbulent energy
 * dissipation rate, and $L$ is turbulent mixing length.
 *
 */

class RansEpsilonTurbulentMixingLengthInletProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    /// Pointer definition of RansEpsilonTurbulentMixingLengthInletProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansEpsilonTurbulentMixingLengthInletProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansEpsilonTurbulentMixingLengthInletProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "turbulent_mixing_length" : 0.005,
            "c_mu"                    : 0.09,
            "echo_level"              : 0,
            "is_fixed"                : true
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mTurbulentMixingLength = mrParameters["turbulent_mixing_length"].GetDouble();
        mIsConstrained = mrParameters["is_fixed"].GetBool();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mCmu_75 = std::pow(mrParameters["c_mu"].GetDouble(), 0.75);

        KRATOS_ERROR_IF(mTurbulentMixingLength < std::numeric_limits<double>::epsilon())
            << "turbulent_mixing_length should be greater than zero.\n";

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansEpsilonTurbulentMixingLengthInletProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override
    {
        if (mIsConstrained)
        {
            ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

            const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
                r_node.Fix(TURBULENT_ENERGY_DISSIPATION_RATE);
            }

            KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
                << "Fixed TURBULENT_ENERGY_DISSIPATION_RATE dofs in "
                << mModelPartName << ".\n";
        }
    }

    void Execute() override
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            CalculateTurbulentValues(r_node);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Applied epsilon values to " << mModelPartName << ".\n";

        KRATOS_CATCH("");
    }

    int Check() override
    {
        // Checking variable definitions
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);

        RansCheckUtilities rans_check_utilities;

        rans_check_utilities.CheckIfModelPartExists(mrModel, mModelPartName);

        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();

        rans_check_utilities.CheckIfVariableExistsInNodesContainer(
            r_nodes, TURBULENT_KINETIC_ENERGY);
        rans_check_utilities.CheckIfVariableExistsInNodesContainer(
            r_nodes, TURBULENT_ENERGY_DISSIPATION_RATE);

        return 0;
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
        return std::string("RansEpsilonTurbulentMixingLengthInletProcess");
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

    double mTurbulentMixingLength;
    double mCmu_75;
    bool mIsConstrained;
    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateTurbulentValues(NodeType& rNode)
    {
        const double tke = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) =
            mCmu_75 * std::pow(std::max(tke, 0.0), 1.5) / mTurbulentMixingLength;
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
    RansEpsilonTurbulentMixingLengthInletProcess& operator=(
        RansEpsilonTurbulentMixingLengthInletProcess const& rOther);

    /// Copy constructor.
    RansEpsilonTurbulentMixingLengthInletProcess(RansEpsilonTurbulentMixingLengthInletProcess const& rOther);

    ///@}

}; // Class RansEpsilonTurbulentMixingLengthInletProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEpsilonTurbulentMixingLengthInletProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EPSILON_TURBULENT_MIXING_LENGTH_INLET_PROCESS_H_INCLUDED defined
