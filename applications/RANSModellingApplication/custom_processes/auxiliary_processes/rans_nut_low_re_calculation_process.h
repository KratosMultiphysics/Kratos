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

#if !defined(KRATOS_RANS_NUT_LOW_RE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_LOW_RE_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
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

class RansNutLowReCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    using NodesContainerType = ModelPart::NodesContainerType;

    /// Pointer definition of RansNutLowReCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutLowReCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutLowReCalculationProcess(Model& rModel, Parameters& rParameters)
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
    ~RansNutLowReCalculationProcess() override
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
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        }

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const double c_mu = r_process_info[TURBULENCE_RANS_C_MU];

        NodesContainerType& r_nodes = r_model_part.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double epsilon =
                r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
            const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
                c_mu, tke, epsilon, f_mu);

            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t;
        }

        const double nu_t_min = r_process_info[TURBULENT_VISCOSITY_MIN];
        const double nu_t_max = r_process_info[TURBULENT_VISCOSITY_MAX];

        unsigned int lower_number_of_nodes, higher_number_of_nodes, total_selected_nodes;
        RansVariableUtils().ClipScalarVariable(
            lower_number_of_nodes, higher_number_of_nodes, total_selected_nodes,
            nu_t_min, nu_t_max, TURBULENT_VISCOSITY, r_model_part);

        KRATOS_WARNING_IF(this->Info(),
                          (lower_number_of_nodes > 0 || higher_number_of_nodes > 0) &&
                              this->mEchoLevel > 0)
            << "TURBULENT_VISCOSITY out of bounds. [ " << lower_number_of_nodes
            << " nodes < " << std::scientific << nu_t_min << ", "
            << higher_number_of_nodes << " nodes > " << std::scientific << nu_t_max
            << "  out of total number of " << total_selected_nodes << " nodes. ]\n";

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << "Calculated nu_t for nodes in" << mModelPartName << "\n";

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
        return std::string("RansNutLowReCalculationProcess");
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
    RansNutLowReCalculationProcess& operator=(RansNutLowReCalculationProcess const& rOther);

    /// Copy constructor.
    RansNutLowReCalculationProcess(RansNutLowReCalculationProcess const& rOther);

    ///@}

}; // Class RansNutLowReCalculationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansNutLowReCalculationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_LOW_RE_CALCULATION_PROCESS_H_INCLUDED defined
