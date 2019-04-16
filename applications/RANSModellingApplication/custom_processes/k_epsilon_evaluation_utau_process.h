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

#if !defined(KRATOS_RANS_K_EPSILON_EVALUATION_UTAU_PROCESS_H_INCLUDED)
#define KRATOS_RANS_K_EPSILON_EVALUATION_UTAU_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
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

/// Auxiliary process to set Boussinesq buoyancy forces in variable temperature flows.
/** This process modifies the BODY_FORCE variable according to the Boussinesq hypothesis
    so that the fluid element can take natural convection into account.

    This process makes use of the following data:
    - TEMPERATURE from the nodal solution step data: current temperature for the node (mandatory).
    - AMBIENT_TEMPERATURE from ProcessInfo: The reference temperature for the simulation (mandatory).
    - gravity from the Parameters passed in the constructor: an array that defines the gravity vector (mandatory).
    - thermal_expansion_coefficient from the Parameters: a double defining the thermal expansion coefficient for the fluid (optional).

    With this, the process calculates the Boussinesq force and assings it to the BODY_FORCE solution step variable of each node.
    The force is set to (1 + thermal_expansion_coefficient*(temperature - ambient_temperature) ) * g

    If the thermal expansion coefficient is not provided, it is assumed to be (1/ambient_temperature).
    This is the usual value for perfect gases (if the temperature is given in Kelvin).
 */

class RansKEpsilonEvaluationUtauProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    /// Pointer definition of RansKEpsilonEvaluationUtauProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansKEpsilonEvaluationUtauProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansKEpsilonEvaluationUtauProcess(Model& rModel, Parameters& rParameters, Process& rRansYPlusModel)
        : mrModel(rModel), mrParameters(rParameters), mrRansYPlusModel(rRansYPlusModel)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "constants": {
                "von_karman"  : 0.41,
                "c_mu"        : 0.09
            },
            "variable_minimums": {
                "turbulent_kinetic_energy"          : 1e-10,
                "turbulent_energy_dissipation_rate" : 1e-10
            },
            "is_initialization_process": false
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mVonKarman = mrParameters["constants"]["von_karman"].GetDouble();
        mCmu = mrParameters["constants"]["c_mu"].GetDouble();
        mTurbulentKineticEnergyMin =
            mrParameters["variable_minimums"]["turbulent_kinetic_energy"].GetDouble();
        mTurbulentEnergyDissipationMin =
            mrParameters["variable_minimums"]
                        ["turbulent_energy_dissipation_rate"]
                            .GetDouble();

        mModelPartName = mrParameters["model_part_name"].GetString();

        mIsInitializationProcess = mrParameters["is_initialization_process"].GetBool();
        mIsValuesAssigned = false;

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansKEpsilonEvaluationUtauProcess() override
    {
        // delete mpDistanceCalculator;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY

        // Calculates y plus values;
        mrRansYPlusModel.Execute();

        ModelPart& r_modelpart = mrModel.GetModelPart(mModelPartName);
        NodesContainerType& r_nodes = r_modelpart.Nodes();

        if (mIsInitializationProcess && !mIsValuesAssigned)
        {
            ApplyValuesForFreeNodes(r_nodes);
            mIsValuesAssigned = true;
        }
        else if (!mIsInitializationProcess)
            ApplyValuesForAllNodes(r_nodes);

        KRATOS_CATCH("");
    }

    int Check() override
    {
        // Checking variable definitions
        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);

        if (!mrModel.HasModelPart(mModelPartName))
        {
            const std::vector<std::string>& r_model_part_names =
                mrModel.GetModelPartNames();

            std::string msg;
            msg = mrModel.Info() + " doesn't have " + mModelPartName +
                  ". Available model parts are: \n";
            for (std::string model_part_name : r_model_part_names)
                msg += "     " + model_part_name + "\n";

            KRATOS_ERROR << msg;
        }

        ModelPart& r_modelpart = mrModel.GetModelPart(mModelPartName);
        ModelPart::NodesContainerType& r_nodes = r_modelpart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        }

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
        return std::string("RansKEpsilonEvaluationUtauProcess");
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
    Process& mrRansYPlusModel;

    double mVonKarman;
    double mCmu;
    double mTurbulentKineticEnergyMin;
    double mTurbulentEnergyDissipationMin;
    std::string mModelPartName;

    bool mIsInitializationProcess;
    bool mIsValuesAssigned;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ApplyValuesForAllNodes(NodesContainerType& rNodes)
    {
        int number_of_nodes = rNodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(rNodes.begin() + i_node);
            CalculateTurbulentValues(r_node);
            r_node.Fix(TURBULENT_KINETIC_ENERGY);
            r_node.Fix(TURBULENT_ENERGY_DISSIPATION_RATE);
        }
    }

    void ApplyValuesForFreeNodes(NodesContainerType& rNodes)
    {
        int number_of_nodes = rNodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(rNodes.begin() + i_node);

            if (!r_node.IsFixed(TURBULENT_KINETIC_ENERGY) && !r_node.IsFixed(TURBULENT_ENERGY_DISSIPATION_RATE))
                CalculateTurbulentValues(r_node);
        }
    }

    void CalculateTurbulentValues(NodeType& rNode)
    {
        const double y_plus = rNode.FastGetSolutionStepValue(RANS_Y_PLUS);
        const double nu = rNode.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double wall_distance = rNode.FastGetSolutionStepValue(DISTANCE);

        double& tke = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        double& epsilon = rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);

        const double u_tau = y_plus * nu / wall_distance;
        tke = std::max(mTurbulentKineticEnergyMin, std::pow(u_tau, 2) / std::sqrt(mCmu));
        epsilon = std::max(mTurbulentEnergyDissipationMin,
                           std::pow(u_tau, 3) / (mVonKarman * wall_distance));
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
    RansKEpsilonEvaluationUtauProcess& operator=(RansKEpsilonEvaluationUtauProcess const& rOther);

    /// Copy constructor.
    RansKEpsilonEvaluationUtauProcess(RansKEpsilonEvaluationUtauProcess const& rOther);

    ///@}

}; // Class RansKEpsilonEvaluationUtauProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansKEpsilonEvaluationUtauProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_K_EPSILON_EVALUATION_UTAU_PROCESS_H_INCLUDED  defined
