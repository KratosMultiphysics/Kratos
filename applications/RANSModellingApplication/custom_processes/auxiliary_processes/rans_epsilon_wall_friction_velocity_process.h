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

#if !defined(KRATOS_RANS_EPSILON_WALL_FRICTION_VELOCITY_PROCESS_H_INCLUDED)
#define KRATOS_RANS_EPSILON_WALL_FRICTION_VELOCITY_PROCESS_H_INCLUDED

// System includes
#include <cmath>
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

class RansEpsilonWallFrictionVelocityProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    using NodesContainerType = ModelPart::NodesContainerType;

    /// Pointer definition of RansEpsilonWallFrictionVelocityProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansEpsilonWallFrictionVelocityProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansEpsilonWallFrictionVelocityProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "von_karman"              : 0.41,
            "echo_level"              : 0,
            "is_fixed"                : true,
            "min_epsilon_value"       : 1e-18
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mVonKarman = mrParameters["von_karman"].GetDouble();
        mMinEpsilonValue = mrParameters["min_epsilon_value"].GetDouble();
        mIsConstrained = mrParameters["is_fixed"].GetBool();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansEpsilonWallFrictionVelocityProcess() override
    {
        // delete mpDistanceCalculator;
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
        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
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
        return std::string("RansEpsilonWallFrictionVelocityProcess");
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

    double mMinEpsilonValue;
    double mVonKarman;
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
        const double y_plus = rNode.FastGetSolutionStepValue(RANS_Y_PLUS);
        const double nu = rNode.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double wall_velocity_magnitude = norm_2(rNode.FastGetSolutionStepValue(VELOCITY));

        const double u_tau = wall_velocity_magnitude / (std::log(y_plus)/0.41 + 5.2);

        double& epsilon = rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
        epsilon = std::max(std::pow(u_tau, 4) / (mVonKarman * nu * y_plus), mMinEpsilonValue);
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
    RansEpsilonWallFrictionVelocityProcess& operator=(RansEpsilonWallFrictionVelocityProcess const& rOther);

    /// Copy constructor.
    RansEpsilonWallFrictionVelocityProcess(RansEpsilonWallFrictionVelocityProcess const& rOther);

    ///@}

}; // Class RansEpsilonWallFrictionVelocityProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEpsilonWallFrictionVelocityProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EPSILON_WALL_FRICTION_VELOCITY_PROCESS_H_INCLUDED defined
