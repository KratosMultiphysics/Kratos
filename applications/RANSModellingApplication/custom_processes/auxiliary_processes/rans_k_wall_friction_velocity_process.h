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

#if !defined(KRATOS_RANS_K_WALL_FRICTION_VELOCITY_PROCESS_H_INCLUDED)
#define KRATOS_RANS_K_WALL_FRICTION_VELOCITY_PROCESS_H_INCLUDED

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

class RansKWallFrictionVelocityProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    /// Pointer definition of RansKWallFrictionVelocityProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansKWallFrictionVelocityProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansKWallFrictionVelocityProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"     : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "c_mu"                : 0.09,
            "echo_level"          : 0,
            "is_fixed"            : true,
            "min_k_value"         : 1e-18
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mCmu = mrParameters["c_mu"].GetDouble();
        mIsConstrained = mrParameters["is_fixed"].GetBool();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mMinValue = mrParameters["min_k_value"].GetDouble();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansKWallFrictionVelocityProcess() override
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
            ModelPart::NodesContainerType& r_nodes =
                mrModel.GetModelPart(mModelPartName).Nodes();
            int number_of_nodes = r_nodes.size();
#pragma omp parallel for
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_node = *(r_nodes.begin() + i_node);
                r_node.Fix(TURBULENT_KINETIC_ENERGY);
            }

            KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
                << "Fixed TURBULENT_KINETIC_ENERGY dofs in " << mModelPartName << ".\n";
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
            << "Applied k values to " << mModelPartName << ".\n";

        KRATOS_CATCH("");
    }

    int Check() override
    {
        // Checking variable definitions
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);

        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
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
        return std::string("RansKWallFrictionVelocityProcess");
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

    double mCmu;
    double mMinValue;
    int mEchoLevel;

    bool mIsConstrained;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateTurbulentValues(NodeType& rNode)
    {
        const double y_plus = rNode.FastGetSolutionStepValue(RANS_Y_PLUS);
        const double wall_velocity_magnitude = norm_2(rNode.FastGetSolutionStepValue(VELOCITY));

        const double u_tau = wall_velocity_magnitude / (std::log(y_plus)/0.41 + 5.2);

        rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) =
            std::max(mMinValue, std::pow(u_tau, 2) / std::sqrt(mCmu));
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
    RansKWallFrictionVelocityProcess& operator=(RansKWallFrictionVelocityProcess const& rOther);

    /// Copy constructor.
    RansKWallFrictionVelocityProcess(RansKWallFrictionVelocityProcess const& rOther);

    ///@}

}; // Class RansKWallFrictionVelocityProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansKWallFrictionVelocityProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_K_WALL_FRICTION_VELOCITY_PROCESS_H_INCLUDED defined
