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

#if !defined(KRATOS_RANS_K_TURBULENT_INTENSITY_INLET_PROCESS_H_INCLUDED)
#define KRATOS_RANS_K_TURBULENT_INTENSITY_INLET_PROCESS_H_INCLUDED

// System includes
#include <cmath>
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
 * @brief Set turbulent kinetic energy value based on the given turbulent intensity
 *
 * This process sets turbulent kinetic energy of a given model part based on the
 * following equation
 *
 * \[
 *
 *     k = \frac{3}{2}\left(I||\underline{u}||\right)^2
 *
 * \]
 *
 * $k$ is the turbulent kinetic energy, $||\underline{u}||$ is the velocity magnitude,
 * $I$ is the turbulent intensity. If the velocity magnitude is zero, then $k_{min}$ is
 * assigned as the turbulent kinetic energy.
 *
 */

class RansKTurbulentIntensityInletProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    using NodesContainerType = ModelPart::NodesContainerType;

    /// Pointer definition of RansKTurbulentIntensityInletProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansKTurbulentIntensityInletProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansKTurbulentIntensityInletProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"     : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "turbulent_intensity" : 0.05,
            "echo_level"          : 0,
            "is_fixed"            : true,
            "min_k_value"         : 1e-18
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mTurbulentIntensity = mrParameters["turbulent_intensity"].GetDouble();
        mIsConstrained = mrParameters["is_fixed"].GetBool();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mMinValue = mrParameters["min_k_value"].GetDouble();

        KRATOS_ERROR_IF(mTurbulentIntensity < 0.0)
            << "Turbulent intensity needs to be positive in the modelpart "
            << mModelPartName << "\n.";
        KRATOS_ERROR_IF(mMinValue < 0.0)
            << "Minimum turbulent kinetic energy needs to be positive in the "
               "modelpart "
            << mModelPartName << "\n.";

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansKTurbulentIntensityInletProcess() override
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

    void ExecuteInitializeSolutionStep() override
    {
        Execute();
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
        KRATOS_TRY

        // Checking variable definitions
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);

        RansCheckUtilities rans_check_utilities;

        rans_check_utilities.CheckIfModelPartExists(mrModel, mModelPartName);

        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();

        rans_check_utilities.CheckIfVariableExistsInNodesContainer(
            r_nodes, TURBULENT_KINETIC_ENERGY);
        rans_check_utilities.CheckIfVariableExistsInNodesContainer(r_nodes, VELOCITY);

        return 0;

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
        return std::string("RansKTurbulentIntensityInletProcess");
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

    double mTurbulentIntensity;
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
        const array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        double velocity_magnitude = norm_2(r_velocity);

        // if (velocity_magnitude > std::numeric_limits<double>::epsilon())
        // {
            rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) =
                1.5 * std::pow(mTurbulentIntensity * velocity_magnitude, 2);
        // }
        // else
        // {
        //     rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = mMinValue;
        // }
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
    RansKTurbulentIntensityInletProcess& operator=(RansKTurbulentIntensityInletProcess const& rOther);

    /// Copy constructor.
    RansKTurbulentIntensityInletProcess(RansKTurbulentIntensityInletProcess const& rOther);

    ///@}

}; // Class RansKTurbulentIntensityInletProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansKTurbulentIntensityInletProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_K_TURBULENT_INTENSITY_INLET_PROCESS_H_INCLUDED defined
