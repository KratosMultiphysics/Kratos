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

#if !defined(KRATOS_RANS_Y_PLUS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_Y_PLUS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/linear_solver_factory.h"
#include "includes/model_part.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/process.h"
#include "rans_modelling_application_variables.h"
#include "utilities/normal_calculation_utils.h"

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

class RansYPlusWallDistanceCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansYPlusWallDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansYPlusWallDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansYPlusWallDistanceCalculationProcess(ModelPart& rModelPart, Parameters& rParameters)
        : mrModelPart(rModelPart), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "echo_level"            : 0,
            "velocity_variable"     : "WALL_VELOCITY",
            "limit_y_plus"          : 11.06,
            "von_karman"            : 0.41,
            "beta"                  : 5.2,
            "minimum_wall_distance" : 1e-5
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = mrParameters["echo_level"].GetInt();
        mVonKarman = mrParameters["von_karman"].GetDouble();
        mBeta = mrParameters["beta"].GetDouble();
        mMinWallDistance = mrParameters["minimum_wall_distance"].GetDouble();
        mWallVelocityVariableName = mrParameters["velocity_variable"].GetString();
        mLimitYPlus = RansCalculationUtilities().CalculateLogarithmicYPlusLimit(
            mVonKarman, mBeta);

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansYPlusWallDistanceCalculationProcess() override
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

        const Variable<array_1d<double, 3>> velocity_variable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mWallVelocityVariableName);

        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(velocity_variable);
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(velocity_variable, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        }

        return 0;

        KRATOS_CATCH("");
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        const Variable<array_1d<double, 3>> velocity_variable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mWallVelocityVariableName);

        const int number_of_nodes = mrModelPart.NumberOfNodes();
        const double inv_von_karman = 1.0 / mVonKarman;

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);

            double& r_distance = r_node.FastGetSolutionStepValue(DISTANCE);

            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            const array_1d<double, 3>& r_velocity =
                r_node.FastGetSolutionStepValue(velocity_variable);
            const double velocity_magnitude = norm_2(r_velocity);

            if (velocity_magnitude > std::numeric_limits<double>::epsilon())
            {
                if (y_plus < mLimitYPlus)
                {
                    r_distance = std::pow(y_plus, 2) * nu / velocity_magnitude;
                }
                else
                {
                    const double u_tau = velocity_magnitude /
                                         (inv_von_karman * std::log(y_plus) + mBeta);
                    r_distance = y_plus * nu / u_tau;
                }
            }
        }

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
        return std::string("RansYPlusWallDistanceCalculationProcess");
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

    ModelPart& mrModelPart;
    Parameters& mrParameters;

    int mEchoLevel;
    double mLimitYPlus;
    double mVonKarman;
    double mBeta;
    double mMinWallDistance;

    std::string mWallVelocityVariableName;

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
    RansYPlusWallDistanceCalculationProcess& operator=(RansYPlusWallDistanceCalculationProcess const& rOther);

    /// Copy constructor.
    RansYPlusWallDistanceCalculationProcess(RansYPlusWallDistanceCalculationProcess const& rOther);

    ///@}

}; // Class RansYPlusWallDistanceCalculationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansYPlusWallDistanceCalculationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_Y_PLUS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED defined
