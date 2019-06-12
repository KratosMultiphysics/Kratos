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

#if !defined(KRATOS_RANS_LOGARITHMIC_Y_PLUS_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_LOGARITHMIC_Y_PLUS_CALCULATION_PROCESS_H_INCLUDED

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

class RansLogarithmicYPlusCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansLogarithmicYPlusCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansLogarithmicYPlusCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansLogarithmicYPlusCalculationProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "step"            : 0,
            "max_iterations"  : 100,
            "tolerance"       : 1e-6,
            "constants": {
                "von_karman"  : 0.41,
                "beta"        : 5.2
            }
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mModelPartName = mrParameters["model_part_name"].GetString();

        mEchoLevel = mrParameters["echo_level"].GetInt();
        mStep = mrParameters["step"].GetInt();

        mMaxIterations = mrParameters["max_iterations"].GetInt();
        mTolerance = mrParameters["tolerance"].GetDouble();

        mVonKarman = mrParameters["constants"]["von_karman"].GetDouble();
        mBeta = mrParameters["constants"]["beta"].GetDouble();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansLogarithmicYPlusCalculationProcess() override
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
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);

        const ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
        }

        return 0;
    }

    void Execute() override
    {
        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();

        const int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            CalculateLogarithmicWallLawYplus(r_node, VELOCITY);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "RANS_Y_PLUS calculated for nodes in " << mModelPartName << ".\n";
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
        return std::string("RansLogarithmicYPlusCalculationProcess");
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

    unsigned int mEchoLevel;
    unsigned int mStep;

    unsigned int mMaxIterations;
    double mTolerance;

    double mVonKarman;
    double mBeta;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateLogarithmicWallLawYplus(NodeType& rNode,
                                          const Variable<array_1d<double, 3>>& rVelocityVariable)
    {
        KRATOS_TRY

        const double kinematic_viscosity =
            rNode.FastGetSolutionStepValue(KINEMATIC_VISCOSITY, mStep);
        const double wall_distance = rNode.FastGetSolutionStepValue(DISTANCE, mStep);
        const array_1d<double, 3>& velocity =
            rNode.FastGetSolutionStepValue(rVelocityVariable, mStep);
        const double velocity_magnitude = norm_2(velocity);

        if (velocity_magnitude < std::numeric_limits<double>::epsilon())
        {
            rNode.FastGetSolutionStepValue(RANS_Y_PLUS) = 0.0;
            return;
        }

        KRATOS_ERROR_IF(wall_distance < std::numeric_limits<double>::epsilon())
            << "DISTANCE at node " << rNode.Coordinates()
            << " with id=" << rNode.Id() << " has zero value. Please specify DISTANCE value > 0.0 for all the nodes in "
            << mModelPartName << " to calculate RANS_Y_PLUS.\n";

        // linear region
        double utau = sqrt(velocity_magnitude * kinematic_viscosity / wall_distance);
        double yplus = wall_distance * utau / kinematic_viscosity;

        const double limit_yplus = 11.06;
        const double inv_von_karman = 1.0 / mVonKarman;

        // log region
        if (yplus > limit_yplus)
        {
            unsigned int iter = 0;
            double dx = 1e10;
            const double tol = mTolerance;
            double uplus = inv_von_karman * log(yplus) + mBeta;

            while (iter < mMaxIterations && fabs(dx) > tol * utau)
            {
                // Newton-Raphson iteration
                double f = utau * uplus - velocity_magnitude;
                double df = uplus + inv_von_karman;
                dx = f / df;

                // Update variables
                utau -= dx;
                yplus = wall_distance * utau / kinematic_viscosity;
                uplus = inv_von_karman * log(yplus) + mBeta;
                ++iter;
            }

            KRATOS_WARNING_IF("RansLogarithmicYPlusCalculationProcess", iter == mMaxIterations && mEchoLevel > 0)
                << "Y plus calculation in Wall (logarithmic region) "
                   "Newton-Raphson did not converge. "
                   "residual > tolerance [ "
                << std::scientific << dx << " > " << std::scientific << tol << " ]\n";
        }
        rNode.FastGetSolutionStepValue(RANS_Y_PLUS) = yplus;

        KRATOS_CATCH("");
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
    RansLogarithmicYPlusCalculationProcess& operator=(RansLogarithmicYPlusCalculationProcess const& rOther);

    /// Copy constructor.
    RansLogarithmicYPlusCalculationProcess(RansLogarithmicYPlusCalculationProcess const& rOther);

    ///@}

}; // namespace Kratos

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansLogarithmicYPlusCalculationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_LOGARITHMIC_Y_PLUS_CALCULATION_PROCESS_H_INCLUDED  defined
