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

#include "custom_utilities/rans_calculation_utilities.h"
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
 * @brief Calculates y_plus value base on the logarithmic law
 *
 * This process calculates $y^+$ value based on the following formula:
 *
 * \[
 * 	u^+ = \frac{||\vvel||}{\vel_\tau} =
 * \begin{cases}
 *	\frac{1}{\kappa}ln\left(y^+\right) + \beta &\text{ for } y^+ > y^+_{limit} \\
 *  y^+ &\text{ for } y^+ \leq y^+_{limit}
 * \end{cases}
 * \]
 * Where,
 * \[
 *	y^+ = \frac{\vel_\tau y}{\nu}
 * \]
 * \[
 *	y^+_{limit} = \frac{1}{\kappa}ln\left(y^+_{limit}\right) + \beta = \vel^+
 * \]
 *
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

        mMaxIterations = mrParameters["max_iterations"].GetInt();
        mTolerance = mrParameters["tolerance"].GetDouble();

        mVonKarman = mrParameters["constants"]["von_karman"].GetDouble();
        mBeta = mrParameters["constants"]["beta"].GetDouble();

        mLimitYPlus = RansCalculationUtilities().CalculateLogarithmicYPlusLimit(
            mVonKarman, mBeta);

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansLogarithmicYPlusCalculationProcess() override
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

        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);

        RansCheckUtilities rans_check_utilities;

        rans_check_utilities.CheckIfModelPartExists(mrModel, mModelPartName);

        const ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();

        rans_check_utilities.CheckIfVariableExistsInNodesContainer(r_nodes, VELOCITY);
        rans_check_utilities.CheckIfVariableExistsInNodesContainer(r_nodes, DISTANCE);
        rans_check_utilities.CheckIfVariableExistsInNodesContainer(r_nodes, KINEMATIC_VISCOSITY);
        rans_check_utilities.CheckIfVariableExistsInNodesContainer(r_nodes, RANS_Y_PLUS);

        return 0;

        KRATOS_CATCH("");
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

    unsigned int mMaxIterations;
    double mTolerance;

    double mVonKarman;
    double mBeta;
    double mLimitYPlus;

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
            rNode.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double wall_distance = rNode.FastGetSolutionStepValue(DISTANCE);
        const array_1d<double, 3>& velocity =
            rNode.FastGetSolutionStepValue(rVelocityVariable);
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
        const double inv_von_karman = 1.0 / mVonKarman;

        // log region
        if (yplus > mLimitYPlus)
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
