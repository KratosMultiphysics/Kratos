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

#if !defined(KRATOS_RANS_LOGARITHMIC_Y_PLUS_VELOCITY_SENSITIVITIES_PROCESS_H_INCLUDED)
#define KRATOS_RANS_LOGARITHMIC_Y_PLUS_VELOCITY_SENSITIVITIES_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "processes/process.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"
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

/**
 * @brief Calculates $y^+$ sensitivities w.r.t. velocity
 *
 * This process calculates sensitivities of $y^+$ w.r.t. velocity, In here
 * The $y^+$ should be calculated based on logarithmic law.
 *
 * @see RansLogarithmicYPlusCalculationProcess
 */

class RansLogarithmicYPlusVelocitySensitivitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::ElementType ElementType;

    typedef Geometry<NodeType> GeometryType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Pointer definition of RansLogarithmicYPlusVelocitySensitivitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansLogarithmicYPlusVelocitySensitivitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansLogarithmicYPlusVelocitySensitivitiesProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "constants": {
                "von_karman"  : 0.41,
                "beta"        : 5.2
            }
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mModelPartName = mrParameters["model_part_name"].GetString();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        mVonKarman = mrParameters["constants"]["von_karman"].GetDouble();
        mBeta = mrParameters["constants"]["beta"].GetDouble();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansLogarithmicYPlusVelocitySensitivitiesProcess() override
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

        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);
        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

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

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        this->Execute();
    }

    void Execute() override
    {
        KRATOS_TRY

        ModelPart& model_part = mrModel.GetModelPart(mModelPartName);

        ModelPart::ElementsContainerType& r_elements = model_part.Elements();

        const int number_of_elements = r_elements.size();

        const int domain_size = model_part.GetProcessInfo()[DOMAIN_SIZE];

        const double inv_kappa = 1.0 / mVonKarman;

#pragma omp parallel for
        for (int i_element = 0; i_element < number_of_elements; ++i_element)
        {
            ElementType& r_element = *(r_elements.begin() + i_element);
            GeometryType& r_geometry = r_element.GetGeometry();
            const int number_of_nodes = r_geometry.PointsNumber();

            Matrix r_adjoint_y_plus_matrix(number_of_nodes, domain_size);
            r_adjoint_y_plus_matrix.clear();

            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                const NodeType& r_node = r_geometry[i_node];
                const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
                const double wall_distance = r_node.FastGetSolutionStepValue(DISTANCE);
                const array_1d<double, 3> velocity =
                    r_node.FastGetSolutionStepValue(VELOCITY);
                const double velocity_magnitude = norm_2(velocity);
                const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

                double value = 0.0;

                if (y_plus > 11.06)
                {
                    value = (inv_kappa * std::log(y_plus) + mBeta);
                    value = value / (std::pow(value, 2) + velocity_magnitude * wall_distance *
                                                              inv_kappa / (nu * y_plus));
                }
                else
                {
                    value = 1.0 / (2.0 * y_plus);
                }
                if (velocity_magnitude > std::numeric_limits<double>::epsilon() &&
                    y_plus > std::numeric_limits<double>::epsilon())
                {
                    for (int i_dim = 0; i_dim < domain_size; ++i_dim)
                    {
                        r_adjoint_y_plus_matrix(i_node, i_dim) =
                            (wall_distance / nu) * value * velocity[i_dim] / velocity_magnitude;
                    }
                }
            }
            r_element.SetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES, r_adjoint_y_plus_matrix);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "RANS_Y_PLUS_VELOCITY_DERIVATIVES calculated for "
            << r_elements.size() << " elements in " << mModelPartName << ".\n";

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
        return std::string("RansLogarithmicYPlusVelocitySensitivitiesProcess");
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

    double mVonKarman;
    double mBeta;

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
    RansLogarithmicYPlusVelocitySensitivitiesProcess& operator=(
        RansLogarithmicYPlusVelocitySensitivitiesProcess const& rOther);

    /// Copy constructor.
    RansLogarithmicYPlusVelocitySensitivitiesProcess(
        RansLogarithmicYPlusVelocitySensitivitiesProcess const& rOther);

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
                                const RansLogarithmicYPlusVelocitySensitivitiesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_LOGARITHMIC_Y_PLUS_VELOCITY_SENSITIVITIES_PROCESS_H_INCLUDED defined
