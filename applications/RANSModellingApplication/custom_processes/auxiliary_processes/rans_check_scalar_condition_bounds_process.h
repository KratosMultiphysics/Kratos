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

#if !defined(KRATOS_RANS_CHECK_SCALAR_CONDITION_BOUNDS_PROCESS_H_INCLUDED)
#define KRATOS_RANS_CHECK_SCALAR_CONDITION_BOUNDS_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "containers/model.h"
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

class RansCheckScalarConditionBoundsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansCheckScalarConditionBoundsProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansCheckScalarConditionBoundsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansCheckScalarConditionBoundsProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"            : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "gauss_integration_method" : "GAUSS_1",
            "echo_level"      : 0
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = mrParameters["variable_name"].GetString();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansCheckScalarConditionBoundsProcess() override
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

        const Variable<double> scalar_variable =
            KratosComponents<Variable<double>>::Get(mVariableName);

        KRATOS_CHECK_VARIABLE_KEY(scalar_variable);

        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(scalar_variable, r_node);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Check passed for " << mModelPartName << " with variable "
            << scalar_variable.Name() << ".\n";

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        KRATOS_TRY

        const Variable<double> scalar_variable =
            KratosComponents<Variable<double>>::Get(mVariableName);

        const ModelPart::ConditionsContainerType& r_conditions =
            mrModel.GetModelPart(mModelPartName).Conditions();

        RansCalculationUtilities rans_calculation_utilities;

        double scalar_min = 0.0;
        double scalar_max = 0.0;
        bool initialized = false;

        // #pragma omp parallel for
        const int number_of_conditions = r_conditions.size();
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
        {
            const ModelPart::ConditionType& r_condition = *(r_conditions.begin() + i_cond);

            KRATOS_ERROR_IF(!r_condition.Has(PARENT_ELEMENT))
                << "Parent element not found for condition id=" << r_condition.Id()
                << ". Please run \"FindConditionParentProcess\" for "
                << mModelPartName << ".\n";

            const ModelPart::ElementType::WeakPointer& p_parent_element =
                r_condition.GetValue(PARENT_ELEMENT);
            const Element::GeometryType& r_parent_element_geometry =
                p_parent_element->GetGeometry();

            Vector gauss_weights;
            Matrix shape_functions;
            ModelPart::ElementType::GeometryType::ShapeFunctionsGradientsType shape_function_derivatives;

            rans_calculation_utilities.CalculateGeometryData(
                r_parent_element_geometry, GeometryData::GI_GAUSS_1,
                gauss_weights, shape_functions, shape_function_derivatives);

            const Vector& gauss_shape_functions = row(shape_functions, 0);
            const double current_value = rans_calculation_utilities.EvaluateInPoint(
                r_parent_element_geometry, scalar_variable, gauss_shape_functions);

            if (!initialized)
            {
                scalar_max = current_value;
                scalar_min = current_value;
                initialized = true;
            }

            if (scalar_max < current_value)
                scalar_max = current_value;
            if (scalar_min > current_value)
                scalar_min = current_value;

        }

        KRATOS_INFO(this->Info())
            << scalar_variable.Name() << " is bounded between [ " << scalar_min
            << ", " << scalar_max << " ] in " << mModelPartName << ".\n";

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
        return std::string("RansCheckScalarConditionBoundsProcess");
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
    std::string mVariableName;
    int mEchoLevel;

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
    RansCheckScalarConditionBoundsProcess& operator=(RansCheckScalarConditionBoundsProcess const& rOther);

    /// Copy constructor.
    RansCheckScalarConditionBoundsProcess(RansCheckScalarConditionBoundsProcess const& rOther);

    ///@}

}; // Class RansCheckScalarConditionBoundsProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansCheckScalarConditionBoundsProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_CHECK_SCALAR_CONDITION_BOUNDS_PROCESS_H_INCLUDED defined
