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

#if !defined(KRATOS_RANS_SCALAR_CELL_CENTER_AVERAGING_PROCESS_H_INCLUDED)
#define KRATOS_RANS_SCALAR_CELL_CENTER_AVERAGING_PROCESS_H_INCLUDED

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
#include "includes/model_part.h"
#include "processes/find_nodal_neighbours_process.h"
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

class RansScalarCellCenterAveragingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansScalarCellCenterAveragingProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansScalarCellCenterAveragingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansScalarCellCenterAveragingProcess(ModelPart& rModelPart, Parameters& rParameters)
        : mrModelPart(rModelPart), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "average_neighbour_elements"      : 10,
            "average_neighbour_nodes"         : 10,
            "echo_level"                 : 0,
            "flag_variable_name"         : "PLEASE_SPECIFY_FLAG_VARIABLE_NAME",
            "flag_variable_value"        : true,
            "averaging_variables" : {
                "input_variable_name"  : "PLEASE_SPECIFY_AVERAGING_VARIABLE",
                "output_variable_name" : "PLEASE_SPECIFY_OUTPUT_VARIABLE"
            }
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mFlagVariableName = mrParameters["flag_variable_name"].GetString();
        mAveragingInputVariableName =
            mrParameters["averaging_variables"]["input_variable_name"].GetString();
        mAveragingOutputVariableName =
            mrParameters["averaging_variables"]["output_variable_name"].GetString();
        mFlagVariableValue = mrParameters["flag_variable_value"].GetBool();
        mAverageNeighbourElements = mrParameters["average_neighbour_elements"].GetInt();
        mAverageNeighbourNodes = mrParameters["average_neighbour_nodes"].GetInt();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansScalarCellCenterAveragingProcess() override
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

        const Variable<double>& r_input_variable =
            KratosComponents<Variable<double>>::Get(mAveragingInputVariableName);
        const Variable<double>& r_output_variable =
            KratosComponents<Variable<double>>::Get(mAveragingOutputVariableName);

        KRATOS_CHECK_VARIABLE_KEY(NEIGHBOUR_ELEMENTS);
        KRATOS_CHECK_VARIABLE_KEY(PARENT_ELEMENT);
        KRATOS_CHECK_VARIABLE_KEY(NUMBER_OF_NEIGHBOUR_CONDITIONS);
        KRATOS_CHECK_VARIABLE_KEY(r_input_variable);
        KRATOS_CHECK_VARIABLE_KEY(r_output_variable);

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_input_variable, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_output_variable, r_node);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << "Check passed for " << mrModelPart.Name() << ".\n";

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        SetConditionsFlag();
        FindNodalNeighbours();
        ClearExistingData();
        CalculateCellCenterAverage();
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
        return std::string("RansScalarCellCenterAveragingProcess");
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

    int mAverageNeighbourElements;
    int mAverageNeighbourNodes;
    int mEchoLevel;

    std::string mFlagVariableName;
    bool mFlagVariableValue;
    std::string mAveragingInputVariableName;
    std::string mAveragingOutputVariableName;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ClearExistingData()
    {
        KRATOS_TRY

        RansVariableUtils rans_variable_utilities;
        const Variable<double>& r_variable =
            KratosComponents<Variable<double>>::Get(mAveragingOutputVariableName);
        const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);
        rans_variable_utilities.SetScalarVarForFlag(
            r_variable, 0.0, mrModelPart.Nodes(), r_flag, mFlagVariableValue);
        rans_variable_utilities.SetNonHistoricalVariable(r_variable, 0.0,
                                                          mrModelPart.Nodes());

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << "Existing variable " << r_variable.Name()
            << " is cleared for selected flags in " << mrModelPart.Name() << ".\n";

        KRATOS_CATCH("");
    }

    void CalculateCellCenterAverage()
    {
        KRATOS_TRY

        const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);
        const Variable<double>& r_input_variable =
            KratosComponents<Variable<double>>::Get(mAveragingInputVariableName);
        const Variable<double>& r_output_variable =
            KratosComponents<Variable<double>>::Get(mAveragingOutputVariableName);

        RansCalculationUtilities rans_calculation_utilities;

        const int number_of_conditions = mrModelPart.NumberOfConditions();
#pragma omp parallel for
        for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
        {
            Condition& r_condition = *(mrModelPart.ConditionsBegin() + i_condition);
            Condition::GeometryType& r_condition_geometry = r_condition.GetGeometry();
            const int number_of_nodes = r_condition_geometry.PointsNumber();

            if (r_condition.Is(r_flag) != mFlagVariableValue)
                continue;

            SetConditionParent(r_condition);
            const Element& r_parent_element = *r_condition.GetValue(PARENT_ELEMENT);
            const Element::GeometryType& r_parent_element_geometry =
                r_parent_element.GetGeometry();

            Vector gauss_weights;
            Matrix shape_functions;
            Element::GeometryType::ShapeFunctionsGradientsType shape_function_derivatives;

            rans_calculation_utilities.CalculateGeometryData(
                r_parent_element_geometry, GeometryData::GI_GAUSS_1,
                gauss_weights, shape_functions, shape_function_derivatives);

            const double gauss_weight = gauss_weights[0];
            const Vector& gauss_shape_functions = row(shape_functions, 0);
            const double value = rans_calculation_utilities.EvaluateInPoint(
                r_parent_element_geometry, r_input_variable, gauss_shape_functions);

            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_node = r_condition_geometry[i_node];

                const double current_output_value = r_node.GetValue(r_output_variable);
                r_node.SetLock();
                r_node.SetValue(r_output_variable, current_output_value + value);
                r_node.UnSetLock();
            }
        }

        const int number_of_nodes = mrModelPart.NumberOfNodes();
        unsigned int number_of_calculated_nodes = 0;
#pragma omp parallel for reduction(+ : number_of_calculated_nodes)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);
            if (r_node.Is(r_flag) == mFlagVariableValue)
            {
                const double number_of_neighbouring_conditions =
                    static_cast<double>(r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS));
                r_node.FastGetSolutionStepValue(r_output_variable) =
                    r_node.GetValue(r_output_variable) / number_of_neighbouring_conditions;
                number_of_calculated_nodes++;
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Cell centered averaging is done for " << number_of_calculated_nodes
            << " nodes in " << mrModelPart.Name() << ".\n";

        KRATOS_CATCH("");
    }

    void FindNodalNeighbours()
    {
        KRATOS_TRY

        FindNodalNeighboursProcess find_nodal_neighbours_process(
            mrModelPart, mAverageNeighbourElements, mAverageNeighbourNodes);
        find_nodal_neighbours_process.Execute();

        const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

        RansVariableUtils().SetNonHistoricalScalarVar(
            NUMBER_OF_NEIGHBOUR_CONDITIONS, 0, mrModelPart.Nodes());

        const int number_of_conditions = mrModelPart.NumberOfConditions();
#pragma omp parallel for
        for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
        {
            Condition& r_condition = *(mrModelPart.ConditionsBegin() + i_condition);

            if (r_condition.Is(r_flag) != mFlagVariableValue)
                continue;

            Condition::GeometryType& r_condition_geometry = r_condition.GetGeometry();
            const int number_of_nodes = r_condition_geometry.PointsNumber();

            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_node = r_condition_geometry[i_node];
                r_node.SetLock();
                const int current_number_of_neighbour_condtions =
                    r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS);
                r_node.SetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS,
                                current_number_of_neighbour_condtions + 1);
                r_node.UnSetLock();
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << "Nodal neighbours found in " << mrModelPart.Name() << ".\n";

        KRATOS_CATCH("");
    }

    void SetConditionsFlag()
    {
        KRATOS_TRY

        const int number_of_conditions = mrModelPart.NumberOfConditions();

        const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

        RansCalculationUtilities rans_calculation_utilities;

#pragma omp parallel for
        for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
        {
            Condition& r_condition = *(mrModelPart.ConditionsBegin() + i_condition);
            Condition::GeometryType& r_condition_geometry = r_condition.GetGeometry();
            const int number_of_nodes = r_condition_geometry.PointsNumber();

            bool condition_flag = true;
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                if (!r_condition_geometry[i_node].Is(r_flag))
                {
                    condition_flag = false;
                    break;
                }
            }
            r_condition.Set(r_flag, condition_flag);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << "Condition flags set in " << mrModelPart.Name() << ".\n";

        KRATOS_CATCH("");
    }

    void SetConditionParent(Condition& rCondition)
    {
        KRATOS_TRY

        GlobalPointersVector<Element> element_candidates;
        const Condition::GeometryType& r_condition_geometry = rCondition.GetGeometry();
        const int number_of_condition_nodes = r_condition_geometry.PointsNumber();

        std::vector<int> node_ids(number_of_condition_nodes), element_node_ids;

        for (int i_node = 0; i_node < number_of_condition_nodes; ++i_node)
        {
            const NodeType& r_node = r_condition_geometry[i_node];
            const GlobalPointersVector<Element>& r_node_element_candidates =
                r_node.GetValue(NEIGHBOUR_ELEMENTS);
            for (int j_element = 0;
                 j_element < static_cast<int>(r_node_element_candidates.size()); ++j_element)
            {
                element_candidates.push_back(r_node_element_candidates(j_element));
            }
            node_ids[i_node] = r_node.Id();
        }

        std::sort(node_ids.begin(), node_ids.end());

        for (int i_element = 0;
             i_element < static_cast<int>(element_candidates.size()); ++i_element)
        {
            const Element::GeometryType& r_geometry =
                element_candidates[i_element].GetGeometry();
            const int number_of_element_candidate_nodes = r_geometry.PointsNumber();
            if (element_node_ids.size() != number_of_element_candidate_nodes)
                element_node_ids.resize(number_of_element_candidate_nodes);

            for (int i_node = 0; i_node < number_of_element_candidate_nodes; ++i_node)
            {
                element_node_ids[i_node] = r_geometry[i_node].Id();
            }

            std::sort(element_node_ids.begin(), element_node_ids.end());
            if (std::includes(element_node_ids.begin(), element_node_ids.end(),
                              node_ids.begin(), node_ids.end()))
            {
                rCondition.SetValue(PARENT_ELEMENT, element_candidates(i_element));
                return;
            }
        }

        KRATOS_ERROR << "Parent element for condition id=" << rCondition.Id() << " not found.\n";
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
    RansScalarCellCenterAveragingProcess& operator=(RansScalarCellCenterAveragingProcess const& rOther);

    /// Copy constructor.
    RansScalarCellCenterAveragingProcess(RansScalarCellCenterAveragingProcess const& rOther);

    ///@}

}; // Class RansScalarCellCenterAveragingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansScalarCellCenterAveragingProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_SCALAR_CELL_CENTER_AVERAGING_PROCESS_H_INCLUDED defined
