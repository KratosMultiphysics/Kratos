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

#if !defined(KRATOS_RANS_CLIP_SCALAR_VARIABLE_BY_NEIGHBOUR_AVERAGING_PROCESS_H_INCLUDED)
#define KRATOS_RANS_CLIP_SCALAR_VARIABLE_BY_NEIGHBOUR_AVERAGING_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "containers/model.h"
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

class RansClipScalarVariableByNeighbourAveragingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansClipScalarVariableByNeighbourAveragingProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansClipScalarVariableByNeighbourAveragingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansClipScalarVariableByNeighbourAveragingProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"      : 0,
            "min_value"       : 1e-18,
            "max_value"       : 1e+30,
            "time_step"       : 0
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = mrParameters["variable_name"].GetString();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mMinValue = mrParameters["min_value"].GetDouble();
        mMaxValue = mrParameters["max_value"].GetDouble();
        mTimeStep = mrParameters["time_step"].GetInt();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansClipScalarVariableByNeighbourAveragingProcess() override
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

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
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

        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();

        RansVariableUtils rans_variable_utils;
        rans_variable_utils.SetNonHistoricalVariableToZero(scalar_variable, r_nodes);

        const int number_of_nodes = r_nodes.size();
        int number_of_nodes_below{0}, number_of_nodes_above{0};

#pragma omp parallel for reduction(+:number_of_nodes_below, number_of_nodes_above)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            const double value = r_node.FastGetSolutionStepValue(scalar_variable);

            if (value < mMinValue)
            {
                const GlobalPointersVector<NodeType>& r_neighbour_nodes =
                    r_node.GetValue(NEIGHBOUR_NODES);
                const int number_of_neighbour_nodes = r_neighbour_nodes.size();
                double number_of_nodes_aggregated = 0.0;
                double aggregated_value = 0.0;
                for (int j_node = 0; j_node < number_of_neighbour_nodes; ++j_node)
                {
                    const double neighbour_value =
                        r_neighbour_nodes[j_node].FastGetSolutionStepValue(scalar_variable, mTimeStep);
                    if (neighbour_value > mMinValue)
                    {
                        aggregated_value += neighbour_value;
                        number_of_nodes_aggregated += 1.0;
                    }
                }
                const double current_value =
                    std::max(mMinValue, (number_of_nodes_aggregated > 0.0)
                                            ? aggregated_value / number_of_nodes_aggregated
                                            : 0.0);
                number_of_nodes_below++;
                r_node.SetValue(scalar_variable, current_value);
            }
            else if (value > mMaxValue)
            {
                const GlobalPointersVector<NodeType>& r_neighbour_nodes =
                    r_node.GetValue(NEIGHBOUR_NODES);
                const int number_of_neighbour_nodes = r_neighbour_nodes.size();
                double number_of_nodes_aggregated = 0.0;
                double aggregated_value = 0.0;
                for (int j_node = 0; j_node < number_of_neighbour_nodes; ++j_node)
                {
                    const double neighbour_value =
                        r_neighbour_nodes[j_node].FastGetSolutionStepValue(scalar_variable, mTimeStep);
                    if (neighbour_value < mMaxValue)
                    {
                        aggregated_value += neighbour_value;
                        number_of_nodes_aggregated += 1.0;
                    }
                }
                const double current_value =
                    std::min(mMaxValue, (number_of_nodes_aggregated > 0.0)
                                            ? aggregated_value / number_of_nodes_aggregated
                                            : 0.0);
                number_of_nodes_above++;
                r_node.SetValue(scalar_variable, current_value);
            }
            else
            {
                r_node.SetValue(scalar_variable, value);
            }
        }

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            r_node.FastGetSolutionStepValue(scalar_variable) =
                r_node.GetValue(scalar_variable);
        }

        const double min_value =
            rans_variable_utils.GetMinimumScalarValue(r_nodes, scalar_variable);
        const double max_value =
            rans_variable_utils.GetMaximumScalarValue(r_nodes, scalar_variable);

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0 && (number_of_nodes_below > 0 || number_of_nodes_above > 0))
            << mVariableName << " is clipped between [ " << min_value << ", " << max_value
            << " ]. [ " << number_of_nodes_below << " nodes < " << mMinValue << " and "
            << number_of_nodes_above << " nodes > " << mMaxValue << " out of "
            << r_nodes.size() << " total nodes in " << mModelPartName << " ].\n";

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
        return std::string("RansClipScalarVariableByNeighbourAveragingProcess");
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
    int mTimeStep;

    double mMinValue;
    double mMaxValue;

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
    RansClipScalarVariableByNeighbourAveragingProcess& operator=(
        RansClipScalarVariableByNeighbourAveragingProcess const& rOther);

    /// Copy constructor.
    RansClipScalarVariableByNeighbourAveragingProcess(
        RansClipScalarVariableByNeighbourAveragingProcess const& rOther);

    ///@}

}; // Class RansClipScalarVariableByNeighbourAveragingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansClipScalarVariableByNeighbourAveragingProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_CLIP_SCALAR_VARIABLE_BY_NEIGHBOUR_AVERAGING_PROCESS_H_INCLUDED defined
