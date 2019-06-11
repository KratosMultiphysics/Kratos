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

#if !defined(KRATOS_RANS_SCALAR_NEIGHBOUR_AVERAGING_PROCESS_H_INCLUDED)
#define KRATOS_RANS_SCALAR_NEIGHBOUR_AVERAGING_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
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

class RansScalarNeighbourAveragingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansScalarNeighbourAveragingProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansScalarNeighbourAveragingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansScalarNeighbourAveragingProcess(ModelPart& rModelPart, Parameters& rParameters)
        : mrModelPart(rModelPart), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "echo_level"    : 0,
            "variable_name" : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "weight"        : 1.0
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = mrParameters["echo_level"].GetInt();

        mWeight = mrParameters["weight"].GetDouble();
        mVariableName = mrParameters["variable_name"].GetString();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansScalarNeighbourAveragingProcess() override
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
        KRATOS_CHECK_VARIABLE_KEY(NEIGHBOUR_NODES);

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(scalar_variable, r_node);
        }

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        KRATOS_TRY

        const Variable<double> scalar_variable =
            KratosComponents<Variable<double>>::Get(mVariableName);

        const int number_of_nodes = mrModelPart.NumberOfNodes();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);
            r_node.FastGetSolutionStepValue(scalar_variable) = 0.0;
        }

        int number_of_modified_nodes = 0;
#pragma omp parallel for reduction(+ : number_of_modified_nodes)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);
            const GlobalPointersVector<NodeType>& r_neighbour_nodes =
                r_node.GetValue(NEIGHBOUR_NODES);
            const int number_of_neighbour_nodes = r_neighbour_nodes.size();

            double value = 0.0;
            double number_of_aggregation_nodes = 0.0;
            for (int j_node = 0; j_node < number_of_neighbour_nodes; ++j_node)
            {
                const NodeType& r_neighbour_node = r_neighbour_nodes[j_node];
                value += r_neighbour_node.FastGetSolutionStepValue(scalar_variable);
                number_of_aggregation_nodes += 1.0;
            }

            if (number_of_aggregation_nodes > 0.0)
            {
                r_node.SetValue(scalar_variable, value * mWeight / number_of_aggregation_nodes);
                number_of_modified_nodes++;
            }
        }

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);
            r_node.FastGetSolutionStepValue(scalar_variable) =
                r_node.GetValue(scalar_variable);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << scalar_variable.Name() << " nodal neighbour averaging  done for "
            << number_of_modified_nodes << " nodes out of total "
            << mrModelPart.NumberOfNodes() << " nodes in " << mrModelPart.Name() << "\n";

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
        return std::string("RansScalarNeighbourAveragingProcess");
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

    double mWeight;
    std::string mVariableName;

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
    RansScalarNeighbourAveragingProcess& operator=(RansScalarNeighbourAveragingProcess const& rOther);

    /// Copy constructor.
    RansScalarNeighbourAveragingProcess(RansScalarNeighbourAveragingProcess const& rOther);

    ///@}

}; // Class RansScalarNeighbourAveragingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansScalarNeighbourAveragingProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_SCALAR_NEIGHBOUR_AVERAGING_PROCESS_H_INCLUDED defined
