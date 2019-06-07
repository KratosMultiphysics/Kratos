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

#if !defined(KRATOS_RANS_WALL_VELOCITY_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_WALL_VELOCITY_CALCULATION_PROCESS_H_INCLUDED

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

class RansWallVelocityCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansWallVelocityCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansWallVelocityCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansWallVelocityCalculationProcess(ModelPart& rModelPart, Parameters& rParameters)
        : mrModelPart(rModelPart), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "average_neighbour_elements" : 10,
            "average_neighbour_nodes"    : 10,
            "weight"                     : 1.0,
            "echo_level"                 : 0
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mAverageNeighbourElements = mrParameters["average_neighbour_elements"].GetInt();
        mAverageNeighbourNodes = mrParameters["average_neighbour_nodes"].GetInt();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mWeight = mrParameters["weight"].GetDouble();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansWallVelocityCalculationProcess() override
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
        KRATOS_CHECK_VARIABLE_KEY(WALL_VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(NEIGHBOUR_NODES);

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(WALL_VELOCITY, r_node);
        }

        return 0;
    }

    void Execute() override
    {
        CalculateWallVelocities();
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
        return std::string("RansWallVelocityCalculationProcess");
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

    double mWeight;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateNormals()
    {
        NormalCalculationUtils().CalculateOnSimplex(
            mrModelPart, mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);

        // Check for corner nodes where INLET/OUTLET and STRUCTURE is true
        const int number_of_nodes = mrModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);

            if (r_node.Is(STRUCTURE) && (r_node.Is(INLET) || r_node.Is(OUTLET)))
            {
                const GlobalPointersVector<NodeType>& r_neighbour_nodes =
                    r_node.GetValue(NEIGHBOUR_NODES);
                const int number_of_neighbour_nodes = r_neighbour_nodes.size();

                double number_of_aggregation_nodes = 0.0;
                array_1d<double, 3> averaged_normal;
                averaged_normal.clear();
                for (int j_node = 0; j_node < number_of_neighbour_nodes; ++j_node)
                {
                    const NodeType& r_neighbour_node = r_neighbour_nodes[j_node];
                    if (r_neighbour_node.Is(STRUCTURE) &&
                        (!r_neighbour_node.Is(INLET) && !r_neighbour_node.Is(OUTLET)))
                    {
                        averaged_normal +=
                            r_neighbour_node.FastGetSolutionStepValue(NORMAL);
                        number_of_aggregation_nodes += 1.0;
                    }
                }

                if (number_of_aggregation_nodes > 0.0)
                {
                    r_node.FastGetSolutionStepValue(NORMAL) =
                        averaged_normal * (1.0 / number_of_aggregation_nodes);
                }
            }
        }
    }

    void CalculateWallVelocities()
    {
        KRATOS_TRY

        FindNodalNeighbours();
        CalculateNormals();

        const int number_of_nodes = mrModelPart.NumberOfNodes();

        unsigned int number_of_velocity_based_node_calculations = 0;
        unsigned int number_of_model_based_node_calculations = 0;

#pragma omp parallel for reduction(+: number_of_velocity_based_node_calculations, number_of_model_based_node_calculations)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);

            if (r_node.Is(STRUCTURE))
            {
                const array_1d<double, 3>& r_normal =
                    r_node.FastGetSolutionStepValue(NORMAL);
                const double normal_magnitude = norm_2(r_normal);

                KRATOS_ERROR_IF(normal_magnitude < std::numeric_limits<double>::epsilon())
                    << "NORMAL magnitude is zero in node " << r_node.Id() << " at "
                    << r_node.Coordinates() << " in " << mrModelPart.Name()
                    << " [ NORMAL = " << r_normal << " ].\n";

                const array_1d<double, 3> r_unit_normal = r_normal / norm_2(r_normal);

                if (r_node.Is(SLIP))
                {
                    const array_1d<double, 3>& r_velocity =
                        r_node.FastGetSolutionStepValue(VELOCITY);
                    r_node.FastGetSolutionStepValue(WALL_VELOCITY) =
                        r_velocity - r_unit_normal * inner_prod(r_velocity, r_unit_normal);
                    number_of_velocity_based_node_calculations++;
                }
                else
                {
                    const GlobalPointersVector<NodeType>& r_neighbour_nodes =
                        r_node.GetValue(NEIGHBOUR_NODES);
                    const int number_of_neighbour_nodes = r_neighbour_nodes.size();
                    array_1d<double, 3> wall_velocity;
                    wall_velocity.clear();
                    double number_of_aggregation_nodes = 0.0;
                    for (int j_node = 0; j_node < number_of_neighbour_nodes; ++j_node)
                    {
                        const NodeType& r_neighbour_node = r_neighbour_nodes[j_node];
                        if (!r_neighbour_node.Is(STRUCTURE))
                        {
                            const array_1d<double, 3>& r_velocity =
                                r_neighbour_node.FastGetSolutionStepValue(VELOCITY);
                            wall_velocity +=
                                (r_velocity - r_unit_normal * inner_prod(r_velocity, r_unit_normal));
                            number_of_aggregation_nodes += 1.0;
                        }
                    }

                    if (number_of_aggregation_nodes > 0.0)
                    {
                        wall_velocity *= (mWeight / number_of_aggregation_nodes);
                        r_node.FastGetSolutionStepValue(WALL_VELOCITY) = wall_velocity;
                        number_of_model_based_node_calculations++;
                    }
                }
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Wall velocities calculated in " << mrModelPart.Name()
            << " using VELOCITY / MODEL [ " << number_of_velocity_based_node_calculations
            << " / " << number_of_model_based_node_calculations << " ] nodes.\n";

        KRATOS_CATCH("");
    }

    void FindNodalNeighbours()
    {
        FindNodalNeighboursProcess find_nodal_neighbours_process(
            mrModelPart, mAverageNeighbourElements, mAverageNeighbourNodes);
        find_nodal_neighbours_process.Execute();
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
    RansWallVelocityCalculationProcess& operator=(RansWallVelocityCalculationProcess const& rOther);

    /// Copy constructor.
    RansWallVelocityCalculationProcess(RansWallVelocityCalculationProcess const& rOther);

    ///@}

}; // Class RansWallVelocityCalculationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansWallVelocityCalculationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_WALL_VELOCITY_CALCULATION_PROCESS_H_INCLUDED defined
