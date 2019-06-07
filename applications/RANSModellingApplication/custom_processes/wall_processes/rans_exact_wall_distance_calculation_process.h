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

#if !defined(KRATOS_RANS_EXACT_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_EXACT_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED

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
#include "processes/variational_distance_calculation_process.h"
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

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class RansExactWallDistanceCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansExactWallDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansExactWallDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansExactWallDistanceCalculationProcess(ModelPart& rModelPart, Parameters& rParameters)
        : mrModelPart(rModelPart), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "max_iterations"             : 10,
            "average_neighbour_elements" : 10,
            "average_neighbour_nodes"    : 10,
            "weight"                     : 1.0,
            "echo_level"                 : 0,
            "linear_solver_settings"     : {
                "solver_type"     : "amgcl"
            }
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(
            mrParameters["linear_solver_settings"]);

        mMaxIterations = mrParameters["max_iterations"].GetInt();
        mAverageNeighbourElements = mrParameters["average_neighbour_elements"].GetInt();
        mAverageNeighbourNodes = mrParameters["average_neighbour_nodes"].GetInt();
        mWeight = mrParameters["weight"].GetDouble();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansExactWallDistanceCalculationProcess() override
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
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(Y_WALL);
        KRATOS_CHECK_VARIABLE_KEY(NEIGHBOUR_NODES);

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
        }

        return 0;
    }

    void Execute() override
    {
        CalculateWallDistances();
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
        return std::string("RansExactWallDistanceCalculationProcess");
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

    typename TLinearSolver::Pointer mpLinearSolver;

    int mMaxIterations;
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

    void InitializeWallNodes()
    {
        RansVariableUtils().SetScalarVar(DISTANCE, 1.0, mrModelPart.Nodes());
        RansVariableUtils().SetScalarVarForFlag(DISTANCE, 0.0, mrModelPart.Nodes(), STRUCTURE);
    }

    void CalculateWallDistances()
    {
        KRATOS_TRY

        InitializeWallNodes();

        VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> distance_calculation_process(
            mrModelPart, mpLinearSolver, mMaxIterations);

        distance_calculation_process.Execute();

        FindNodalNeighbours();

        const int number_of_nodes = mrModelPart.NumberOfNodes();
        int number_of_model_based_node_calculations = 0;

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);
            if (r_node.Is(STRUCTURE))
            {
                const GlobalPointersVector<NodeType>& r_neighbour_nodes =
                    r_node.GetValue(NEIGHBOUR_NODES);
                const int number_of_neighbour_nodes = r_neighbour_nodes.size();
                double distance = 0.0;
                for (int j_node = 0; j_node < number_of_neighbour_nodes; ++j_node)
                {
                    const NodeType& r_neighbour_node = r_neighbour_nodes[j_node];
                    distance += r_neighbour_node.FastGetSolutionStepValue(DISTANCE);
                }

                distance *= (mWeight / static_cast<double>(number_of_neighbour_nodes));
                r_node.SetValue(Y_WALL, distance);
            }
        }

// TODO: Remove this, after migrating all the calculations in the elements from DISTANCE to Y_WALL
#pragma omp parallel for reduction(+ : number_of_model_based_node_calculations)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(mrModelPart.NodesBegin() + i_node);
            if (r_node.Is(STRUCTURE))
            {
                r_node.FastGetSolutionStepValue(DISTANCE) = r_node.GetValue(Y_WALL);
                number_of_model_based_node_calculations++;
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Wall distances calculated in " << mrModelPart.Name() << " using DISTANCE / MODEL [ "
            << number_of_nodes - number_of_model_based_node_calculations
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
    RansExactWallDistanceCalculationProcess& operator=(RansExactWallDistanceCalculationProcess const& rOther);

    /// Copy constructor.
    RansExactWallDistanceCalculationProcess(RansExactWallDistanceCalculationProcess const& rOther);

    ///@}

}; // Class RansExactWallDistanceCalculationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansExactWallDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EXACT_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED defined
