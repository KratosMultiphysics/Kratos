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

#if !defined(KRATOS_RANS_K_TURBULENT_INTENSITY_EVALUATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_K_TURBULENT_INTENSITY_EVALUATION_PROCESS_H_INCLUDED

// System includes
#include <cmath>
#include <string>

// External includes

// Project includes
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

class RansKTurbulentIntensityEvaluationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    /// Pointer definition of RansKTurbulentIntensityEvaluationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansKTurbulentIntensityEvaluationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansKTurbulentIntensityEvaluationProcess(ModelPart& rModelPart,
                                             Parameters& rParameters,
                                             const bool IsConstrained)
        : mrModelPart(rModelPart), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "turbulent_intensity"       : 0.05
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mTurbulentIntensity = mrParameters["turbulent_intensity"].GetDouble();
        mIsConstrained = IsConstrained;

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansKTurbulentIntensityEvaluationProcess() override
    {
        // delete mpDistanceCalculator;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    void ExecuteInitializeSolutionStep() override
    {
        Execute();
    }

    void Execute() override
    {
        KRATOS_TRY

        mrModelPart.GetProcessInfo()[RANS_PROCESS_WALL_VELOCITY]->Execute();

        KRATOS_INFO(this->Info()) << "Applying k values to " << mrModelPart.Name() << ".\n";
        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            CalculateTurbulentValues(r_node);
            if (mIsConstrained)
                r_node.Fix(TURBULENT_KINETIC_ENERGY);
        }

        KRATOS_CATCH("");
    }

    int Check() override
    {
        // Checking variable definitions
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        }

        return 0;
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
        return std::string("RansKTurbulentIntensityEvaluationProcess");
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

    double mTurbulentIntensity;

    bool mIsConstrained;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateTurbulentValues(NodeType& rNode)
    {
        const array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        double velocity_magnitude = norm_2(r_velocity);

        if (velocity_magnitude < std::numeric_limits<double>::epsilon())
        {
            const array_1d<double, 3>& r_wall_velocity =
                rNode.FastGetSolutionStepValue(WALL_VELOCITY);
            velocity_magnitude = norm_2(r_wall_velocity);
        }

        rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) =
            1.5 * std::pow(mTurbulentIntensity * velocity_magnitude, 2);
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
    RansKTurbulentIntensityEvaluationProcess& operator=(RansKTurbulentIntensityEvaluationProcess const& rOther);

    /// Copy constructor.
    RansKTurbulentIntensityEvaluationProcess(RansKTurbulentIntensityEvaluationProcess const& rOther);

    ///@}

}; // Class RansKTurbulentIntensityEvaluationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansKTurbulentIntensityEvaluationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_K_TURBULENT_INTENSITY_EVALUATION_PROCESS_H_INCLUDED defined
