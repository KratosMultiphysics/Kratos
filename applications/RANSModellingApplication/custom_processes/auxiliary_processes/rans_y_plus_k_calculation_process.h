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

#if !defined(KRATOS_RANS_Y_PLUS_K_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_Y_PLUS_K_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "custom_utilities/rans_calculation_utilities.h"
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

class RansYPlusKCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansYPlusKCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansYPlusKCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansYPlusKCalculationProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09,
            "min_y_plus"      : 12.0
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mModelPartName = mrParameters["model_part_name"].GetString();

        mEchoLevel = mrParameters["echo_level"].GetInt();
        mMinYPlus = mrParameters["min_y_plus"].GetDouble();
        mCmu25 = std::pow(mrParameters["c_mu"].GetDouble(), 0.25);

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansYPlusKCalculationProcess() override
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
            CalculateWallLawYplus(r_node, VELOCITY);
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
        return std::string("RansYPlusKCalculationProcess");
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
    double mMinYPlus;
    double mCmu25;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateWallLawYplus(NodeType& rNode, const Variable<array_1d<double, 3>>& rVelocityVariable)
    {
        KRATOS_TRY

        const double nu = rNode.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double tke = rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double wall_distance = rNode.FastGetSolutionStepValue(DISTANCE);

        const double u_tau = mCmu25 * std::sqrt(std::max(tke, 0.0));
        rNode.FastGetSolutionStepValue(RANS_Y_PLUS) =
            std::max(u_tau * wall_distance / nu, mMinYPlus);

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
    RansYPlusKCalculationProcess& operator=(RansYPlusKCalculationProcess const& rOther);

    /// Copy constructor.
    RansYPlusKCalculationProcess(RansYPlusKCalculationProcess const& rOther);

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
                                const RansYPlusKCalculationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_Y_PLUS_K_CALCULATION_PROCESS_H_INCLUDED  defined
