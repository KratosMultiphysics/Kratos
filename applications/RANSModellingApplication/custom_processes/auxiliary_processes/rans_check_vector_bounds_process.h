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

#if !defined(KRATOS_RANS_CHECK_VECTOR_BOUNDS_PROCESS_H_INCLUDED)
#define KRATOS_RANS_CHECK_VECTOR_BOUNDS_PROCESS_H_INCLUDED

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

class RansCheckVectorBoundsProcess : public Process
{
    enum VectorComponent
    {
        Magnitude,
        X,
        Y,
        Z
    };

public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansCheckVectorBoundsProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansCheckVectorBoundsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansCheckVectorBoundsProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"      : 0,
            "component_type"       : "magnitude"
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = mrParameters["variable_name"].GetString();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        if (mrParameters["component_type"].GetString() == "magnitude")
            mVectorComponent = VectorComponent::Magnitude;
        else if (mrParameters["component_type"].GetString() == "x")
            mVectorComponent = VectorComponent::X;
        else if (mrParameters["component_type"].GetString() == "y")
            mVectorComponent = VectorComponent::Y;
        else if (mrParameters["component_type"].GetString() == "z")
            mVectorComponent = VectorComponent::Z;
        else
            KRATOS_ERROR
                << "Vector component_type type not found. [ component_type = "
                << mrParameters["component_type"].GetString() << " ]\n.";

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansCheckVectorBoundsProcess() override
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

        const Variable<array_1d<double, 3>> vector_variable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

        KRATOS_CHECK_VARIABLE_KEY(vector_variable);

        ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(vector_variable, r_node);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Check passed for " << mModelPartName << " with variable "
            << vector_variable.Name() << ".\n";

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        KRATOS_TRY
        const Variable<array_1d<double, 3>> vector_variable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

        const ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();

        array_1d<double, 3> vector_weights;

        switch (mVectorComponent)
        {
        case VectorComponent::Magnitude:
            vector_weights[0] = 1.0;
            vector_weights[1] = 1.0;
            vector_weights[2] = 1.0;
            break;
        case VectorComponent::X:
            vector_weights[0] = 1.0;
            vector_weights[1] = 0.0;
            vector_weights[2] = 0.0;
            break;
        case VectorComponent::Y:
            vector_weights[0] = 0.0;
            vector_weights[1] = 1.0;
            vector_weights[2] = 0.0;
            break;
        case VectorComponent::Z:
            vector_weights[0] = 0.0;
            vector_weights[1] = 0.0;
            vector_weights[2] = 1.0;
            break;
        }

        double min_value = 0.0;
        double max_value = 0.0;
        bool is_initialized = false;

        const int number_of_nodes = r_nodes.size();
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const ModelPart::NodeType& r_node = *(r_nodes.begin() + i_node);
            const array_1d<double, 3>& vector_value =
                r_node.FastGetSolutionStepValue(vector_variable);

            double current_value = 0.0;
            for (int dim = 0; dim < 3; ++dim)
                current_value += std::pow(vector_value[dim] * vector_weights[dim], 2);

            current_value = std::pow(current_value, 0.5);

            if (!is_initialized)
            {
                min_value = current_value;
                max_value = current_value;
                is_initialized = true;
            }

            if (min_value > current_value)
                min_value = current_value;
            if (max_value < current_value)
                max_value = current_value;
        }

        KRATOS_INFO(this->Info())
            << vector_variable.Name() << " is bounded between [ " << min_value
            << ", " << max_value << " ] in " << mModelPartName << ".\n";

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
        return std::string("RansCheckVectorBoundsProcess");
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

    VectorComponent mVectorComponent;

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
    RansCheckVectorBoundsProcess& operator=(RansCheckVectorBoundsProcess const& rOther);

    /// Copy constructor.
    RansCheckVectorBoundsProcess(RansCheckVectorBoundsProcess const& rOther);

    ///@}

}; // Class RansCheckVectorBoundsProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansCheckVectorBoundsProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_CHECK_VECTOR_BOUNDS_PROCESS_H_INCLUDED defined
