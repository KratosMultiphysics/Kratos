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

#if !defined(KRATOS_RANS_APPLY_FLAG_PROCESS_H_INCLUDED)
#define KRATOS_RANS_APPLY_FLAG_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

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

class RansApplyFlagProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansApplyFlagProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansApplyFlagProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansApplyFlagProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"                : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"                     : 0,
            "flag_variable_name"             : "PLEASE_SPECIFY_FLAG_VARIABLE_NAME",
            "flag_variable_value"            : true,
            "apply_to_model_part_conditions" : "all"
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mModelPartName = mrParameters["model_part_name"].GetString();
        mFlagVariableName = mrParameters["flag_variable_name"].GetString();
        mFlagVariableValue = mrParameters["flag_variable_value"].GetBool();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        if (mrParameters["apply_to_model_part_conditions"].GetString() == "all")
        {
            mIsAllModelPartsConsideredForConditionFlags = true;
        }
        else
        {
            mIsAllModelPartsConsideredForConditionFlags = false;
        }

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansApplyFlagProcess() override
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

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << "Check passed for " << mModelPartName << ".\n";

        return 0;

        KRATOS_CATCH("");
    }

    void ExecuteInitialize() override
    {
        ApplyNodeFlags();

        if (mIsAllModelPartsConsideredForConditionFlags)
        {
            const std::vector<std::string>& r_model_part_names =
                mrModel.GetModelPartNames();
            for (std::string model_part_name : r_model_part_names)
            {
                ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);
                ApplyConditionFlags(r_model_part);
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << mFlagVariableName << " condition flag set for " << mModelPartName << ".\n";
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
        return std::string("RansApplyFlagProcess");
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
    int mEchoLevel;

    std::string mModelPartName;

    std::string mFlagVariableName;
    bool mFlagVariableValue;
    bool mIsAllModelPartsConsideredForConditionFlags;
    std::vector<std::string> mModelPartsForConditionFlags;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ApplyNodeFlags()
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
        const int number_of_nodes = r_model_part.NumberOfNodes();

        const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            (r_model_part.NodesBegin() + i_node)->Set(r_flag, mFlagVariableValue);

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << mFlagVariableName << " is set to nodes " << mFlagVariableValue
            << " in " << mModelPartName << ".\n";

        KRATOS_CATCH("");
    }

    void ApplyConditionFlags(ModelPart& rModelPart)
    {
        KRATOS_TRY

        const int number_of_conditions = rModelPart.NumberOfConditions();

        const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

#pragma omp parallel for
        for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
        {
            Condition& r_condition = *(rModelPart.ConditionsBegin() + i_condition);
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
            << mFlagVariableName << " flags set for conditions in "
            << rModelPart.Name() << ".\n";

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
    RansApplyFlagProcess& operator=(RansApplyFlagProcess const& rOther);

    /// Copy constructor.
    RansApplyFlagProcess(RansApplyFlagProcess const& rOther);

    ///@}

}; // Class RansApplyFlagProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const RansApplyFlagProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_APPLY_FLAG_PROCESS_H_INCLUDED defined
