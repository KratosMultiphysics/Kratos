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

#include "custom_utilities/rans_check_utilities.h"

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

/**
 * @brief Apply a specific flag for nodes and conditions
 *
 * This process apply a given flag to nodes of the modelpart.
 * Then, if preferred, applies to all conditions in the given model, which has
 * the given flag applied to all the nodes in the specific condition.
 *
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
            "apply_to_model_part_conditions" : ["ALL_MODEL_PARTS"]
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mModelPartName = mrParameters["model_part_name"].GetString();
        mFlagVariableName = mrParameters["flag_variable_name"].GetString();
        mFlagVariableValue = mrParameters["flag_variable_value"].GetBool();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartsForConditionFlags =
            mrParameters["apply_to_model_part_conditions"].GetStringArray();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansApplyFlagProcess() override
    {
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

        RansCheckUtilities().CheckIfModelPartExists(mrModel, mModelPartName);

        return 0;

        KRATOS_CATCH("");
    }

    void ExecuteInitialize() override
    {
        ApplyNodeFlags();

        if (mModelPartsForConditionFlags.size() == 1 &&
            mModelPartsForConditionFlags[0] == "ALL_MODEL_PARTS")
        {
            mModelPartsForConditionFlags.clear();
            for (const std::string model_part_name : mrModel.GetModelPartNames())
            {
                mModelPartsForConditionFlags.push_back(model_part_name);
            }
        }

        for (std::string model_part_name : mModelPartsForConditionFlags)
        {
            ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);
            ApplyConditionFlags(r_model_part);
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
