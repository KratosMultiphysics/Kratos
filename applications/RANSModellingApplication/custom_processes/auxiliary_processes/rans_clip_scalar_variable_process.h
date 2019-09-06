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

#if !defined(KRATOS_RANS_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED)
#define KRATOS_RANS_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED

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
#include "custom_utilities/rans_variable_utils.h"

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
 * @brief Clips given scalar variable to a range
 *
 * This process clips a given scalar variable to a range in all nodes in the model part.
 *
 */

class RansClipScalarVariableProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansClipScalarVariableProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansClipScalarVariableProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansClipScalarVariableProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"      : 0,
            "min_value"       : 1e-18,
            "max_value"       : 1e+30
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = mrParameters["variable_name"].GetString();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mEchoLevel = mrParameters["echo_level"].GetInt();
        mMinValue = mrParameters["min_value"].GetDouble();
        mMaxValue = mrParameters["max_value"].GetDouble();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansClipScalarVariableProcess() override
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

        const Variable<double>& r_scalar_variable =
            KratosComponents<Variable<double>>::Get(mVariableName);

        RansCheckUtilities rans_check_utilities;

        rans_check_utilities.CheckIfModelPartExists(mrModel, mModelPartName);
        rans_check_utilities.CheckIfVariableExistsInModelPart(
            mrModel.GetModelPart(mModelPartName), r_scalar_variable);

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        KRATOS_TRY

        const Variable<double>& r_scalar_variable =
            KratosComponents<Variable<double>>::Get(mVariableName);

        RansVariableUtils rans_variable_utils;

        unsigned int nodes_below, nodes_above, total_nodes;

        rans_variable_utils.ClipScalarVariable(
            nodes_below, nodes_above, total_nodes, mMinValue, mMaxValue,
            r_scalar_variable, mrModel.GetModelPart(mModelPartName));

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0 && (nodes_below > 0 || nodes_above > 0))
            << mVariableName << " is clipped between [ " << mMinValue << ", "
            << mMaxValue << " ]. [ " << nodes_below << " nodes < " << mMinValue
            << " and " << nodes_above << " nodes > " << mMaxValue << " out of "
            << total_nodes << " total nodes in " << mModelPartName << " ].\n";

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
        return std::string("RansClipScalarVariableProcess");
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
    RansClipScalarVariableProcess& operator=(RansClipScalarVariableProcess const& rOther);

    /// Copy constructor.
    RansClipScalarVariableProcess(RansClipScalarVariableProcess const& rOther);

    ///@}

}; // Class RansClipScalarVariableProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansClipScalarVariableProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_CLIP_SCALAR_VARIABLE_PROCESS_H_INCLUDED defined
