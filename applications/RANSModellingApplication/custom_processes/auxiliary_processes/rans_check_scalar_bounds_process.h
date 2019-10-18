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

#if !defined(KRATOS_RANS_CHECK_SCALAR_BOUNDS_PROCESS_H_INCLUDED)
#define KRATOS_RANS_CHECK_SCALAR_BOUNDS_PROCESS_H_INCLUDED

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
 * @brief Checks lower and upper bounds of a scalar
 *
 * This process checks lower and upper bounds of a variable in nodes in a given modelpart
 *
 */

class RansCheckScalarBoundsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansCheckScalarBoundsProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansCheckScalarBoundsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansCheckScalarBoundsProcess(Model& rModel, Parameters& rParameters)
        : Process(), mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE"
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = mrParameters["variable_name"].GetString();
        mModelPartName = mrParameters["model_part_name"].GetString();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansCheckScalarBoundsProcess() override
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

        const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        RansVariableUtils rans_variable_utils;
        const double min_value =
            rans_variable_utils.GetMinimumScalarValue(r_model_part, r_scalar_variable);
        const double max_value =
            rans_variable_utils.GetMaximumScalarValue(r_model_part, r_scalar_variable);

        KRATOS_INFO(this->Info())
            << r_scalar_variable.Name() << " is bounded between [ " << min_value
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
        return std::string("RansCheckScalarBoundsProcess");
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
    RansCheckScalarBoundsProcess& operator=(RansCheckScalarBoundsProcess const& rOther);

    /// Copy constructor.
    RansCheckScalarBoundsProcess(RansCheckScalarBoundsProcess const& rOther);

    ///@}

}; // Class RansCheckScalarBoundsProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansCheckScalarBoundsProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_CHECK_SCALAR_BOUNDS_PROCESS_H_INCLUDED defined
