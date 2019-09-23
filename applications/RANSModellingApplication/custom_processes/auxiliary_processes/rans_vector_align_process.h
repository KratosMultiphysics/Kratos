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

#if !defined(KRATOS_RANS_VECTOR_ALIGN_PROCESS_H_INCLUDED)
#define KRATOS_RANS_VECTOR_ALIGN_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "custom_utilities/rans_variable_utils.h"
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

class RansVectorAlignProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansVectorAlignProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansVectorAlignProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansVectorAlignProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_variable_name"     : "PLEASE_SPECIFY_INPUT_VARIABLE_NAME",
            "output_variable_name"    : "PLEASE_SPECIFY_OUTPUT_VARIABLE_NAME",
            "alignment_variable_name" : "NORMAL",
            "is_tangential_alignment" : true,
            "echo_level"              : 0
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mModelPartName = mrParameters["model_part_name"].GetString();
        mInputVariableName = mrParameters["input_variable_name"].GetString();
        mOutputVariableName = mrParameters["output_variable_name"].GetString();
        mAlignmentVariableName = mrParameters["alignment_variable_name"].GetString();
        mIsTangentialAlignment = mrParameters["is_tangential_alignment"].GetBool();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansVectorAlignProcess() override
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
        return 0;
    }

    void Execute() override
    {
        CalculateAlignmentVectors();
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
        return std::string("RansVectorAlignProcess");
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
    std::string mAlignmentVariableName;
    std::string mInputVariableName;
    std::string mOutputVariableName;

    bool mIsTangentialAlignment;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateAlignmentVectors()
    {
        KRATOS_TRY

        const Variable<array_1d<double, 3>>& r_input_variable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mInputVariableName);
        const Variable<array_1d<double, 3>>& r_output_variable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mOutputVariableName);
        const Variable<array_1d<double, 3>>& r_alignment_variable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mAlignmentVariableName);

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const int number_of_nodes = r_model_part.NumberOfNodes();

        if (mIsTangentialAlignment)
        {
#pragma omp parallel for
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_node = *(r_model_part.NodesBegin() + i_node);

                const array_1d<double, 3>& r_normal =
                    r_node.FastGetSolutionStepValue(r_alignment_variable);
                const double normal_magnitude = norm_2(r_normal);

                KRATOS_ERROR_IF(normal_magnitude < std::numeric_limits<double>::epsilon())
                    << "NORMAL magnitude is zero in node " << r_node.Id()
                    << " at " << r_node.Coordinates() << " in "
                    << mModelPartName << " [ NORMAL = " << r_normal << " ].\n";

                const array_1d<double, 3> r_unit_normal = r_normal / norm_2(r_normal);

                const array_1d<double, 3>& r_velocity =
                    r_node.FastGetSolutionStepValue(r_input_variable);
                r_node.FastGetSolutionStepValue(r_output_variable) =
                    r_velocity - r_unit_normal * inner_prod(r_velocity, r_unit_normal);
            }
        }
        else
        {
#pragma omp parallel for
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_node = *(r_model_part.NodesBegin() + i_node);

                const array_1d<double, 3>& r_normal =
                    r_node.FastGetSolutionStepValue(r_alignment_variable);
                const double normal_magnitude = norm_2(r_normal);

                KRATOS_ERROR_IF(normal_magnitude < std::numeric_limits<double>::epsilon())
                    << "NORMAL magnitude is zero in node " << r_node.Id()
                    << " at " << r_node.Coordinates() << " in "
                    << mModelPartName << " [ NORMAL = " << r_normal << " ].\n";

                const array_1d<double, 3> r_unit_normal = r_normal / norm_2(r_normal);

                const array_1d<double, 3>& r_velocity =
                    r_node.FastGetSolutionStepValue(r_input_variable);
                r_node.FastGetSolutionStepValue(r_output_variable) =
                    r_unit_normal * inner_prod(r_velocity, r_unit_normal);
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << r_input_variable.Name() << " is aligned using "
            << r_alignment_variable.Name() << " and stored at "
            << r_output_variable.Name() << " for nodes in " << mModelPartName << ".\n";

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
    RansVectorAlignProcess& operator=(RansVectorAlignProcess const& rOther);

    /// Copy constructor.
    RansVectorAlignProcess(RansVectorAlignProcess const& rOther);

    ///@}

}; // Class RansVectorAlignProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const RansVectorAlignProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_VECTOR_ALIGN_PROCESS_H_INCLUDED defined
