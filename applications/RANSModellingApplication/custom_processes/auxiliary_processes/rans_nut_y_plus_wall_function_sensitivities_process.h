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

#if !defined(KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_SENSITIVITIES_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_SENSITIVITIES_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "rans_modelling_application_variables.h"
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

class RansNutYPlusWallFunctionSensitivitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansNutYPlusWallFunctionSensitivitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutYPlusWallFunctionSensitivitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutYPlusWallFunctionSensitivitiesProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansNutYPlusWallFunctionSensitivitiesProcess() override
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

        RansCheckUtilities().CheckIfModelPartExists(mrModel, mModelPartName);

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const int number_of_nodes = r_model_part.NumberOfNodes();

        unsigned int number_of_modified_nu_t_nodes = 0;

#pragma omp parallel for reduction(+ : number_of_modified_nu_t_nodes)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            Vector nut_partial_derivatives(2);
            nut_partial_derivatives.clear();
            r_node.SetValue(RANS_NUT_SCALAR_PARTIAL_DERIVATIVES, nut_partial_derivatives);
            number_of_modified_nu_t_nodes++;
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Applied nu_t y_plus wall function sensitivities to "
            << number_of_modified_nu_t_nodes << " of total "
            << r_model_part.NumberOfNodes() << " nodes in " << mModelPartName << "\n";

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
        return std::string("RansNutYPlusWallFunctionSensitivitiesProcess");
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

    int mEchoLevel;

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
    RansNutYPlusWallFunctionSensitivitiesProcess& operator=(
        RansNutYPlusWallFunctionSensitivitiesProcess const& rOther);

    /// Copy constructor.
    RansNutYPlusWallFunctionSensitivitiesProcess(RansNutYPlusWallFunctionSensitivitiesProcess const& rOther);

    ///@}

}; // Class RansNutYPlusWallFunctionSensitivitiesProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansNutYPlusWallFunctionSensitivitiesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_SENSITIVITIES_PROCESS_H_INCLUDED defined
