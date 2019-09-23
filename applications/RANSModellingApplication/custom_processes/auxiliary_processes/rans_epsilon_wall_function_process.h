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

#if !defined(KRATOS_RANS_EPSILON_WALL_FUNCTION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_EPSILON_WALL_FUNCTION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
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


class RansEpsilonWallFunctionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansEpsilonWallFunctionProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansEpsilonWallFunctionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansEpsilonWallFunctionProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09,
            "von_karman"      : 0.41,
            "is_fixed"        : true
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mCmu = mrParameters["c_mu"].GetDouble();
        mVonKarman = mrParameters["von_karman"].GetDouble();
        mIsConstrained = mrParameters["is_fixed"].GetBool();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansEpsilonWallFunctionProcess() override
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

        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_KINETIC_ENERGY);
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_ENERGY_DISSIPATION_RATE);

        const ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        }

        return 0;

        KRATOS_CATCH("");
    }

    void ExecuteInitialize() override
    {
        if (mIsConstrained)
        {
            ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

            const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
                r_node.Fix(TURBULENT_ENERGY_DISSIPATION_RATE);
            }

            KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
                << "Fixed TURBULENT_ENERGY_DISSIPATION_RATE dofs in "
                << mModelPartName << ".\n";
        }
    }

    void Execute() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const int number_of_nodes = r_model_part.NumberOfNodes();

        const double c_mu_75 = std::pow(mCmu, 0.75);
        const double inv_von_karman = 1.0 / mVonKarman;

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double wall_distance = r_node.FastGetSolutionStepValue(DISTANCE);
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) =
                c_mu_75 * inv_von_karman * std::pow(std::max(tke, 0.0), 1.5) / wall_distance;
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Applied epsilon wall function to " << mModelPartName << "\n";

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
        return std::string("RansEpsilonWallFunctionProcess");
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

    double mCmu;
    double mVonKarman;
    bool mIsConstrained;

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
    RansEpsilonWallFunctionProcess& operator=(RansEpsilonWallFunctionProcess const& rOther);

    /// Copy constructor.
    RansEpsilonWallFunctionProcess(RansEpsilonWallFunctionProcess const& rOther);

    ///@}

}; // Class RansEpsilonWallFunctionProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEpsilonWallFunctionProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EPSILON_WALL_FUNCTION_PROCESS_H_INCLUDED defined
