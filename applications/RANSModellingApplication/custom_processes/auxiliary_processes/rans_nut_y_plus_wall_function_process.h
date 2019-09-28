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

#if !defined(KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_PROCESS_H_INCLUDED

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

class RansNutYPlusWallFunctionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansNutYPlusWallFunctionProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutYPlusWallFunctionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutYPlusWallFunctionProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09,
            "von_karman"      : 0.41,
            "beta"            : 5.2
        })");

        mrParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = mrParameters["echo_level"].GetInt();
        mModelPartName = mrParameters["model_part_name"].GetString();
        mCmu = mrParameters["c_mu"].GetDouble();
        mVonKarman = mrParameters["von_karman"].GetDouble();
        mBeta = mrParameters["beta"].GetDouble();
        mLimitYPlus = RansCalculationUtilities().CalculateLogarithmicYPlusLimit(
            mVonKarman, mBeta);

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~RansNutYPlusWallFunctionProcess() override
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
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(TURBULENT_VISCOSITY);

        const ModelPart::NodesContainerType& r_nodes =
            mrModel.GetModelPart(mModelPartName).Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        }

        return 0;

        KRATOS_CATCH("");
    }

    void Execute() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const int number_of_nodes = r_model_part.NumberOfNodes();

        const double nu_t_min = r_model_part.GetProcessInfo()[TURBULENT_VISCOSITY_MIN];

        unsigned int number_of_modified_nu_t_nodes = 0;

#pragma omp parallel for reduction(+ : number_of_modified_nu_t_nodes)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

            if (y_plus > mLimitYPlus)
            {
                r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) =
                    mVonKarman * y_plus * nu;
                number_of_modified_nu_t_nodes++;
            }
            else
            {
                r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t_min;
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Applied nu_t y_plus wall function to " << number_of_modified_nu_t_nodes
            << " of total " << r_model_part.NumberOfNodes() << " nodes in "
            << mModelPartName << "\n";

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
        return std::string("RansNutYPlusWallFunctionProcess");
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

    double mLimitYPlus;

    double mCmu;
    double mVonKarman;
    double mBeta;

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
    RansNutYPlusWallFunctionProcess& operator=(RansNutYPlusWallFunctionProcess const& rOther);

    /// Copy constructor.
    RansNutYPlusWallFunctionProcess(RansNutYPlusWallFunctionProcess const& rOther);

    ///@}

}; // Class RansNutYPlusWallFunctionProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansNutYPlusWallFunctionProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_PROCESS_H_INCLUDED defined
