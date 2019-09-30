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

#if !defined(KRATOS_RANS_APPLY_EXACT_NODAL_PERIODIC_CONDITION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_APPLY_EXACT_NODAL_PERIODIC_CONDITION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "custom_utilities/periodic_condition_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "processes/reorder_and_optimize_modelpart_process.h"
#include "utilities/openmp_utils.h"

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
 * @brief A process to create periodic conditions between two boundaries
 *
 * This process does node to node matching to create periodic boundary
 * conditions between the master and slave model parts, to be used in ResidualBasedBlockBuilderAndSolverPeriodic.
 *
 * This process requires two boundaries to have maching meshes.
 * Make sure to avoid calling "ReplaceElementsAndConditions" method call, and directly impose element
 * and condition name in the "*.mdpa" file in order to correctly use this process alongside with
 * ResidualBasedBlockBuilderAndSolverPeriodic
 *
 * @see class ResidualBasedBlockBuilderAndSolverPeriodic
 */
class RansApplyExactNodalPeriodicConditionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansApplyExactNodalPeriodicConditionProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansApplyExactNodalPeriodicConditionProcess);

    using NodeType = ModelPart::NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansApplyExactNodalPeriodicConditionProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "base_model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "master_model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "slave_model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_names_list"            : [],
            "tolerance"                      : 1e-9,
            "translation_settings" :
            {
                "direction" : [0.0, 0.0, 0.0],
                "magnitude" : 0.0
            },
            "rotation_settings" :
            {
                "axis"   : [0.0, 0.0, 0.0],
                "center" : [0.0, 0.0, 0.0],
                "angle"  : 0.0
            },
            "echo_level"                     : 0,
            "reorder"                        : true
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mBaseModelPartName = mrParameters["base_model_part_name"].GetString();
        mMasterModelPartName = mrParameters["master_model_part_name"].GetString();
        mSlaveModelPartName = mrParameters["slave_model_part_name"].GetString();
        mVariablesList = mrParameters["variable_names_list"].GetStringArray();
        mReorder = mrParameters["reorder"].GetBool();
        mTolerance = mrParameters["tolerance"].GetDouble();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        const double eps = std::numeric_limits<double>::epsilon();

        KRATOS_ERROR_IF(mVariablesList.size() == 0)
            << "No variables are provided. Please specify at least one "
               "variable in \"variable_names_list\" to use periodic "
               "conditions.\n";

        // translation settings
        noalias(mTranslationDirection) =
            mrParameters["translation_settings"]["direction"].GetVector();
        const double translation_direction_norm = norm_2(mTranslationDirection);
        if (translation_direction_norm > eps)
            noalias(mTranslationDirection) = mTranslationDirection / translation_direction_norm;
        mTranslationMagnitude =
            mrParameters["translation_settings"]["magnitude"].GetDouble();

        // rotation settings
        noalias(mRotationCenter) =
            mrParameters["rotation_settings"]["center"].GetVector();
        noalias(mRotationAxis) = mrParameters["rotation_settings"]["axis"].GetVector();
        const double rotation_axis_norm = norm_2(mRotationAxis);
        if (rotation_axis_norm > eps)
            noalias(mRotationAxis) = mRotationAxis / rotation_axis_norm;
        mRotationAngle = mrParameters["rotation_settings"]["angle"].GetDouble();

        KRATOS_ERROR_IF(mRotationAngle < eps && mTranslationMagnitude < eps)
            << "Either translation or rotation should be defined.\n";
        KRATOS_ERROR_IF(mRotationAngle > eps && rotation_axis_norm < eps)
            << "Rotation angle defined, but rotation axis is not defined.\n";
        KRATOS_ERROR_IF(mTranslationMagnitude > eps && translation_direction_norm < eps) << "Translation magnitude defined, but translation direction is not defined.\n";

        const int domain_size =
            mrModel.GetModelPart(mBaseModelPartName).GetProcessInfo()[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size == 2 && (mRotationAxis[0] > eps || mRotationAxis[1] > eps))
            << "2D rotation axis should be [0.0, 0.0, 1.0].\n";

        if (mRotationAngle > eps && mTranslationMagnitude > eps)
            noalias(mRotationCenter) += mTranslationDirection * mTranslationMagnitude;

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansApplyExactNodalPeriodicConditionProcess() override
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

        RansCheckUtilities rans_check_utilities;

        rans_check_utilities.CheckIfModelPartExists(mrModel, mBaseModelPartName);
        rans_check_utilities.CheckIfModelPartExists(mrModel, mMasterModelPartName);
        rans_check_utilities.CheckIfModelPartExists(mrModel, mSlaveModelPartName);

        ModelPart::NodesContainerType& r_base_model_part_nodes =
            mrModel.GetModelPart(mBaseModelPartName).Nodes();

        for (std::string variable_name : mVariablesList)
        {
            if (KratosComponents<Variable<double>>::Has(variable_name))
            {
                const Variable<double>& variable =
                    KratosComponents<Variable<double>>::Get(variable_name);
                rans_check_utilities.CheckIfVariableExistsInNodesContainer(
                    r_base_model_part_nodes, variable);
            }
            else if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(
                         variable_name))
            {
                const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& variable =
                    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                        variable_name);
                rans_check_utilities.CheckIfVariableExistsInNodesContainer(
                    r_base_model_part_nodes, variable);
            }
            else
            {
                KRATOS_ERROR << "Variable " << variable_name << " not found.\n";
            }
        }

        return 0;

        KRATOS_CATCH("");
    }

    void ExecuteInitialize() override
    {
        CreatePeriodicConditions();
        if (mReorder)
        {
            ModelPart& r_model_part = mrModel.GetModelPart(mBaseModelPartName);
            Parameters default_params(R"({})");
            ReorderAndOptimizeModelPartProcess reoder_process(r_model_part, default_params);
            reoder_process.Execute();
        }
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
        return std::string("RansApplyExactNodalPeriodicConditionProcess");
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

    std::string mBaseModelPartName;
    std::string mMasterModelPartName;
    std::string mSlaveModelPartName;

    std::vector<std::string> mVariablesList;
    double mTolerance;

    // translation settings
    array_1d<double, 3> mTranslationDirection;
    double mTranslationMagnitude;

    // rotation settings
    array_1d<double, 3> mRotationAxis;
    array_1d<double, 3> mRotationCenter;
    double mRotationAngle;

    bool mReorder;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CreatePeriodicConditions()
    {
        KRATOS_TRY

        ModelPart& r_base_model_part = mrModel.GetModelPart(mBaseModelPartName);
        int condition_id = r_base_model_part.NumberOfConditions();
        Properties::Pointer p_properties = r_base_model_part.CreateNewProperties(
            r_base_model_part.NumberOfProperties() + 1);

        for (std::string variable_name : mVariablesList)
        {
            if (KratosComponents<Variable<double>>::Has(variable_name))
            {
                const Variable<double>& variable =
                    KratosComponents<Variable<double>>::Get(variable_name);
                p_properties->GetValue(PERIODIC_VARIABLES).Add(variable);
            }
            else if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(
                         variable_name))
            {
                const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& variable =
                    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                        variable_name);
                p_properties->GetValue(PERIODIC_VARIABLES).Add(variable);
            }
            else
            {
                KRATOS_ERROR << "Variable " << variable_name << " not found.\n";
            }
        }

        ModelPart::NodesContainerType& r_master_model_part_nodes =
            mrModel.GetModelPart(mMasterModelPartName).Nodes();
        ModelPart::NodesContainerType& r_slave_model_part_nodes =
            mrModel.GetModelPart(mSlaveModelPartName).Nodes();

        KRATOS_ERROR_IF(r_master_model_part_nodes.size() !=
                        r_slave_model_part_nodes.size())
            << this->Info() << " failed. Master and slave model part nodes mismatch.\n";

        const int number_of_nodes = r_master_model_part_nodes.size();
        const double eps = std::numeric_limits<double>::epsilon();
        const array_1d<double, 3> translation = mTranslationDirection * mTranslationMagnitude;

        if (mRotationAngle > eps && mTranslationMagnitude > eps)
        {
#pragma omp parallel for shared(condition_id)
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_master_node = *(r_master_model_part_nodes.begin() + i_node);
                const array_1d<double, 3>& master_initial_position =
                    r_master_node.GetInitialPosition().Coordinates();
                const array_1d<double, 3>& master_final_position =
                    CalculateRotatedPosition(master_initial_position + translation);
                int& r_master_patch_index =
                    r_master_node.FastGetSolutionStepValue(PATCH_INDEX);
                for (int j_node = 0; j_node < number_of_nodes; ++j_node)
                {
                    const NodeType& r_slave_node =
                        *(r_slave_model_part_nodes.begin() + j_node);
                    if (norm_2(master_final_position -
                               r_slave_node.GetInitialPosition().Coordinates()) < mTolerance)
                    {
                        r_master_patch_index = r_slave_node.Id();
                        break;
                    }
                }

                KRATOS_ERROR_IF(r_master_patch_index == 0)
                    << "Slave node is not found in " << mSlaveModelPartName
                    << " near " << master_final_position << " with tolerance "
                    << mTolerance << " for master node id=" << r_master_node.Id()
                    << " at " << r_master_node.GetInitialPosition().Coordinates()
                    << " in " << mMasterModelPartName << ".\n";

                if (static_cast<std::size_t>(r_master_patch_index) !=
                    r_master_node.Id())
                {
                    std::vector<std::size_t> node_id_list = {
                        r_master_node.Id(), static_cast<std::size_t>(r_master_patch_index)};
#pragma omp critical
                    {
                        condition_id++;
                        Condition::Pointer p_condition = r_base_model_part.CreateNewCondition(
                            "PeriodicCondition", condition_id, node_id_list, p_properties);
                        p_condition->Set(PERIODIC);
                        r_master_node.Set(PERIODIC);
                    }
                }
            }
        }
        else if (mRotationAngle > eps)
        {
#pragma omp parallel for shared(condition_id)
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_master_node = *(r_master_model_part_nodes.begin() + i_node);
                const array_1d<double, 3>& master_initial_position =
                    r_master_node.GetInitialPosition().Coordinates();
                const array_1d<double, 3>& master_final_position =
                    CalculateRotatedPosition(master_initial_position);
                int& r_master_patch_index =
                    r_master_node.FastGetSolutionStepValue(PATCH_INDEX);
                for (int j_node = 0; j_node < number_of_nodes; ++j_node)
                {
                    const NodeType& r_slave_node =
                        *(r_slave_model_part_nodes.begin() + j_node);
                    if (norm_2(master_final_position -
                               r_slave_node.GetInitialPosition().Coordinates()) < mTolerance)
                    {
                        r_master_patch_index = r_slave_node.Id();
                        break;
                    }
                }

                KRATOS_ERROR_IF(r_master_patch_index == 0)
                    << "Slave node is not found in " << mSlaveModelPartName
                    << " near " << master_final_position << " with tolerance "
                    << mTolerance << " for master node id=" << r_master_node.Id()
                    << " at " << r_master_node.GetInitialPosition().Coordinates()
                    << " in " << mMasterModelPartName << ".\n";

                if (static_cast<std::size_t>(r_master_patch_index) !=
                    r_master_node.Id())
                {
                    std::vector<std::size_t> node_id_list = {
                        r_master_node.Id(), static_cast<std::size_t>(r_master_patch_index)};
#pragma omp critical
                    {
                        condition_id++;
                        Condition::Pointer p_condition = r_base_model_part.CreateNewCondition(
                            "PeriodicCondition", condition_id, node_id_list, p_properties);
                        p_condition->Set(PERIODIC);
                        r_master_node.Set(PERIODIC);
                    }
                }
            }
        }
        else if (mTranslationMagnitude > eps)
        {
#pragma omp parallel for shared(condition_id)
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                NodeType& r_master_node = *(r_master_model_part_nodes.begin() + i_node);
                const array_1d<double, 3>& master_initial_position =
                    r_master_node.GetInitialPosition().Coordinates();
                const array_1d<double, 3>& master_final_position =
                    master_initial_position + translation;
                int& r_master_patch_index =
                    r_master_node.FastGetSolutionStepValue(PATCH_INDEX);
                for (int j_node = 0; j_node < number_of_nodes; ++j_node)
                {
                    const NodeType& r_slave_node =
                        *(r_slave_model_part_nodes.begin() + j_node);
                    if (norm_2(master_final_position -
                               r_slave_node.GetInitialPosition().Coordinates()) < mTolerance)
                    {
                        r_master_patch_index = static_cast<int>(r_slave_node.Id());
                        break;
                    }
                }

                KRATOS_ERROR_IF(r_master_patch_index == 0)
                    << "Slave node is not found in " << mSlaveModelPartName
                    << " near " << master_final_position << " with tolerance "
                    << mTolerance << " for master node id=" << r_master_node.Id()
                    << " at " << r_master_node.GetInitialPosition().Coordinates()
                    << " in " << mMasterModelPartName << ".\n";

                if (static_cast<std::size_t>(r_master_patch_index) !=
                    r_master_node.Id())
                {
                    std::vector<std::size_t> node_id_list = {
                        r_master_node.Id(), static_cast<std::size_t>(r_master_patch_index)};
#pragma omp critical
                    {
                        condition_id++;
                        Condition::Pointer p_condition = r_base_model_part.CreateNewCondition(
                            "PeriodicCondition", condition_id, node_id_list, p_properties);
                        p_condition->Set(PERIODIC);
                        r_master_node.Set(PERIODIC);
                    }
                }
            }
        }

        const int number_of_conditions = r_base_model_part.NumberOfConditions();
#pragma omp parallel for
        for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
        {
            Condition& r_condition = *(r_base_model_part.ConditionsBegin() + i_condition);
            if (r_condition.Is(PERIODIC))
            {
                const NodeType& r_node_master = r_condition.GetGeometry()[0];
                NodeType& r_node_slave = r_condition.GetGeometry()[1];

                r_node_slave.SetLock();
                r_node_slave.FastGetSolutionStepValue(PATCH_INDEX) =
                    r_node_master.Id();
                r_node_slave.Set(PERIODIC);
                r_node_slave.UnSetLock();
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Created periodic conditions between " << mMasterModelPartName << " and "
            << mSlaveModelPartName << " in " << mBaseModelPartName << ".\n";

        KRATOS_CATCH("");
    }

    array_1d<double, 3> CalculateRotatedPosition(const array_1d<double, 3>& rInitialPosition) const
    {
        BoundedMatrix<double, 3, 3> rotation_matrix;
        const array_1d<double, 3> local_coord = rInitialPosition - mRotationCenter;
        CalculateRotationMatrix(rotation_matrix);
        const array_1d<double, 3> global_coord =
            mRotationCenter + prod(rotation_matrix, local_coord);

        return global_coord;
    }

    void CalculateRotationMatrix(BoundedMatrix<double, 3, 3>& rOutput) const
    {
        const double ux = mRotationAxis[0];
        const double uy = mRotationAxis[1];
        const double uz = mRotationAxis[2];
        const double c2 = std::cos(mRotationAngle);
        const double s2 = std::sin(mRotationAngle);
        const double s1 = std::sin(mRotationAngle / 2.);

        rOutput(0, 0) = c2 + 2. * ux * ux * s1 * s1;
        rOutput(1, 1) = c2 + 2. * uy * uy * s1 * s1;
        rOutput(2, 2) = c2 + 2. * uz * uz * s1 * s1;
        rOutput(0, 1) = 2. * ux * uy * s1 * s1 - uz * s2;
        rOutput(1, 0) = 2. * ux * uy * s1 * s1 + uz * s2;
        rOutput(0, 2) = 2. * ux * uz * s1 * s1 + uy * s2;
        rOutput(2, 0) = 2. * ux * uz * s1 * s1 - uy * s2;
        rOutput(1, 2) = 2. * uy * uz * s1 * s1 - ux * s2;
        rOutput(2, 1) = 2. * uy * uz * s1 * s1 + ux * s2;
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
    RansApplyExactNodalPeriodicConditionProcess& operator=(
        RansApplyExactNodalPeriodicConditionProcess const& rOther);

    /// Copy constructor.
    RansApplyExactNodalPeriodicConditionProcess(RansApplyExactNodalPeriodicConditionProcess const& rOther);

    ///@}

}; // Class RansApplyExactNodalPeriodicConditionProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansApplyExactNodalPeriodicConditionProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_APPLY_EXACT_NODAL_PERIODIC_CONDITION_PROCESS_H_INCLUDED defined
