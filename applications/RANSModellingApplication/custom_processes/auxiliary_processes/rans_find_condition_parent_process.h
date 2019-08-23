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

#if !defined(KRATOS_RANS_FIND_CONDITION_PARENT_PROCESS_H_INCLUDED)
#define KRATOS_RANS_FIND_CONDITION_PARENT_PROCESS_H_INCLUDED

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

/**
 * @brief Finds parent element of conditions
 *
 * This process finds parent element of conditions
 *
 */

class RansFindConditionParentProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansFindConditionParentProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansFindConditionParentProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansFindConditionParentProcess(Model& rModel, Parameters& rParameters)
        : mrModel(rModel), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"                : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"                     : 0
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mModelPartName = mrParameters["model_part_name"].GetString();
        mEchoLevel = mrParameters["echo_level"].GetInt();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansFindConditionParentProcess() override
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
        KRATOS_TRY

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
        for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
            SetConditionParent(*(r_model_part.ConditionsBegin() + i_condition));

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Found parents for conditions in " << mModelPartName << ".\n";

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
        return std::string("RansFindConditionParentProcess");
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void SetConditionParent(Condition& rCondition)
    {
        KRATOS_TRY

        GlobalPointersVector<Element> element_candidates;
        const Condition::GeometryType& r_condition_geometry = rCondition.GetGeometry();
        const int number_of_condition_nodes = r_condition_geometry.PointsNumber();

        std::vector<int> node_ids(number_of_condition_nodes), element_node_ids;

        for (int i_node = 0; i_node < number_of_condition_nodes; ++i_node)
        {
            const NodeType& r_node = r_condition_geometry[i_node];
            const GlobalPointersVector<Element>& r_node_element_candidates =
                r_node.GetValue(NEIGHBOUR_ELEMENTS);

            KRATOS_ERROR_IF(r_node_element_candidates.size() == 0)
                << "No neighbour elements found for node with id=" << r_node.Id()
                << " in condition with condition id=" << rCondition.Id()
                << " belongs to " << mModelPartName << ".\n";

            for (int j_element = 0;
                 j_element < static_cast<int>(r_node_element_candidates.size()); ++j_element)
            {
                element_candidates.push_back(r_node_element_candidates(j_element));
            }
            node_ids[i_node] = r_node.Id();
        }

        std::sort(node_ids.begin(), node_ids.end());

        for (int i_element = 0;
             i_element < static_cast<int>(element_candidates.size()); ++i_element)
        {
            const Element::GeometryType& r_geometry =
                element_candidates[i_element].GetGeometry();
            const int number_of_element_candidate_nodes = r_geometry.PointsNumber();
            if (static_cast<int>(element_node_ids.size()) != number_of_element_candidate_nodes)
                element_node_ids.resize(number_of_element_candidate_nodes);

            for (int i_node = 0; i_node < number_of_element_candidate_nodes; ++i_node)
            {
                element_node_ids[i_node] = r_geometry[i_node].Id();
            }

            std::sort(element_node_ids.begin(), element_node_ids.end());
            if (std::includes(element_node_ids.begin(), element_node_ids.end(),
                              node_ids.begin(), node_ids.end()))
            {
                rCondition.SetValue(PARENT_ELEMENT, element_candidates(i_element));
                return;
            }
        }

        KRATOS_ERROR << "Parent element for condition id=" << rCondition.Id()
                     << " not found in " << mModelPartName << ".\n";
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
    RansFindConditionParentProcess& operator=(RansFindConditionParentProcess const& rOther);

    /// Copy constructor.
    RansFindConditionParentProcess(RansFindConditionParentProcess const& rOther);

    ///@}

}; // Class RansFindConditionParentProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansFindConditionParentProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_FIND_CONDITION_PARENT_PROCESS_H_INCLUDED defined
