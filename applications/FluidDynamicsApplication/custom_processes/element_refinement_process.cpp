//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <string>
#include <vector>
#include <functional>
#include <algorithm>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "processes/process.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"

// Application incldues
#include "meshing_application_variables.h"

// Include base h
#include "element_refinement_process.h"

namespace Kratos
{

namespace ElementRefinementProcessHelperUtilities
{

using IndexType = std::size_t;

using ElementType = ModelPart::ElementType;

using ConditionType = ModelPart::ConditionType;

template<class TEntityType>
void CopyEntityData(
    TEntityType& rOutput,
    const TEntityType& rInput)
{
    KRATOS_TRY

    rOutput.SetFlags(rInput.GetFlags());
    rOutput.SetData(rInput.GetData());
    rOutput.SetProperties(rInput.pGetProperties());
    rOutput.Set(ACTIVE, true);

    KRATOS_CATCH("");
}

void FindElementSurfaceConditions(ModelPart& rModelPart)
{
    KRATOS_TRY

    block_for_each(rModelPart.Elements(), [](ElementType& rElement) {
        rElement.SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<ConditionType>());
    });

    block_for_each(rModelPart.Conditions(), [](ConditionType& rCondition) {
        if (rCondition.Has(NEIGHBOUR_ELEMENTS)) {
            auto& parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];

#pragma omp critical
            {
                parent_element.GetValue(NEIGHBOUR_CONDITIONS).push_back(GlobalPointer<ConditionType>(&rCondition));
            }
        }
    });

    KRATOS_CATCH("");
}

void AddVariableToList(
    std::vector<const Variable<double>*>& rVariablesList,
    const std::string& rVariableName)
{
    KRATOS_TRY

    if (KratosComponents<Variable<double>>::Has(rVariableName)) {
        rVariablesList.push_back(&KratosComponents<Variable<double>>::Get(rVariableName));
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)) {
        rVariablesList.push_back(&KratosComponents<Variable<double>>::Get(rVariableName + "_X"));
        rVariablesList.push_back(&KratosComponents<Variable<double>>::Get(rVariableName + "_Y"));
        rVariablesList.push_back(&KratosComponents<Variable<double>>::Get(rVariableName + "_Z"));
    } else {
        KRATOS_ERROR << "Unsupported variable type provided. Only scalar and "
                        "3D vectors are supported. [ variable_name = "
                     << rVariableName << " ].\n";
    }

    KRATOS_CATCH("");
}

// template instantiations
template void CopyEntityData<ConditionType>(ConditionType&, const ConditionType&);
template void CopyEntityData<ElementType>(ElementType&, const ElementType&);

}

const Parameters ElementRefinementProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"                     : "PLEASE_SPECIFY_INPUT_MODEL_PART_NAME",
            "refinement_model_part_prefix"        : "Thread_",
            "refinement_level"                    : 1,
            "echo_level"                          : 0,
            "refined_element_name"                : "PLEASE_SPECIFY_REFINED_ELEMENT_NAME",
            "refined_condition_name"              : "PLEASE_SPECIFY_REFINED_CONDITION_NAME",
            "nodal_interpolation_settings": {
                "historical_variables_list"    : ["ALL_VARIABLES_FROM_VARIABLES_LIST"],
                "non_hitsorical_variables_list": [],
                "flags_list"                   : [],
                "historical_buffer_size"       : 1
            }
        })");

    return default_parameters;
}

ElementRefinementProcess::ElementRefinementProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.RecursivelyValidateAndAssignDefaults(GetDefaultParameters());

    mModelPartName = rParameters["model_part_name"].GetString();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mRefinedModelPartNamePrefix = mModelPartName + "_" + rParameters["refinement_model_part_prefix"].GetString();
    std::replace(mRefinedModelPartNamePrefix.begin(), mRefinedModelPartNamePrefix.end(), '.', ':');

    mRefinementLevel = rParameters["refinement_level"].GetInt();

    KRATOS_ERROR_IF(mRefinementLevel <= 0)
        << "Refinement level should be greater than zero. [ refinement_level = " << mRefinementLevel
        << " ].\n";

    mRefinedElementName = rParameters["refined_element_name"].GetString();
    mRefinedConditionName = rParameters["refined_condition_name"].GetString();

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const auto& nodal_interpolation_settings = rParameters["nodal_interpolation_settings"];
    auto historical_variables_names_list = nodal_interpolation_settings["historical_variables_list"].GetStringArray();
    mHistoricalVariableInterpolationBufferSize = nodal_interpolation_settings["historical_buffer_size"].GetInt();

    if (historical_variables_names_list.size() == 1 && historical_variables_names_list[0] == "ALL_VARIABLES_FROM_VARIABLES_LIST") {
        historical_variables_names_list.clear();
        for (const auto& r_variable : r_model_part.GetNodalSolutionStepVariablesList()) {
            historical_variables_names_list.push_back(r_variable.Name());
        }
    }

    for (const auto& r_var_name : historical_variables_names_list) {
        ElementRefinementProcessHelperUtilities::AddVariableToList(mNodalHistoricalVariablesList, r_var_name);
    }

    for (const auto& r_var_name : nodal_interpolation_settings["non_hitsorical_variables_list"].GetStringArray()) {
        ElementRefinementProcessHelperUtilities::AddVariableToList(mNodalNonHistoricalVariablesList, r_var_name);
    }

    for (const auto& r_flag_name : nodal_interpolation_settings["flags_list"].GetStringArray()) {
        if (KratosComponents<Flags>::Has(r_flag_name)) {
            mNodalFlagsList.push_back(&KratosComponents<Flags>::Get(r_flag_name));
        } else {
            KRATOS_ERROR << "Flag name not defined. [ flag name = " << r_flag_name << " ].\n";
        }
    }

    for (const auto p_var : mNodalHistoricalVariablesList) {
        KRATOS_ERROR_IF_NOT(r_model_part.HasNodalSolutionStepVariable(*p_var))
            << p_var->Name() << " not found in nodal historical variables list in "
            << r_model_part.Name() << ".\n";
    }

    std::stringstream msg;
    msg << "Interpolation Settings:";
    msg << "\n        Historical variable interpolation buffer size: " << mHistoricalVariableInterpolationBufferSize;
    msg << "\n        Historical variables:";
    for (const auto p_var : mNodalHistoricalVariablesList){
        msg << "\n            " << p_var->Name();
    }
    msg << "\n        NonHistorical variables:";
    for (const auto p_var : mNodalNonHistoricalVariablesList){
        msg << "\n            " << p_var->Name();
    }
    msg << "\n        Flags:";
    for (const auto& r_flag_name : nodal_interpolation_settings["flags_list"].GetStringArray()){
        msg << "\n            " << r_flag_name;
    }

    KRATOS_INFO(this->Info()) << msg.str() << std::endl;

    KRATOS_CATCH("");
}

ElementRefinementProcess::ThreadLocalStorage& ElementRefinementProcess::GetThreadLocalStorage()
{
    return mThreadLocalStorage[OpenMPUtils::ThisThread()];
}

const std::string ElementRefinementProcess::GetThreadLocalModelPartName(const IndexType ThreadId) const
{
    return mRefinedModelPartNamePrefix + std::to_string(ThreadId + 1);
}

void ElementRefinementProcess::CreateSurfaceCondition(
    ThreadLocalStorage& rTLS,
    const IndexType SurfaceIndex,
    Properties::Pointer pProperties,
    Element::Pointer pParentElement,
    const IndexType ConditionId,
    const std::vector<IndexType>& rNodeIds)
{
    KRATOS_TRY

    auto p_surface_condition = rTLS.pRefinedModelPart->CreateNewCondition(mRefinedConditionName, ConditionId, rNodeIds, pProperties);
    rTLS.mRefinedModelPartSurfaceParentElements[SurfaceIndex].push_back(std::make_pair(&(*p_surface_condition), pParentElement));

    auto p_current_surface_model_part = rTLS.mSurfaceModelPartList[SurfaceIndex];
    p_current_surface_model_part->AddCondition(p_surface_condition);
    p_current_surface_model_part->AddNodes(rNodeIds);

    KRATOS_CATCH("");
}

void ElementRefinementProcess::CreateSurfaceAuxiliaries(ThreadLocalStorage& rTLS)
{
    KRATOS_TRY

    const IndexType number_of_surfaces = mSurfacesNodeIndexOrdering.size();

    rTLS.mSurfaceModelPartList.resize(number_of_surfaces);
    rTLS.mRefinedModelPartSurfaceParentElements.resize(number_of_surfaces);

    for (IndexType i = 0; i < number_of_surfaces; ++i) {
        rTLS.mSurfaceModelPartList[i] = &rTLS.pRefinedModelPart->CreateSubModelPart("Surface_" + std::to_string(i + 1));
    }

    KRATOS_CATCH("");
}

template <>
void ElementRefinementProcess::RefineElement<2, 3, 3>(
    ThreadLocalStorage& rTLS,
    const IndexType RefinementLevel,
    IndexType& rCurrentNodeId,
    IndexType& rCurrentConditionId,
    IndexType& rCurrentElementId,
    const BMatrixNM<3, 3>& rNodalInterpolationValues,
    const SurfaceIndicesArray<3>& rSurfaceIndices)
{
    KRATOS_TRY

    if (RefinementLevel == 0) {
        // add all refined elements, conditions and nodes
        const auto& r_geometry = mReferenceCoarseElement->GetGeometry();

        // interpolate nodal locations
        BMatrixNM<3, 3> reference_coordinates, interpolated_coordinates;
        row(reference_coordinates, 0) = r_geometry[0].Coordinates();
        row(reference_coordinates, 1) = r_geometry[1].Coordinates();
        row(reference_coordinates, 2) = r_geometry[2].Coordinates();

        noalias(interpolated_coordinates) = prod(rNodalInterpolationValues, reference_coordinates);

        // create nodes
        std::vector<IndexType> node_ids(3, 0);
        for (IndexType i = 0; i < 3; ++i) {
            auto& r_current_id = node_ids[i];
            const Array3D& current_coordinates = row(interpolated_coordinates, i);

            for (auto& r_node : rTLS.pRefinedModelPart->Nodes()) {
                if (norm_2(r_node.Coordinates() - current_coordinates) <= std::numeric_limits<double>::epsilon()) {
                    r_current_id = r_node.Id();
                    break;
                }
            }

            if (r_current_id == 0) {
                r_current_id = rCurrentNodeId++;
                auto p_node = rTLS.pRefinedModelPart->CreateNewNode(r_current_id, current_coordinates[0], current_coordinates[1], current_coordinates[2]);
                p_node->SetValue(NODAL_INTERPOLATION_VALUES, row(rNodalInterpolationValues, i));
            }
        }

        auto p_properties = mReferenceCoarseElement->pGetProperties();

        // create elements
        auto p_element = &*(rTLS.pRefinedModelPart->CreateNewElement(mRefinedElementName, rCurrentElementId++, node_ids, p_properties));

        // create surface conditions
        for (IndexType i = 0; i < 3; ++i) {
            // surface node indices should make the normal outward pointing
            const IndexType current_surface_index = rSurfaceIndices[i];

            if (current_surface_index < 3) {
                CreateSurfaceCondition(rTLS, current_surface_index, p_properties,
                                       p_element, rCurrentConditionId++,
                                       {node_ids[i % 3], node_ids[(i + 1) % 3]});
            }
        }
    } else {
        /*       a(0)
                  *
                 / \
                / 1 \
         S0   2 *-----* 1  S2
              / \ 4 / \
             / 2 \ / 3 \
             *----*-----*
           b(1)   0     c(2)
                  S1
        */
        // rows index local node, cols coarse mesh node index (a, b, c)
        BMatrixNM<3, 3> sub_triangle_interpolation_values;

        // triangle 1
        row(sub_triangle_interpolation_values, 0) = row(rNodalInterpolationValues, 0);
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 2)) * 0.5;
        this->RefineElement<2, 3, 3>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<3>({rSurfaceIndices[0], 4, rSurfaceIndices[2]}));

        // triangle 2
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = row(rNodalInterpolationValues, 1);
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        this->RefineElement<2, 3, 3>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<3>({rSurfaceIndices[0], rSurfaceIndices[1], 4}));

        // triangle 3
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = row(rNodalInterpolationValues, 2);
        this->RefineElement<2, 3, 3>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<3>({4, rSurfaceIndices[1], rSurfaceIndices[2]}));

        // triangle 4
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 0)) * 0.5;
        this->RefineElement<2, 3, 3>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<3>({4, 4, 4}));
    }

    KRATOS_CATCH("");
}

template <>
void ElementRefinementProcess::RefineElement<3, 4, 4>(
    ThreadLocalStorage& rTLS,
    const IndexType RefinementLevel,
    IndexType& rCurrentNodeId,
    IndexType& rCurrentConditionId,
    IndexType& rCurrentElementId,
    const BMatrixNM<4, 4>& rNodalInterpolationValues,
    const SurfaceIndicesArray<4>& rSurfaceIndices)
{
    KRATOS_TRY

    if (RefinementLevel == 0) {
        // add all refined elements, conditions and nodes
        const auto& r_geometry = mReferenceCoarseElement->GetGeometry();

        // interpolate nodal locations
        BMatrixNM<4, 3> reference_coordinates, interpolated_coordinates;
        row(reference_coordinates, 0) = r_geometry[0].Coordinates();
        row(reference_coordinates, 1) = r_geometry[1].Coordinates();
        row(reference_coordinates, 2) = r_geometry[2].Coordinates();
        row(reference_coordinates, 3) = r_geometry[3].Coordinates();

        noalias(interpolated_coordinates) = prod(rNodalInterpolationValues, reference_coordinates);

        // create nodes
        std::vector<IndexType> node_ids(4, 0);
        for (IndexType i = 0; i < 4; ++i) {
            auto& r_current_id = node_ids[i];
            const Array3D& current_coordinates = row(interpolated_coordinates, i);

            for (auto& r_node : rTLS.pRefinedModelPart->Nodes()) {
                if (norm_2(r_node.Coordinates() - current_coordinates) <= std::numeric_limits<double>::epsilon()) {
                    r_current_id = r_node.Id();
                    break;
                }
            }

            if (r_current_id == 0) {
                r_current_id = rCurrentNodeId++;
                auto p_node = rTLS.pRefinedModelPart->CreateNewNode(r_current_id, current_coordinates[0], current_coordinates[1], current_coordinates[2]);
                p_node->SetValue(NODAL_INTERPOLATION_VALUES, row(rNodalInterpolationValues, i));
            }
        }

        auto p_properties = mReferenceCoarseElement->pGetProperties();

        // create elements
        auto p_element = &*(rTLS.pRefinedModelPart->CreateNewElement(mRefinedElementName, rCurrentElementId++, node_ids, p_properties));

        // create surface conditions
        // create surface 1
        IndexType surface_index = rSurfaceIndices[0];
        if (surface_index < 4) {
            // surface node indices should make the normal outward pointing
            CreateSurfaceCondition(rTLS, surface_index, p_properties, p_element,
                                   rCurrentConditionId++,
                                   {node_ids[0], node_ids[1], node_ids[3]});
        }

        // create surface 2
        surface_index = rSurfaceIndices[1];
        if (surface_index < 4) {
            // surface node indices should make the normal outward pointing
            CreateSurfaceCondition(rTLS, surface_index, p_properties, p_element,
                                   rCurrentConditionId++,
                                   {node_ids[3], node_ids[1], node_ids[2]});
        }

        // create surface 3
        surface_index = rSurfaceIndices[2];
        if (surface_index < 4) {
            // surface node indices should make the normal outward pointing
            CreateSurfaceCondition(rTLS, surface_index, p_properties, p_element,
                                   rCurrentConditionId++,
                                   {node_ids[0], node_ids[3], node_ids[2]});
        }

        // create surface 4
        surface_index = rSurfaceIndices[3];
        if (surface_index < 4) {
            // surface node indices should make the normal outward pointing
            CreateSurfaceCondition(rTLS, surface_index, p_properties, p_element,
                                   rCurrentConditionId++,
                                   {node_ids[1], node_ids[0], node_ids[2]});
        }
    } else {
        /*                   0
                             *
                           / | \
                          /  |  \    (d - surface)
                         /   |   \
                        /  a |  c \
                       /     *     \
                      /      3      \
                     /               \
                    /        b        \
                   * _________________ *
                   1                    2
        */
        // above, a (side),b(bottom),c (side),d (outer surface) are surfaces and 0,1,2,3 are node index order
        // rows index local node, cols coarse mesh node index (a, b, c)
        BMatrixNM<4, 4> sub_triangle_interpolation_values;

        // refine for four corner tetrahedra
        // corner tetrahedra 1 (top corner)
        row(sub_triangle_interpolation_values, 0) = row(rNodalInterpolationValues, 0);
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 3) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 3)) * 0.5;
        this->RefineElement<3, 4, 4>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<4>({rSurfaceIndices[0], 5, rSurfaceIndices[2], rSurfaceIndices[3]}));

        // corner tetrahedra 2 (bottom front left)
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 0)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = row(rNodalInterpolationValues, 1);
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 3) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 3)) * 0.5;
        this->RefineElement<3, 4, 4>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<4>({rSurfaceIndices[0], rSurfaceIndices[1], 5, rSurfaceIndices[3]}));

        // corner tetrahedra 3 (bottom front right)
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 2) + row(rNodalInterpolationValues, 0)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 2) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = row(rNodalInterpolationValues, 2);
        row(sub_triangle_interpolation_values, 3) = (row(rNodalInterpolationValues, 2) + row(rNodalInterpolationValues, 3)) * 0.5;
        this->RefineElement<3, 4, 4>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<4>({5, rSurfaceIndices[1], rSurfaceIndices[2], rSurfaceIndices[3]}));

        // corner tetrahedra 4 (bottom back)
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 3) + row(rNodalInterpolationValues, 0)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 3) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 3) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 3) = row(rNodalInterpolationValues, 3);
        this->RefineElement<3, 4, 4>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<4>({rSurfaceIndices[0], rSurfaceIndices[1], rSurfaceIndices[2], 5}));

        // create sub triangles by dividing octagon
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 3)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 3) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 1)) * 0.5;
        this->RefineElement<3, 4, 4>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<4>({5, 5, rSurfaceIndices[3], 5}));

        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 3)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 3)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 0)) * 0.5;
        row(sub_triangle_interpolation_values, 3) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        this->RefineElement<3, 4, 4>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<4>({5, 5, 5, rSurfaceIndices[0]}));

        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 3) + row(rNodalInterpolationValues, 0)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 2) + row(rNodalInterpolationValues, 3)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 2) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 3) = (row(rNodalInterpolationValues, 2) + row(rNodalInterpolationValues, 0)) * 0.5;
        this->RefineElement<3, 4, 4>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<4>({rSurfaceIndices[2], 5, 5, 5}));

        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 2) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 3) + row(rNodalInterpolationValues, 0)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 3) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 3) = (row(rNodalInterpolationValues, 3) + row(rNodalInterpolationValues, 2)) * 0.5;
        this->RefineElement<3, 4, 4>(rTLS, RefinementLevel - 1, rCurrentNodeId, rCurrentConditionId, rCurrentElementId, sub_triangle_interpolation_values, SurfaceIndicesArray<4>({5, 5, rSurfaceIndices[1], 5}));
    }

    KRATOS_CATCH("");
}

void ElementRefinementProcess::ComputeSurfaceMap()
{
    KRATOS_TRY

    using namespace ElementRefinementProcessHelperUtilities;

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    FindElementSurfaceConditions(r_model_part);

    // first create all the maps for each element
    for (const auto& r_element : r_model_part.Elements()) {
        mElementSurfaceIndicesMap[r_element.Id()];
    }

    const IndexType number_of_surfaces = mSurfacesNodeIndexOrdering.size();

    block_for_each(r_model_part.Elements(), [&](Element& rCoarseElement) {
        auto& r_element_surface_map = mElementSurfaceIndicesMap.find(rCoarseElement.Id())->second;
        for (auto p_surface_condition : rCoarseElement.GetValue(NEIGHBOUR_CONDITIONS).GetContainer()) {
            const IndexType surface_id = p_surface_condition->Id();
            const auto& r_element_geometry = rCoarseElement.GetGeometry();
            const auto& r_surface_geometry = p_surface_condition->GetGeometry();

            IndexType surface_index = -1;
            for (surface_index = 0; surface_index < number_of_surfaces; ++surface_index) {
                IndexType found_id_count = 0;
                const auto& node_indices = mSurfacesNodeIndexOrdering[surface_index];

                for (const auto& r_surface_node : r_surface_geometry) {
                    for (const auto node_index : node_indices) {
                        found_id_count += (r_element_geometry[node_index].Id() == r_surface_node.Id());
                    }
                }

                if (found_id_count == node_indices.size()) {
                    // found the correct surface
                    break;
                }
            }

            KRATOS_ERROR_IF(surface_index == number_of_surfaces)
                << "No matching surface found for coarse surface condition "
                   "with id "
                << surface_id << " in element " << rCoarseElement.Id() << ".\n";

            r_element_surface_map[surface_id] = surface_index;
        }
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Surface map computed for " << r_model_part.Name() << ".\n";

    KRATOS_CATCH("");
}

void ElementRefinementProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const IndexType number_of_processes = ParallelUtilities::GetNumThreads();
    mThreadLocalStorage.resize(number_of_processes);

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Model Summary [Before]:" << std::endl << mrModel.Info();

    // create thread model parts
    for (IndexType i = 0; i < number_of_processes; ++i) {
        const std::string thread_local_model_part_name = GetThreadLocalModelPartName(i);
        KRATOS_ERROR_IF(mrModel.HasModelPart(thread_local_model_part_name))
            << thread_local_model_part_name
            << " already exists in the model. Please remove it first.\n";
        auto& r_thread_local_model_part = mrModel.CreateModelPart(thread_local_model_part_name);
        r_thread_local_model_part.SetBufferSize(r_model_part.GetBufferSize());
        r_thread_local_model_part.GetNodalSolutionStepVariablesList() = r_model_part.GetNodalSolutionStepVariablesList();
        r_thread_local_model_part.GetProcessInfo() = r_model_part.GetProcessInfo();
        mThreadLocalStorage[i].pRefinedModelPart = &r_thread_local_model_part;
    }

    const IndexType domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

    mReferenceCoarseElement = &r_model_part.Elements().front();
    const IndexType number_of_nodes = mReferenceCoarseElement->GetGeometry().size();

    std::function<void(ThreadLocalStorage&)> execute_refinement;
    if (domain_size == 2) {
        if (number_of_nodes == 3) {
            mSurfacesNodeIndexOrdering = {{0, 1}, {1, 2}, {2, 0}};
            execute_refinement = [&](ThreadLocalStorage& rTLS) {
                IndexType current_node_id{1}, current_condition_id{1}, current_element_id{1};
                BMatrixNM<3, 3> interpolation_matrix = IdentityMatrix(3);
                SurfaceIndicesArray<3> surface_indices({0, 1, 2});
                this->RefineElement<2, 3, 3>(
                    rTLS, mRefinementLevel, current_node_id, current_condition_id,
                    current_element_id, interpolation_matrix, surface_indices);
            };
        } else {
            KRATOS_ERROR
                << "Unsupported geometry type found [ given geometry "
                    "dimension = "
                << domain_size << ", geometry number of nodes = " << number_of_nodes
                << "]. Supported geometries are:\n\t Triangle2D3N\n";
        }
    } else if (domain_size == 3) {
        if (number_of_nodes == 4) {
            // create outer surface model_parts;
            mSurfacesNodeIndexOrdering = {{0, 1, 3}, {3, 1, 2}, {0, 3, 2}, {0, 1, 2}};
            execute_refinement = [&](ThreadLocalStorage& rTLS) {
                IndexType current_node_id{1}, current_condition_id{1}, current_element_id{1};
                BMatrixNM<4, 4> interpolation_matrix = IdentityMatrix(4);
                SurfaceIndicesArray<4> surface_indices({0, 1, 2, 3});
                this->RefineElement<3, 4, 4>(
                    rTLS, mRefinementLevel, current_node_id, current_condition_id,
                    current_element_id, interpolation_matrix, surface_indices);
            };
        } else {
            KRATOS_ERROR
                << "Unsupported geometry type found [ given geometry "
                    "dimension = "
                << domain_size << ", geometry number of nodes = " << number_of_nodes
                << "]. Supported geometries are:\n\t Tetrahedra3D4N\n";
        }
    } else {
        KRATOS_ERROR << "Unsupported domain size. [ domain_size = " << domain_size
                        << " ].\n";
    }

    const auto& dofs = mReferenceCoarseElement->GetGeometry()[0].GetDofs();
    std::vector<const Variable<double>*> dof_vars_list;
    for (const auto& p_dof : dofs) {
        dof_vars_list.push_back(&(KratosComponents<Variable<double>>::Get(p_dof->GetVariable().Name())));
    }

    IndexPartition<IndexType>(number_of_processes).for_each([&](const IndexType) {
        auto& rTLS = GetThreadLocalStorage();

        CreateSurfaceAuxiliaries(rTLS);
        execute_refinement(rTLS);

        IndexType equation_id = 1;
        for (auto& r_node : rTLS.pRefinedModelPart->Nodes()) {
            for (const auto p_dof_var : dof_vars_list) {
                auto& new_dof = r_node.AddDof(*p_dof_var);
                new_dof.SetEquationId(equation_id++);
            }
        }
    });

    ComputeSurfaceMap();

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Model Summary [After]:" << std::endl << mrModel.Info();

    KRATOS_CATCH("");
}

void ElementRefinementProcess::InterpolateAllRefinedMeshesFromCoarseElement(const ElementType& rCoarseElement)
{
    IndexPartition<IndexType>(mThreadLocalStorage.size()).for_each([&](const IndexType CurrentIndex) {
        this->InterpolateThreadLocalRefinedMeshFromCoarseElement(rCoarseElement);
    });
}

void ElementRefinementProcess::InterpolateThreadLocalRefinedMeshFromCoarseElement(const ElementType& rCoarseElement)
{
    KRATOS_TRY

    using namespace ElementRefinementProcessHelperUtilities;

    auto& rTLS = GetThreadLocalStorage();

    const auto& r_geometry = rCoarseElement.GetGeometry();

    for (auto& r_node : rTLS.pRefinedModelPart->Nodes()) {
        r_node.Coordinates() = ZeroVector(3);
        for (IndexType b = 0; b < mHistoricalVariableInterpolationBufferSize; ++b) {
            for (const auto p_variable : mNodalHistoricalVariablesList) {
                r_node.FastGetSolutionStepValue(*p_variable, b) = 0.0;
            }
        }
        for (const auto p_variable : mNodalNonHistoricalVariablesList) {
            r_node.SetValue(*p_variable, 0.0);
        }
        // clears the flags
        r_node.Clear();
    }

    // interpolate nodal data
    for (auto& r_refined_node : rTLS.pRefinedModelPart->Nodes()) {
        const auto& interpolation_values = r_refined_node.GetValue(NODAL_INTERPOLATION_VALUES);
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            const auto& r_coarse_node = r_geometry[i];
            const double interpolation_value = interpolation_values[i];

            // interpolate coordinates
            r_refined_node.Coordinates() += r_coarse_node.Coordinates() * interpolation_value;

            // interpolate historical data
            for (IndexType b = 0; b < mHistoricalVariableInterpolationBufferSize; ++b) {
                for (const auto p_variable : mNodalHistoricalVariablesList) {
                    r_refined_node.FastGetSolutionStepValue(*p_variable, b) += r_coarse_node.FastGetSolutionStepValue(*p_variable, b) * interpolation_value;
                }
            }

            // interpolate non historical data
            for (const auto p_variable : mNodalNonHistoricalVariablesList) {
                r_refined_node.GetValue(*p_variable) += r_coarse_node.GetValue(*p_variable) * interpolation_value;
            }
        }
    }

    // copy element data
    for (auto& r_element : rTLS.pRefinedModelPart->Elements()) {
        CopyEntityData(r_element, rCoarseElement);
    }

    // reset condition data
    for (auto& r_condition : rTLS.pRefinedModelPart->Conditions()) {
        r_condition.Clear();
        r_condition.Data().Clear();
        r_condition.Set(ACTIVE, false);
    }

    const IndexType coarse_element_id = rCoarseElement.Id();
    const auto& r_surface_details = mElementSurfaceIndicesMap[coarse_element_id];

    // copy condition data
    for (auto p_surface_conditions : rCoarseElement.GetValue(NEIGHBOUR_CONDITIONS).GetContainer()) {
        const auto& r_coarse_surface_condition = *p_surface_conditions;
        const IndexType coarse_surface_id = r_coarse_surface_condition.Id();

        const auto p_itr = r_surface_details.find(coarse_surface_id);
        const IndexType surface_index = p_itr->second;
        auto& surface_conditions = rTLS.mRefinedModelPartSurfaceParentElements[surface_index];

        // copy parent element data and flag data
        for (auto& r_pair : surface_conditions) {
            auto& r_condition = *r_pair.first;

            // copy coarse surface condition data to refined condition
            CopyEntityData(r_condition, r_coarse_surface_condition);

            // set the parent element for refined condition properly. This has to be done
            // because we are clearing the data container of r_condition always
            r_condition.SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<ElementType>({r_pair.second}));

            // now set flags on the refined condition nodes if surface condition has defined them
            for (auto& r_node : r_condition.GetGeometry()) {
                for (const auto p_flag : mNodalFlagsList) {
                    if (r_coarse_surface_condition.IsDefined(*p_flag)) {
                        r_node.Set(*p_flag, r_coarse_surface_condition.Is(*p_flag));
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");
}

void ElementRefinementProcess::SetEntityIds(const Variable<int>& rVariable)
{
    KRATOS_TRY

    block_for_each(mThreadLocalStorage, [&](ThreadLocalStorage& rTLS){
        for (auto& r_node : rTLS.pRefinedModelPart->Nodes()) {
            r_node.SetValue(rVariable, r_node.Id());
        }

        for (auto& r_condition : rTLS.pRefinedModelPart->Conditions()) {
            r_condition.SetValue(rVariable, r_condition.Id());
        }

        for (auto& r_element : rTLS.pRefinedModelPart->Elements()) {
            r_element.SetValue(rVariable, r_element.Id());
        }
    });

    KRATOS_CATCH("");
}

void ElementRefinementProcess::SetConditionParentIds(const Variable<int>& rVariable)
{
    block_for_each(mThreadLocalStorage, [&](ThreadLocalStorage& rTLS){
        for (auto& r_condition : rTLS.pRefinedModelPart->Conditions()) {
            if (r_condition.Has(NEIGHBOUR_ELEMENTS)) {
                const auto& neighbors = r_condition.GetValue(NEIGHBOUR_ELEMENTS);

                if (neighbors.size() == 1) {
                    r_condition.SetValue(rVariable, neighbors[0].Id());
                } else {
                    r_condition.SetValue(rVariable, 0);
                }
            } else {
                r_condition.SetValue(rVariable, 0);
            }

        }
    });
}

ModelPart& ElementRefinementProcess::GetThreadLocalModelPart()
{
    return *GetThreadLocalStorage().pRefinedModelPart;
}

std::string ElementRefinementProcess::Info() const
{
    return "ElementRefinementProcess";
}

void ElementRefinementProcess::PrintInfo(std::ostream& rOStream) const
{
}

void ElementRefinementProcess::PrintData(std::ostream& rOStream) const
{
}
} // namespace Kratos