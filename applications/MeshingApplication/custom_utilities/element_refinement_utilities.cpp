// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/data_communicator.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "meshing_application_variables.h"

// Include base h
#include "element_refinement_utilities.h"

namespace Kratos
{

namespace ElementRefinementHelperUtilities
{

using IndexType = std::size_t;

template<class TContainerType>
IndexType FindMaxId(
    const TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const IndexType local_max_id = IndexPartition<IndexType>(rContainer.size()).for_each<MaxReduction<IndexType>>([&](const IndexType Index) {
        return (rContainer.begin() + Index)->Id();
    });

    return rDataCommunicator.MaxAll(local_max_id);

    KRATOS_CATCH("");
}

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

    block_for_each(rModelPart.Elements(), [](Element& rElement) {
        rElement.SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>());
    });

    block_for_each(rModelPart.Conditions(), [](Condition& rCondition) {
        if (rCondition.Has(NEIGHBOUR_ELEMENTS)) {
            auto& parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];

#pragma omp critical
            {
                parent_element.GetValue(NEIGHBOUR_CONDITIONS).push_back(GlobalPointer<Condition>(&rCondition));
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
template IndexType FindMaxId<ModelPart::NodesContainerType>(const ModelPart::NodesContainerType&, const DataCommunicator&);
template IndexType FindMaxId<ModelPart::ConditionsContainerType>(const ModelPart::ConditionsContainerType&, const DataCommunicator&);
template IndexType FindMaxId<ModelPart::ElementsContainerType>(const ModelPart::ElementsContainerType&, const DataCommunicator&);

template void CopyEntityData<ModelPart::ConditionType>(ModelPart::ConditionType&, const ModelPart::ConditionType&);
template void CopyEntityData<ModelPart::ElementType>(ModelPart::ElementType&, const ModelPart::ElementType&);

}

template <>
void ElementRefinementUtilities::RefineElement<2, 3, 3>(
    const IndexType RefinementLevel,
    const BMatrixNN<3>& rNodalInterpolationValues,
    const SurfaceIndicesArray<3>& rSurfaceIndices)
{
    KRATOS_TRY

    if (RefinementLevel == 0) {
        // add all refined elements, conditions and nodes
        const auto& r_geometry = mrReferenceElement.GetGeometry();

        // interpolate nodal locations
        BMatrixNN<3> reference_coordinates, interpolated_coordinates;
        row(reference_coordinates, 0) = r_geometry[0].Coordinates();
        row(reference_coordinates, 1) = r_geometry[1].Coordinates();
        row(reference_coordinates, 2) = r_geometry[2].Coordinates();

        noalias(interpolated_coordinates) = prod(rNodalInterpolationValues, reference_coordinates);

        // create nodes
        std::vector<IndexType> node_ids(3, 0);
        for (IndexType i = 0; i < 3; ++i) {
            auto& r_current_id = node_ids[i];
            const Array3D& current_coordinates = row(interpolated_coordinates, i);

            for (auto& r_node : mpRefinedModelPart->Nodes()) {
                if (norm_2(r_node.Coordinates() - current_coordinates) <= std::numeric_limits<double>::epsilon()) {
                    r_current_id = r_node.Id();
                    break;
                }
            }

            if (r_current_id == 0) {
                r_current_id = mStartNodeId++;
                auto p_node = mpRefinedModelPart->CreateNewNode(r_current_id, current_coordinates[0], current_coordinates[1], current_coordinates[2]);
                p_node->SetValue(NODAL_INTERPOLATION_VALUES, row(rNodalInterpolationValues, i));
            }
        }

        auto p_properties = mrReferenceElement.pGetProperties();

        // create elements
        auto p_element = &*(mpRefinedModelPart->CreateNewElement(mrRefinedElementName, mStartElementId++, node_ids, p_properties));

        // create surface conditions
        for (IndexType i = 0; i < 3; ++i) {
            const IndexType current_surface_index = rSurfaceIndices[i];
            if (current_surface_index < 3) {
                const auto& condition_node_ids = std::vector<IndexType>{node_ids[i % 3], node_ids[(i + 1) % 3]};
                auto p_surface_condition = mpRefinedModelPart->CreateNewCondition(
                    mrRefinedConditionName,
                    mStartConditionId++,
                    condition_node_ids, p_properties);
                mRefinedModelPartSurfaceParentElements[current_surface_index].push_back(std::make_pair(&(*p_surface_condition), p_element));
                mSurfaceModelParts[current_surface_index]->AddCondition(p_surface_condition);
                mSurfaceModelParts[current_surface_index]->AddNodes(condition_node_ids);
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
        BMatrixNN<3> sub_triangle_interpolation_values;

        // triangle 1
        row(sub_triangle_interpolation_values, 0) = row(rNodalInterpolationValues, 0);
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 2)) * 0.5;
        this->RefineElement<2, 3, 3>(RefinementLevel - 1, sub_triangle_interpolation_values, SurfaceIndicesArray<3>({rSurfaceIndices[0], 4, rSurfaceIndices[2]}));

        // triangle 2
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 1)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = row(rNodalInterpolationValues, 1);
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        this->RefineElement<2, 3, 3>(RefinementLevel - 1, sub_triangle_interpolation_values, SurfaceIndicesArray<3>({rSurfaceIndices[0], rSurfaceIndices[1], 4}));

        // triangle 3
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = row(rNodalInterpolationValues, 2);
        this->RefineElement<2, 3, 3>(RefinementLevel - 1, sub_triangle_interpolation_values, SurfaceIndicesArray<3>({4, rSurfaceIndices[1], rSurfaceIndices[2]}));

        // triangle 4
        row(sub_triangle_interpolation_values, 0) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 1) = (row(rNodalInterpolationValues, 0) + row(rNodalInterpolationValues, 2)) * 0.5;
        row(sub_triangle_interpolation_values, 2) = (row(rNodalInterpolationValues, 1) + row(rNodalInterpolationValues, 0)) * 0.5;
        this->RefineElement<2, 3, 3>(RefinementLevel - 1, sub_triangle_interpolation_values, SurfaceIndicesArray<3>({4, 4, 4}));
    }

    KRATOS_CATCH("");
}

template <>
void ElementRefinementUtilities::RefineElement<3, 4, 4>(
    const IndexType RefinementLevel,
    const BMatrixNN<4>& rNodalInterpolationValues,
    const SurfaceIndicesArray<4>& rSurfaceIndices)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

void ElementRefinementUtilities::InterpolateToRefinedMeshFromCoarseElement(
    const ElementType& rCoarseElement)
{
    KRATOS_TRY

    using namespace ElementRefinementHelperUtilities;

    const auto& r_geometry = rCoarseElement.GetGeometry();
    const IndexType coarse_element_number_of_nodes = r_geometry.size();

    for (auto& r_node : mpRefinedModelPart->Nodes()) {
        r_node.Coordinates() = ZeroVector(3);
        for (const auto p_variable : mNodalHistoricalVariablesList) {
            r_node.FastGetSolutionStepValue(*p_variable) = 0.0;
        }
        for (const auto p_variable : mNodalNonHistoricalVariablesList) {
            r_node.SetValue(*p_variable, 0.0);
        }
        // clears the flags
        r_node.Clear();
    }

    // interpolate nodal data
    for (auto& r_refined_node : mpRefinedModelPart->Nodes()) {
        const auto& interpolation_values = r_refined_node.GetValue(NODAL_INTERPOLATION_VALUES);
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            const auto& r_coarse_node = r_geometry[i];
            const double interpolation_value = interpolation_values[i];

            // interpolate coordinates
            r_refined_node.Coordinates() += r_coarse_node.Coordinates() * interpolation_value;

            // interpolate historical data
            for (const auto p_variable : mNodalHistoricalVariablesList) {
                r_refined_node.FastGetSolutionStepValue(*p_variable) += r_coarse_node.FastGetSolutionStepValue(*p_variable) * interpolation_value;
            }

            // interpolate non historical data
            for (const auto p_variable : mNodalNonHistoricalVariablesList) {
                r_refined_node.GetValue(*p_variable) += r_coarse_node.GetValue(*p_variable) * interpolation_value;
            }
        }
    }

    // copy element data
    for (auto& r_element : mpRefinedModelPart->Elements()) {
        CopyEntityData(r_element, rCoarseElement);
    }

    // reset condition data
    for (auto& r_condition : mpRefinedModelPart->Conditions()) {
        r_condition.Clear();
        r_condition.Data().Clear();
        r_condition.Set(ACTIVE, false);
    }

    // copy condition data
    for (auto p_surface_conditions : rCoarseElement.GetValue(NEIGHBOUR_CONDITIONS).GetContainer()) {
        const auto& r_surface_condition = *p_surface_conditions;
        const auto& r_surface_geometry = r_surface_condition.GetGeometry();

        const IndexType start_id = r_surface_geometry[0].Id();
        bool is_reversed = false;

        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            const IndexType coarse_element_node_id = r_geometry[i].Id();

            // found starting node of the condition
            if (start_id == coarse_element_node_id) {
                for (IndexType j = 1; j < r_surface_geometry.size(); ++j) {
                    if (r_surface_geometry[j].Id() != r_geometry[(i + j) % coarse_element_number_of_nodes].Id()) {
                        is_reversed = true;
                        break;
                    }
                }

                const IndexType surface_index =
                    (mRefinedModelPartSurfaceParentElements.size() - is_reversed + i) %
                    mRefinedModelPartSurfaceParentElements.size();

                auto& surface_conditions = mRefinedModelPartSurfaceParentElements[surface_index];

                // copy data
                for (auto& r_pair : surface_conditions) {
                    CopyEntityData(*(r_pair.first), r_surface_condition);
                }

                // reverse the order of nodes if reversed coarse condition is found
                if (is_reversed) {
                    for (auto& r_pair : surface_conditions) {
                        auto& r_geometry = r_pair.first->GetGeometry();
                        for (IndexType i = 0; i < r_geometry.size(); ++i) {
                            r_geometry(i).swap(r_geometry(r_geometry.PointsNumber() - i));
                        }
                    }
                }

                // copy parent element data and flag data
                for (auto& r_pair : surface_conditions) {
                    auto& r_condition = *r_pair.first;
                    r_condition.SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<ElementType>({r_pair.second}));
                    for (auto& r_node : r_condition.GetGeometry()) {
                        for (const auto p_flag : mNodalFlagsList) {
                            if (r_surface_condition.IsDefined(*p_flag)) {
                                r_node.Set(*p_flag, r_surface_condition.Is(*p_flag));
                            }
                        }
                    }
                }

                break;
            }
        }
    }

    KRATOS_CATCH("");
}

ElementRefinementUtilities::ElementRefinementUtilities(
    ModelPart& rModelPart,
    const std::string& rRefinedModelPartName,
    const ElementType& rReferenceElement,
    const IndexType RefinementLevel,
    const std::string rRefinedElementName,
    const std::string rRefinedConditionName,
    const std::vector<std::string>& rHistoricalVariableNamesList,
    const std::vector<std::string>& rNonHistoricalVariableNamesList,
    const std::vector<std::string>& rNodalFlagNamesList)
    : mrReferenceElement(rReferenceElement),
      mRefinementLevel(RefinementLevel),
      mrRefinedElementName(rRefinedElementName),
      mrRefinedConditionName(rRefinedConditionName)
{
    KRATOS_TRY

    using namespace ElementRefinementHelperUtilities;

    FindElementSurfaceConditions(rModelPart);

    // create new local model part and fill it with elements, conditions and nodes
    auto& r_model = rModelPart.GetModel();

    KRATOS_ERROR_IF(r_model.HasModelPart(rRefinedModelPartName))
        << rRefinedModelPartName << " already exists in the model. Please remove it first.\n";

    mStartNodeId = 1; //FindMaxId(r_parent_model_part.Nodes(), r_data_communicator) + 1;
    mStartConditionId = 1; //FindMaxId(r_parent_model_part.Conditions(), r_data_communicator) + 1;
    mStartElementId = 1; //FindMaxId(r_parent_model_part.Elements(), r_data_communicator) + 1;

    mpRefinedModelPart = &r_model.CreateModelPart(rRefinedModelPartName);
    mpRefinedModelPart->GetNodalSolutionStepVariablesList() = rModelPart.GetNodalSolutionStepVariablesList();
    mpRefinedModelPart->GetProcessInfo() = rModelPart.GetProcessInfo();

    const int domain_size = mpRefinedModelPart->GetProcessInfo()[DOMAIN_SIZE];
    const IndexType number_of_nodes = mrReferenceElement.GetGeometry().size();

    if (domain_size == 2) {
        if (number_of_nodes == 3) {
            // create outer surface model_parts;
            CreateSurfaceAuxiliaries(3);
            BMatrixNN<3> interpolation_matrix = IdentityMatrix(3);
            SurfaceIndicesArray<3> surface_indices({0, 1, 2});
            this->RefineElement<2, 3, 3>(RefinementLevel, interpolation_matrix, surface_indices);
        } else {
            KRATOS_ERROR << "Unsupported geometry with " << number_of_nodes << " nodes in 2D provided.\n";
        }
    } else if (domain_size == 3) {
        if (number_of_nodes == 4) {
            // create outer surface model_parts;
            CreateSurfaceAuxiliaries(4);
            BMatrixNN<4> interpolation_matrix = IdentityMatrix(4);
            SurfaceIndicesArray<4> surface_indices({0, 1, 2, 3});
            this->RefineElement<3, 4, 4>(RefinementLevel, interpolation_matrix, surface_indices);
        } else {
            KRATOS_ERROR << "Unsupported geometry with " << number_of_nodes << " nodes in 3D provided.\n";
        }
    } else {
        KRATOS_ERROR << "Unsupported domain size provided. [ domain_size = " << domain_size << " ].\n";
    }

    const auto& dofs = rReferenceElement.GetGeometry()[0].GetDofs();

    std::vector<const Variable<double>*> dof_vars_list;
    for (const auto& p_dof : dofs) {
        dof_vars_list.push_back(&(KratosComponents<Variable<double>>::Get(p_dof->GetVariable().Name())));
    }

    IndexType equation_id = 1;
    for (auto& r_node : mpRefinedModelPart->Nodes()) {
        for (const auto p_dof_var : dof_vars_list) {
            auto& new_dof = r_node.AddDof(*p_dof_var);
            new_dof.SetEquationId(equation_id++);
        }
    }

    for (const auto& r_var_name : rHistoricalVariableNamesList) {
        AddVariableToList(mNodalHistoricalVariablesList, r_var_name);
    }

    for (const auto p_var : mNodalHistoricalVariablesList) {
        KRATOS_ERROR_IF_NOT(mpRefinedModelPart->HasNodalSolutionStepVariable(*p_var))
            << p_var->Name() << " not found in nodal historical variables list in "
            << rModelPart.Name() << ".\n";
    }

    for (const auto& r_var_name : rNonHistoricalVariableNamesList) {
        AddVariableToList(mNodalNonHistoricalVariablesList, r_var_name);
    }

    for (const auto& r_flag_name : rNodalFlagNamesList) {
        if (KratosComponents<Flags>::Has(r_flag_name)) {
            mNodalFlagsList.push_back(&KratosComponents<Flags>::Get(r_flag_name));
        } else {
            KRATOS_ERROR << "Flag name not defined. [ flag name = " << r_flag_name << " ].\n";
        }
    }

    KRATOS_CATCH("");
}


void ElementRefinementUtilities::SetIds(const Variable<int>& rVariable)
{
    KRATOS_TRY

    for (auto& r_node : mpRefinedModelPart->Nodes()) {
        r_node.SetValue(rVariable, r_node.Id());
    }

    for (auto& r_condition : mpRefinedModelPart->Conditions()) {
        r_condition.SetValue(rVariable, r_condition.Id());
    }

    for (auto& r_element : mpRefinedModelPart->Elements()) {
        r_element.SetValue(rVariable, r_element.Id());
    }

    KRATOS_CATCH("");
}

void ElementRefinementUtilities::SetConditionParentIds(const Variable<int>& rVariable)
{
    KRATOS_TRY

    for (auto& r_condition : mpRefinedModelPart->Conditions()) {
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

    KRATOS_CATCH("");
}

void ElementRefinementUtilities::CreateSurfaceAuxiliaries(const IndexType NumberOfSurfaces)
{
    KRATOS_TRY

    mRefinedModelPartSurfaceParentElements.resize(NumberOfSurfaces);
    for (IndexType i = 0; i < NumberOfSurfaces; ++i) {
        mSurfaceModelParts.push_back(&mpRefinedModelPart->CreateSubModelPart("Surface_" + std::to_string(i)));
    }

    KRATOS_CATCH("");
}

ModelPart& ElementRefinementUtilities::GetRefinedModelPart()
{
    return *mpRefinedModelPart;
}

// template instantiations

template void ElementRefinementUtilities::RefineElement<2, 3, 3>(const std::size_t, const BoundedMatrix<double, 3, 3>&, const array_1d<std::size_t, 3>&);
template void ElementRefinementUtilities::RefineElement<3, 4, 4>(const std::size_t, const BoundedMatrix<double, 4, 4>&, const array_1d<std::size_t, 4>&);

} // namespace Kratos
