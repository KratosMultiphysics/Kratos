// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part_io.h"
#include "structural_mechanics_application_variables.h"
#include "custom_processes/shell_to_solid_shell_process.h"
#include "custom_processes/solid_shell_thickness_compute_process.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
template<SizeType TNumNodes>
ShellToSolidShellProcess<TNumNodes>::ShellToSolidShellProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
        mThisParameters(ThisParameters)
{
    KRATOS_TRY

    Parameters default_parameters(GetDefaultParameters());

    // Some initial checks
    if (mThisParameters.Has("collapse_geometry")) {
        if (mThisParameters["collapse_geometry"].GetBool()) {
            const std::string element_name = "Element3D" + std::to_string(TNumNodes) + "N";
            if (mThisParameters.Has("element_name")) {
                Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
                if (r_clone_element.GetGeometry().size() != TNumNodes) {
                    mThisParameters["element_name"].SetString(element_name);
                }
            } else {
                default_parameters["element_name"].SetString(element_name);
            }
        }
    }

    mThisParameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::Execute()
{
    // If we compute the extrussion or the inverse extrussion
    if (mThisParameters["collapse_geometry"].GetBool()) {
        ExecuteCollapse();
    } else {
        ExecuteExtrusion();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::ReorderAllIds(const bool ReorderAccordingShellConnectivity)
{
    if (!ReorderAccordingShellConnectivity) {
        NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();
        for(SizeType i = 0; i < r_nodes_array.size(); ++i)
            (it_node_begin + i)->SetId(i + 1);
    } else {
        // The name of the submodelpart
        const std::string& r_model_part_name = mThisParameters["model_part_name"].GetString();
        ModelPart& r_geometry_model_part = r_model_part_name == "" ? mrThisModelPart : mrThisModelPart.GetSubModelPart(r_model_part_name);

        // Auxiliary values
        NodesArrayType& r_nodes_array = r_geometry_model_part.Nodes();
        const SizeType geometry_number_of_nodes = r_nodes_array.size();
        NodesArrayType& total_nodes_array = mrThisModelPart.Nodes();
        const SizeType total_number_of_nodes = total_nodes_array.size();

        // We reoder first all the nodes
        for(SizeType i = 0; i < total_number_of_nodes; ++i)
            (total_nodes_array.begin() + i)->SetId(total_number_of_nodes + i + 1);

        // We reoder now just the shell the nodes
        const auto it_node_begin = r_nodes_array.begin();
        for(SizeType i = 0; i < geometry_number_of_nodes; ++i) {
            auto it_node = it_node_begin + i;
            it_node->SetId(i + 1);
            it_node->Set(VISITED, true);
        }

        // We reoder the rest of all the nodes
        IndexType aux_index = 0;
        for(SizeType i = 0; i < total_number_of_nodes; ++i) {
            auto it_node = total_nodes_array.begin() + i;
            if (it_node->IsNot(VISITED)) {
                it_node->SetId(geometry_number_of_nodes + aux_index + 1);
                aux_index++;
            } else {
                it_node->Set(VISITED, false);
            }
        }
    }

    ConditionsArrayType& condition_array = mrThisModelPart.Conditions();
    for(SizeType i = 0; i < condition_array.size(); ++i)
        (condition_array.begin() + i)->SetId(i + 1);

    ElementsArrayType& element_array = mrThisModelPart.Elements();
    for(SizeType i = 0; i < element_array.size(); ++i)
        (element_array.begin() + i)->SetId(i + 1);

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::ExecuteExtrusion()
{
    KRATOS_TRY

    Model& r_current_model = mrThisModelPart.GetModel();

    // The name of the submodelpart
    const std::string& r_model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_geometry_model_part = r_model_part_name == "" ? mrThisModelPart : mrThisModelPart.GetSubModelPart(r_model_part_name);

    // We initialize some values use later
    NodeType::Pointer p_node_begin = *(r_geometry_model_part.NodesBegin().base());

    // Auxiliary model part where to store new nodes and elements
    ModelPart& r_auxiliary_model_part = r_current_model.CreateModelPart("Extruded" + r_model_part_name);
    const bool create_submodelparts_external_layers = mThisParameters["create_submodelparts_external_layers"].GetBool();
    const bool append_submodelparts_external_layers = mThisParameters["append_submodelparts_external_layers"].GetBool();

    ModelPart& r_auxiliary_model_part_upper = r_current_model.CreateModelPart("AuxiliaryUpper" + r_model_part_name);
    ModelPart& r_auxiliary_model_part_lower = r_current_model.CreateModelPart("AuxiliaryLower" + r_model_part_name);

    // Auxiliary values
    NodesArrayType& r_nodes_array = r_geometry_model_part.Nodes();
    ElementsArrayType& r_elements_array = r_geometry_model_part.Elements();
    const SizeType geometry_number_of_nodes = r_nodes_array.size();
    const SizeType geometry_number_of_elements = r_elements_array.size();
    const SizeType total_number_of_nodes = mrThisModelPart.Nodes().size();
    const SizeType total_number_of_conditions = mrThisModelPart.Conditions().size();
    const SizeType total_number_of_elements = mrThisModelPart.Elements().size();
    const bool replace_previous_geometry = mThisParameters["replace_previous_geometry"].GetBool();

    // First we reoder the ids
    ReorderAllIds(true);

    // We copy the dof from the first node
    const auto it_node_begin = r_nodes_array.begin();
    NodeType::DofsContainerType& dofs = it_node_begin->GetDofs();

    // We initialize the thickness
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;
        it_node->SetValue(THICKNESS, 0.0);
        it_node->SetValue(NODAL_AREA, 0.0);
    }

    // We initialize the r_normal
    const auto it_elem_begin = r_elements_array.begin();
    const array_1d<double, 3> zero_vector = ZeroVector(3);
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(geometry_number_of_elements); ++i)
        (it_elem_begin + i)->SetValue(NORMAL, zero_vector);

    // Calculate the mean of the r_normal in all the nodes
    ComputeNodesMeanNormalModelPartNonHistorical();

    // We compute the nodal thickness
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(geometry_number_of_elements); ++i) {
        auto it_elem = it_elem_begin + i;

        // We get the thickness
        KRATOS_DEBUG_ERROR_IF_NOT(it_elem->GetProperties().Has(THICKNESS)) << "ERROR:: THICKNESS NOT DEFINED" << std::endl;
        const double thickness = it_elem->GetProperties()[THICKNESS];

        auto r_geometry = it_elem->GetGeometry();
        for (IndexType i = 0; i < TNumNodes; ++i) {
            auto& r_node = r_geometry[i];

            double& node_thickness = r_node.GetValue(THICKNESS);

            AtomicAdd(node_thickness, thickness);

            double& nodal_area = r_node.GetValue(NODAL_AREA);

            AtomicAdd(nodal_area, 1.0);
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(geometry_number_of_nodes); ++i) {
        auto it_node = it_node_begin + i;
        double& thickness = it_node->GetValue(THICKNESS);
        thickness /= it_node->GetValue(NODAL_AREA);
    }

    // We create the new nodes
    const SizeType number_of_layers = mThisParameters["number_of_layers"].GetInt();
    for(IndexType i = 0; i < geometry_number_of_nodes; ++i) {
        auto it_node = it_node_begin + i;

        const array_1d<double, 3>& r_normal = it_node->GetValue(NORMAL);
        const double thickness = it_node->GetValue(THICKNESS);
        array_1d<double, 3> coordinates = it_node->Coordinates() - 0.5 * r_normal * thickness;
        const double delta_thickness = thickness/static_cast<double>(number_of_layers);

        IndexType node_id = total_number_of_nodes + i + 1;
        NodeType::Pointer p_node0 = r_auxiliary_model_part.CreateNewNode(node_id, coordinates[0], coordinates[1], coordinates[2]);

        // Set the DOFs in the nodes
        for (auto it_dof = dofs.begin(); it_dof != dofs.end(); ++it_dof)
            p_node0->pAddDof(**it_dof);

        // We copy the step data
        CopyVariablesList(p_node0, p_node_begin);

        // We add the node to the external layers
        if (create_submodelparts_external_layers) r_auxiliary_model_part_lower.AddNode(p_node0);

        for (IndexType j = 0; j < number_of_layers; ++j) {
            coordinates += r_normal * delta_thickness;
            node_id = (j + 1) * geometry_number_of_nodes + total_number_of_nodes + i + 1;
            NodeType::Pointer p_node1 = r_auxiliary_model_part.CreateNewNode(node_id, coordinates[0], coordinates[1], coordinates[2]);

            // Set the DOFs in the nodes
            for (auto it_dof = dofs.begin(); it_dof != dofs.end(); ++it_dof)
                p_node1->pAddDof(**it_dof);

            // We copy the step data
            CopyVariablesList(p_node1, p_node_begin);

            // We add the node to the external layers
            if (create_submodelparts_external_layers && j == (number_of_layers - 1)) r_auxiliary_model_part_upper.AddNode(p_node1);
        }

        // We set the flag TO_ERASE for later remove the nodes
        it_node->Set(TO_ERASE, replace_previous_geometry);
    }

    const std::string& element_name = mThisParameters["element_name"].GetString();
    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    KRATOS_ERROR_IF_NOT(r_clone_element.GetGeometry().size() == 2 * TNumNodes) << "ERROR: Element " << element_name << " has a different number of nodes to " << 2 * TNumNodes << std::endl;

    // We will save the list of properties ids to later set the CL
    std::unordered_set<IndexType> set_id_properties;

    // We create the new elements
    IndexType element_counter = total_number_of_elements;
    for(IndexType i = 0; i < geometry_number_of_elements; ++i) {
        auto it_elem = r_elements_array.begin() + i;

        auto p_prop = it_elem->pGetProperties();
        set_id_properties.insert(p_prop->Id());
        for (IndexType j = 0; j < number_of_layers; ++j) {
            std::vector<IndexType> element_node_ids (2 * TNumNodes);
            for (IndexType k = 0; k < TNumNodes; ++k) {
                const IndexType index_node = it_elem->GetGeometry()[k].Id();
                element_node_ids[k] = total_number_of_nodes + index_node + j * geometry_number_of_nodes;
                element_node_ids[k + TNumNodes] = total_number_of_nodes + index_node + (j + 1) * geometry_number_of_nodes;
            }
            element_counter++;
            r_auxiliary_model_part.CreateNewElement(element_name, element_counter, element_node_ids, p_prop);
        }

        // We set the flag TO_ERASE for later remove the elements
        it_elem->Set(TO_ERASE, replace_previous_geometry);
    }

    // We create the conditions if necessary
    if (create_submodelparts_external_layers) {
        const std::string condition_name = "SurfaceCondition3D" + std::to_string(TNumNodes) + "N";
        IndexType condition_counter = total_number_of_conditions;
        for(IndexType i = 0; i < geometry_number_of_elements; ++i) {
            auto it_elem = r_elements_array.begin() + i;

            auto p_prop = it_elem->pGetProperties();
            set_id_properties.insert(p_prop->Id());

            std::vector<IndexType> upper_condition_node_ids (TNumNodes);
            std::vector<IndexType> lower_condition_node_ids (TNumNodes);
            for (IndexType k = 0; k < TNumNodes; ++k) {
                const IndexType index_node = it_elem->GetGeometry()[k].Id();
                lower_condition_node_ids[TNumNodes - (k + 1)] = total_number_of_nodes + index_node; // We invert the order for the lower face to have the normal in the right direction
                upper_condition_node_ids[k] = total_number_of_nodes + index_node + number_of_layers * geometry_number_of_nodes;
            }
            condition_counter++;
            r_auxiliary_model_part_lower.CreateNewCondition(condition_name, condition_counter, lower_condition_node_ids, p_prop);
            condition_counter++;
            r_auxiliary_model_part_upper.CreateNewCondition(condition_name, condition_counter, upper_condition_node_ids, p_prop);
        }
    }

    // We reassign a new constitutive law
    ReassignConstitutiveLaw(r_geometry_model_part, set_id_properties);

    // In case we replace the geometry
    if (replace_previous_geometry) {
        ReplacePreviousGeometry(r_geometry_model_part, r_auxiliary_model_part);
    }

    // We copy the external layers
    if (create_submodelparts_external_layers) {
        const std::string name_upper = "Upper_" + r_model_part_name;
        ModelPart& r_upper_model_part = append_submodelparts_external_layers ? r_geometry_model_part.CreateSubModelPart(name_upper) :  mrThisModelPart.CreateSubModelPart(name_upper);
        r_upper_model_part.AddNodes( r_auxiliary_model_part_upper.NodesBegin(), r_auxiliary_model_part_upper.NodesEnd() );
        r_upper_model_part.AddConditions( r_auxiliary_model_part_upper.ConditionsBegin(), r_auxiliary_model_part_upper.ConditionsEnd() );
        const std::string name_lower = "Lower_" + r_model_part_name;
        ModelPart& r_lower_model_part = append_submodelparts_external_layers ? r_geometry_model_part.CreateSubModelPart(name_lower) : mrThisModelPart.CreateSubModelPart(name_lower);
        r_lower_model_part.AddNodes( r_auxiliary_model_part_lower.NodesBegin(), r_auxiliary_model_part_lower.NodesEnd() );
        r_lower_model_part.AddConditions( r_auxiliary_model_part_lower.ConditionsBegin(), r_auxiliary_model_part_lower.ConditionsEnd() );
    }

    // We add to the computing model part if available
    const std::string& computing_model_part_name = mThisParameters["computing_model_part_name"].GetString();
    if (computing_model_part_name != "") {
        ModelPart& r_computing_model_part = mrThisModelPart.GetSubModelPart(computing_model_part_name);
        r_computing_model_part.AddNodes( r_auxiliary_model_part.NodesBegin(), r_auxiliary_model_part.NodesEnd() );
        r_computing_model_part.AddElements( r_auxiliary_model_part.ElementsBegin(), r_auxiliary_model_part.ElementsEnd() );
        r_computing_model_part.AddConditions( r_auxiliary_model_part_upper.ConditionsBegin(), r_auxiliary_model_part_upper.ConditionsEnd() );
        r_computing_model_part.AddConditions( r_auxiliary_model_part_lower.ConditionsBegin(), r_auxiliary_model_part_lower.ConditionsEnd() );
    }

    // Reorder again all the IDs
    ReorderAllIds();

    // We initialize the new elements
    if (mThisParameters["initialize_elements"].GetBool()) {
        InitializeElements();
    }

    // Export to *.mdpa if desired
    if (mThisParameters["export_to_mdpa"].GetBool()) {
        ExportToMDPA();
    }

    // Clean the model
    CleanModel();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::ExecuteCollapse()
{
    KRATOS_TRY

    Model& r_current_model = mrThisModelPart.GetModel();

    // The name of the submodelpart
    const std::string& r_model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_geometry_model_part = r_model_part_name == "" ? mrThisModelPart : mrThisModelPart.GetSubModelPart(r_model_part_name);

    // We initialize some values use later
    NodeType::Pointer p_node_begin = *(r_geometry_model_part.NodesBegin().base());

    // Auxiliary model part where to store new nodes and elements
    ModelPart& r_auxiliary_model_part = r_current_model.CreateModelPart("Collapsed" + r_model_part_name);

    // We compute the thickness
    SolidShellThickComputeProcess thickness_process(r_geometry_model_part);
    thickness_process.Execute();

    // Auxiliary values
    NodesArrayType& r_nodes_array = r_geometry_model_part.Nodes();
    ElementsArrayType& r_elements_array = r_geometry_model_part.Elements();
    const SizeType geometry_number_of_elements = r_elements_array.size();
    const SizeType total_number_of_nodes = mrThisModelPart.Nodes().size();
    const SizeType total_number_of_elements = mrThisModelPart.Elements().size();
    const bool replace_previous_geometry = mThisParameters["replace_previous_geometry"].GetBool();

    // First we reoder the ids
    ReorderAllIds(true);

    // We copy the dof from the first node
    const auto it_node_begin = r_nodes_array.begin();
    NodeType::DofsContainerType& dofs = it_node_begin->GetDofs();

    // Initial check
    const SizeType number_of_layers = mThisParameters["number_of_layers"].GetInt();
    KRATOS_ERROR_IF(number_of_layers > 1) << "Collapsed only compatible with one layer" << std::endl;

    // We get the element name
    const std::string& element_name = mThisParameters["element_name"].GetString();
    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    KRATOS_ERROR_IF_NOT(r_clone_element.GetGeometry().size() == TNumNodes) << "ERROR: Element " << element_name << " has a different number of nodes to " << TNumNodes << std::endl;

    // We will save the list of properties ids to later set the CL
    std::unordered_set<IndexType> set_id_properties;

    // We create the new nodes and elements
    const auto it_elem_begin = r_elements_array.begin();
    IndexType node_id = total_number_of_nodes + 1;
    IndexType element_id = total_number_of_elements + 1;
    std::vector<IndexType> element_node_ids (TNumNodes);
    for(IndexType i = 0; i < geometry_number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;

        // Getting properties and geometry
        auto p_prop = it_elem->pGetProperties();
        auto& r_geometry = it_elem->GetGeometry();

        // Adding property
        set_id_properties.insert(p_prop->Id());

        // Collapsing nodes
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3> coordinates = 0.5 * r_geometry[i_node].Coordinates() + 0.5 * r_geometry[i_node + TNumNodes].Coordinates();
            NodeType::Pointer p_node = r_auxiliary_model_part.CreateNewNode(node_id, coordinates[0], coordinates[1], coordinates[2]);
            p_node->SetValue(THICKNESS, r_geometry[i_node].GetValue(THICKNESS));

            element_node_ids[i_node] = node_id;

            // Set the DOFs in the nodes
            for (auto it_dof = dofs.begin(); it_dof != dofs.end(); ++it_dof)
                p_node->pAddDof(**it_dof);

            // We copy the step data
            CopyVariablesList(p_node, p_node_begin);

            ++node_id;
        }

        // Now we create the new elements
        r_auxiliary_model_part.CreateNewElement(element_name, element_id, element_node_ids, p_prop);

        // We set the flag TO_ERASE for later remove the nodes and elements
        it_elem->Set(TO_ERASE, replace_previous_geometry);
        for (auto& r_node : r_geometry) {
            r_node.Set(TO_ERASE, replace_previous_geometry);
        }

        ++element_id;
    }

    // We reassign a new constitutive law
    ReassignConstitutiveLaw(r_geometry_model_part, set_id_properties);

    // In case we replace the geometry
    if (replace_previous_geometry) {
        ReplacePreviousGeometry(r_geometry_model_part, r_auxiliary_model_part);
    }

    // We add to the computing model part if available
    const std::string& computing_model_part_name = mThisParameters["computing_model_part_name"].GetString();
    if (computing_model_part_name != "") {
        ModelPart& r_computing_model_part = mrThisModelPart.GetSubModelPart(computing_model_part_name);
        r_computing_model_part.AddNodes( r_auxiliary_model_part.NodesBegin(), r_auxiliary_model_part.NodesEnd() );
        r_computing_model_part.AddElements( r_auxiliary_model_part.ElementsBegin(), r_auxiliary_model_part.ElementsEnd() );
    }

    // Reorder again all the IDs
    ReorderAllIds();

    // We initialize the new elements
    if (mThisParameters["initialize_elements"].GetBool()) {
        InitializeElements();
    }

    // Export to *.mdpa if desired
    if (mThisParameters["export_to_mdpa"].GetBool()) {
        ExportToMDPA();
    }

    // Clean the model
    CleanModel();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::ReplacePreviousGeometry(
    ModelPart& rGeometryModelPart,
    ModelPart& rAuxiliaryModelPart
    )
{
    // Finally we remove the old nodes and elements
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // We copy the new model part to the original one
    rGeometryModelPart.AddNodes( rAuxiliaryModelPart.NodesBegin(), rAuxiliaryModelPart.NodesEnd() );
    rGeometryModelPart.AddElements( rAuxiliaryModelPart.ElementsBegin(), rAuxiliaryModelPart.ElementsEnd() );
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::ReassignConstitutiveLaw(
    ModelPart& rGeometryModelPart,
    std::unordered_set<IndexType>& rSetIdProperties
    )
{
    const std::string& new_constitutive_law_name = mThisParameters["new_constitutive_law_name"].GetString();
    if (new_constitutive_law_name != "") {
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get(new_constitutive_law_name).Clone();
        for (auto id_prop : rSetIdProperties) {
            auto p_prop = rGeometryModelPart.pGetProperties(id_prop);
            p_prop->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::InitializeElements()
{
    ElementsArrayType& element_array = mrThisModelPart.Elements();
    const auto& r_process_info = mrThisModelPart.GetProcessInfo();
    for(SizeType i = 0; i < element_array.size(); ++i)
        (element_array.begin() + i)->Initialize(r_process_info);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::ExportToMDPA()
{
    const std::string& output_name = mThisParameters["output_name"].GetString();
    std::ofstream output_file;
    ModelPartIO model_part_io(output_name, IO::WRITE);
    model_part_io.WriteModelPart(mrThisModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::CleanModel()
{
    // The original model part name
    const std::string& r_model_part_name = mThisParameters["model_part_name"].GetString();

    // If we remove the Extruded/Collapsed model part
    const bool replace_previous_geometry = mThisParameters["replace_previous_geometry"].GetBool();
    const bool collapse_geometry = mThisParameters["collapse_geometry"].GetBool();

    // Getting the model
    Model& r_current_model = mrThisModelPart.GetModel();

    // Removing model parts
    if (replace_previous_geometry) {
        if (collapse_geometry) {
            r_current_model.DeleteModelPart("Collapsed" + r_model_part_name);
        } else {
            r_current_model.DeleteModelPart("Extruded" + r_model_part_name);
        }
    }
    r_current_model.DeleteModelPart("AuxiliaryUpper"  + r_model_part_name);
    r_current_model.DeleteModelPart("AuxiliaryLower"  + r_model_part_name);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
inline void ShellToSolidShellProcess<TNumNodes>::ComputeNodesMeanNormalModelPartNonHistorical()
{
    // Tolerance
    const double tolerance = std::numeric_limits<double>::epsilon();

    // The name of the submodelpart
    const std::string& r_model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_geometry_model_part = r_model_part_name == "" ? mrThisModelPart : mrThisModelPart.GetSubModelPart(r_model_part_name);

    // We iterate over the nodes
    NodesArrayType& r_nodes_array = r_geometry_model_part.Nodes();
    const int num_nodes = static_cast<int>(r_nodes_array.size());

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i)
        (r_nodes_array.begin() + i)->SetValue(NORMAL, ZeroVector(3));

    // Sum all the nodes normals
    ElementsArrayType& r_elements_array = r_geometry_model_part.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = r_elements_array.begin() + i;
        GeometryType& this_geometry = it_elem->GetGeometry();

        // Aux coordinates
        array_1d<double, 3> aux_coords;
        aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_geometry.Center());

        it_elem->SetValue(NORMAL, this_geometry.UnitNormal(aux_coords));

        const unsigned int number_nodes = this_geometry.PointsNumber();

        for (unsigned int i = 0; i < number_nodes; ++i) {
            auto& this_node = this_geometry[i];
            aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_node.Coordinates());
            const array_1d<double, 3>& r_normal = this_geometry.UnitNormal(aux_coords);
            array_1d<double, 3>& aux_normal = this_node.GetValue(NORMAL);
            AtomicAdd(aux_normal, r_normal);
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) {
        auto it_node = r_nodes_array.begin() + i;

        array_1d<double, 3>& r_normal = it_node->GetValue(NORMAL);
        const double norm_normal = norm_2(r_normal);
        if (norm_normal > tolerance) r_normal /= norm_normal;
        else KRATOS_ERROR << "WARNING:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
inline void ShellToSolidShellProcess<TNumNodes>::CopyVariablesList(
    NodeType::Pointer pNodeNew,
    NodeType::Pointer pNodeOld
    )
{
    auto& node_data = pNodeNew->SolutionStepData();
    auto& node_data_reference = pNodeOld->SolutionStepData();
    node_data.SetVariablesList(node_data_reference.pGetVariablesList());
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
const Parameters ShellToSolidShellProcess<TNumNodes>::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "element_name"                         : "SolidShellElementSprism3D6N",
        "new_constitutive_law_name"            : "",
        "model_part_name"                      : "",
        "number_of_layers"                     : 1,
        "export_to_mdpa"                       : false,
        "output_name"                          : "output",
        "computing_model_part_name"            : "computing_domain",
        "create_submodelparts_external_layers" : false,
        "append_submodelparts_external_layers" : false,
        "initialize_elements"                  : false,
        "replace_previous_geometry"            : true,
        "collapse_geometry"                    : false
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class ShellToSolidShellProcess<3>;
template class ShellToSolidShellProcess<4>;
// class ShellToSolidShellProcess
} // namespace Kratos
