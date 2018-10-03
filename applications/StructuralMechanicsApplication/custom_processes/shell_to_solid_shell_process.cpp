// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "includes/kernel.h"
#include "containers/model.h"
#include "custom_processes/shell_to_solid_shell_process.h"
#include "structural_mechanics_application_variables.h"
#include "includes/model_part_io.h"

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

    Parameters default_parameters = Parameters(R"(
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
        "initialize_elements"                  : false
    })" );

    mThisParameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::Execute()
{
    KRATOS_TRY

    Model& current_model = mrThisModelPart.GetOwnerModel();

    // The name of the submodelpart
    const std::string& model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& geometry_model_part = model_part_name == "" ? mrThisModelPart : mrThisModelPart.GetSubModelPart(model_part_name);

    // We initialize some values use later
    NodeType::Pointer p_node_begin = *(geometry_model_part.NodesBegin().base());

    // Auxiliar model part where to store new nodes and elements
    ModelPart& auxiliar_model_part = current_model.CreateModelPart("AuxiliarModelPart");
    const bool create_submodelparts_external_layers = mThisParameters["create_submodelparts_external_layers"].GetBool();
    const bool append_submodelparts_external_layers = mThisParameters["append_submodelparts_external_layers"].GetBool();

    
    ModelPart& auxiliar_model_part_upper = current_model.CreateModelPart("upper");
    ModelPart& auxiliar_model_part_lower = current_model.CreateModelPart("lower");

    // Auxiliar values
    NodesArrayType& nodes_array = geometry_model_part.Nodes();
    ElementsArrayType& elements_array = geometry_model_part.Elements();
    const SizeType geometry_number_of_nodes = nodes_array.size();
    const SizeType geometry_number_of_elements = elements_array.size();
    const SizeType total_number_of_nodes = mrThisModelPart.Nodes().size();
    const SizeType total_number_of_conditions = mrThisModelPart.Conditions().size();
    const SizeType total_number_of_elements = mrThisModelPart.Elements().size();

    // First we reoder the ids
    ReorderAllIds(true);

    // We copy the dof from the first node
    NodeType::DofsContainerType dofs = nodes_array.begin()->GetDofs();

    // We initialize the thickness
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->SetValue(THICKNESS, 0.0);
        it_node->SetValue(NODAL_AREA, 0.0);
    }

    // We initialize the normal
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(geometry_number_of_elements); ++i)
        (elements_array.begin() + i)->SetValue(NORMAL, ZeroVector(3));

    // Calculate the mean of the normal in all the nodes
    ComputeNodesMeanNormalModelPartNonHistorical();

    // We compute the nodal thickness
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(geometry_number_of_elements); ++i) {
        auto it_elem = elements_array.begin() + i;

        // We get the thickness
        KRATOS_DEBUG_ERROR_IF_NOT(it_elem->GetProperties().Has(THICKNESS)) << "ERROR:: THICKNESS NOT DEFINED" << std::endl;
        const double thickness = it_elem->GetProperties()[THICKNESS];

        auto geom = it_elem->GetGeometry();
        for (IndexType i = 0; i < TNumNodes; ++i) {
            auto& node = geom[i];
            node.SetLock();
            node.GetValue(THICKNESS) += thickness;
            node.GetValue(NODAL_AREA) += 1.0;
            node.UnSetLock();
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(geometry_number_of_nodes); ++i) {
        auto it_node = nodes_array.begin() + i;
        double& thickness = it_node->GetValue(THICKNESS);
        thickness /= it_node->GetValue(NODAL_AREA);
    }

    // We create the new nodes
    const SizeType number_of_layers = mThisParameters["number_of_layers"].GetInt();
    for(IndexType i = 0; i < geometry_number_of_nodes; ++i) {
        auto it_node = nodes_array.begin() + i;

        const array_1d<double, 3>& normal = it_node->GetValue(NORMAL);
        const double thickness = it_node->GetValue(THICKNESS);
        array_1d<double, 3> coordinates = it_node->Coordinates() - 0.5 * normal * thickness;
        const double delta_thickness = thickness/static_cast<double>(number_of_layers);

        IndexType node_id = total_number_of_nodes + i + 1;
        NodeType::Pointer p_node0 = auxiliar_model_part.CreateNewNode(node_id, coordinates[0], coordinates[1], coordinates[2]);

        // Set the DOFs in the nodes
        for (auto it_dof = dofs.begin(); it_dof != dofs.end(); ++it_dof)
            p_node0->pAddDof(*it_dof);

        // We copy the step data
        CopyVariablesList(p_node0, p_node_begin);

        // We add the node to the external layers
        if (create_submodelparts_external_layers) auxiliar_model_part_lower.AddNode(p_node0);

        for (IndexType j = 0; j < number_of_layers; ++j) {
            coordinates += normal * delta_thickness;
            node_id = (j + 1) * geometry_number_of_nodes + total_number_of_nodes + i + 1;
            NodeType::Pointer p_node1 = auxiliar_model_part.CreateNewNode(node_id, coordinates[0], coordinates[1], coordinates[2]);

            // Set the DOFs in the nodes
            for (auto it_dof = dofs.begin(); it_dof != dofs.end(); ++it_dof)
                p_node1->pAddDof(*it_dof);

            // We copy the step data
            CopyVariablesList(p_node1, p_node_begin);

            // We add the node to the external layers
            if (create_submodelparts_external_layers && j == (number_of_layers - 1)) auxiliar_model_part_upper.AddNode(p_node1);
        }

        // We set the flag TO_ERASE for later remove the nodes
        it_node->Set(TO_ERASE, true);
    }

    const std::string& element_name = mThisParameters["element_name"].GetString();
    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    KRATOS_ERROR_IF_NOT(r_clone_element.GetGeometry().size() == 2 * TNumNodes) << "ERROR: Element " << element_name << " has a different number of nodes to " << 2 * TNumNodes << std::endl;

    // We will save the list of properties ids to later set the CL
    std::unordered_set<IndexType> set_id_properties;

    // We create the new elements
    IndexType element_counter = total_number_of_elements;
    for(IndexType i = 0; i < geometry_number_of_elements; ++i) {
        auto it_elem = elements_array.begin() + i;

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
            auxiliar_model_part.CreateNewElement(element_name, element_counter, element_node_ids, p_prop);
        }

        // We set the flag TO_ERASE for later remove the elements
        it_elem->Set(TO_ERASE, true);
    }

    // We create the conditions if necessary
    if (create_submodelparts_external_layers) {
        const std::string condition_name = "SurfaceCondition3D" + std::to_string(TNumNodes) + "N";
        IndexType condition_counter = total_number_of_conditions;
        for(IndexType i = 0; i < geometry_number_of_elements; ++i) {
            auto it_elem = elements_array.begin() + i;

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
            auxiliar_model_part_lower.CreateNewCondition(condition_name, condition_counter, lower_condition_node_ids, p_prop);
            condition_counter++;
            auxiliar_model_part_upper.CreateNewCondition(condition_name, condition_counter, upper_condition_node_ids, p_prop);
        }
    }

    // We reassign a new constitutive law
    const std::string& new_constitutive_law_name = mThisParameters["new_constitutive_law_name"].GetString();
    if (new_constitutive_law_name != "") {
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get(new_constitutive_law_name).Clone();
        for (auto id_prop : set_id_properties) {
            auto p_prop = geometry_model_part.pGetProperties(id_prop);
            p_prop->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
        }
    }

    // Finally we remove the old nodes and elements
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // We copy the new model part to the original one
    geometry_model_part.AddNodes( auxiliar_model_part.NodesBegin(), auxiliar_model_part.NodesEnd() );
    geometry_model_part.AddElements( auxiliar_model_part.ElementsBegin(), auxiliar_model_part.ElementsEnd() );
    
    KRATOS_ERROR << "something wrong with the creation of the modelparts...the author of this file should take a look" << std::endl;

    // We copy the external layers
/*    if (create_submodelparts_external_layers) {
        const std::string name_upper = "Upper_"+model_part_name;
        ModelPart* p_upper_model_part = append_submodelparts_external_layers ? 
                                &geometry_model_part.CreateSubModelPart(name_upper) : 
                                mrThisModelPart.CreateSubModelPart(name_upper);
        p_upper_model_part->AddNodes( auxiliar_model_part_upper.NodesBegin(), auxiliar_model_part_upper.NodesEnd() );
        p_upper_model_part->AddConditions( auxiliar_model_part_upper.ConditionsBegin(), auxiliar_model_part_upper.ConditionsEnd() );
        const std::string name_lower = "Lower_"+model_part_name;
        ModelPart* p_lower_model_part = append_submodelparts_external_layers ? 
            geometry_model_part.CreateSubModelPart(name_lower) : 
            mrThisModelPart.CreateSubModelPart(name_lower);
        p_lower_model_part->AddNodes( auxiliar_model_part_lower.NodesBegin(), auxiliar_model_part_lower.NodesEnd() );
        p_lower_model_part->AddConditions( auxiliar_model_part_lower.ConditionsBegin(), auxiliar_model_part_lower.ConditionsEnd() );
    }

    // We add to the computing model part if available
    const std::string& computing_model_part_name = mThisParameters["computing_model_part_name"].GetString();
    if (computing_model_part_name != "") {
        ModelPart& computing_model_part = mrThisModelPart.GetSubModelPart(computing_model_part_name);
        computing_model_part.AddNodes( auxiliar_model_part.NodesBegin(), auxiliar_model_part.NodesEnd() );
        computing_model_part.AddElements( auxiliar_model_part.ElementsBegin(), auxiliar_model_part.ElementsEnd() );
        computing_model_part.AddConditions( auxiliar_model_part_upper.ConditionsBegin(), auxiliar_model_part_upper.ConditionsEnd() );
        computing_model_part.AddConditions( auxiliar_model_part_lower.ConditionsBegin(), auxiliar_model_part_lower.ConditionsEnd() );
    }
*/
    // Reorder again all the IDs
    ReorderAllIds();

    // We initialize the new elements
    if (mThisParameters["initialize_elements"].GetBool()) {
        InitializeElements();
    }

    if (mThisParameters["export_to_mdpa"].GetBool()) {
        const std::string& output_name = mThisParameters["output_name"].GetString();
        std::ofstream output_file;
        ModelPartIO model_part_io(output_name, IO::WRITE);
        model_part_io.WriteModelPart(mrThisModelPart);
    }
    
    current_model.DeleteModelPart("AuxiliarModelPart");
    current_model.DeleteModelPart("upper");
    current_model.DeleteModelPart("lower");

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::ReorderAllIds(const bool ReorderAccordingShellConnectivity)
{
    KRATOS_ERROR << "author needs to review the creation of modelparts in this file" << std::endl;
    /*
    if (!ReorderAccordingShellConnectivity) {
        NodesArrayType& nodes_array = mrThisModelPart.Nodes();
        for(SizeType i = 0; i < nodes_array.size(); ++i)
            (nodes_array.begin() + i)->SetId(i + 1);
    } else {
        // The name of the submodelpart
        const std::string& model_part_name = mThisParameters["model_part_name"].GetString();
        ModelPart& geometry_model_part = model_part_name == "" ? mrThisModelPart : mrThisModelPart.GetSubModelPart(model_part_name);

        // Auxiliar model part where to store new nodes and elements
        //ModelPart auxiliar_model_part;
        
        KRATOS_ERROR << "author needs to review the creation of modelparts in this file" << std::endl;

        // Auxiliar values
        NodesArrayType& nodes_array = geometry_model_part.Nodes();
        const SizeType geometry_number_of_nodes = nodes_array.size();
        NodesArrayType& total_nodes_array = mrThisModelPart.Nodes();
        const SizeType total_number_of_nodes = total_nodes_array.size();

        // We reoder first all the nodes
        for(SizeType i = 0; i < total_number_of_nodes; ++i)
            (total_nodes_array.begin() + i)->SetId(total_number_of_nodes + i + 1);

        // We reoder now just the shell the nodes
        for(SizeType i = 0; i < geometry_number_of_nodes; ++i) {
            auto it_node = nodes_array.begin() + i;
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
        */
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::InitializeElements()
{
    ElementsArrayType& element_array = mrThisModelPart.Elements();
    for(SizeType i = 0; i < element_array.size(); ++i)
        (element_array.begin() + i)->Initialize();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
inline void ShellToSolidShellProcess<TNumNodes>::ComputeNodesMeanNormalModelPartNonHistorical()
{
    // Tolerance
    const double tolerance = std::numeric_limits<double>::epsilon();

    // The name of the submodelpart
    const std::string& model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& geometry_model_part = model_part_name == "" ? mrThisModelPart : mrThisModelPart.GetSubModelPart(model_part_name);

    // We iterate over the nodes
    NodesArrayType& nodes_array = geometry_model_part.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i)
        (nodes_array.begin() + i)->SetValue(NORMAL, ZeroVector(3));

    // Sum all the nodes normals
    ElementsArrayType& elements_array = geometry_model_part.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
        auto it_elem = elements_array.begin() + i;
        GeometryType& this_geometry = it_elem->GetGeometry();

        // Aux coordinates
        array_1d<double, 3> aux_coords;
        aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_geometry.Center());

        it_elem->SetValue(NORMAL, this_geometry.UnitNormal(aux_coords));

        const unsigned int number_nodes = this_geometry.PointsNumber();

        for (unsigned int i = 0; i < number_nodes; ++i) {
            auto& this_node = this_geometry[i];
            aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_node.Coordinates());
            const array_1d<double, 3>& normal = this_geometry.UnitNormal(aux_coords);
            auto& aux_normal = this_node.GetValue(NORMAL);
            for (unsigned int index = 0; index < 3; ++index) {
                #pragma omp atomic
                aux_normal[index] += normal[index];
            }
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) {
        auto it_node = nodes_array.begin() + i;

        array_1d<double, 3>& normal = it_node->GetValue(NORMAL);
        const double norm_normal = norm_2(normal);
        if (norm_normal > tolerance) normal /= norm_normal;
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

template class ShellToSolidShellProcess<3>;
template class ShellToSolidShellProcess<4>;
// class ShellToSolidShellProcess
} // namespace Kratos
