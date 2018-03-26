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

// External includes

// Project includes
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
        "element_name"    : "SolidShellElementSprism3D6N",
        "model_part_name" : "",
        "number_of_layers": 1,
        "export_to_mdpa"  : false,
        "output_name"     : "output"
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

    // The name of the submodelpart
    const std::string& model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& geometry_model_part = model_part_name == "" ? mrThisModelPart : mrThisModelPart.GetSubModelPart(model_part_name);

    // Auxiliar model part where to store new nodes and elements
    ModelPart auxiliar_model_part;

    // Auxiliar values
    NodesArrayType& nodes_array = geometry_model_part.Nodes();
    ElementsArrayType& elements_array = geometry_model_part.Elements();
    const SizeType geometry_number_of_nodes = nodes_array.size();
    const SizeType geometry_number_of_elements = elements_array.size();
    const SizeType total_number_of_nodes = mrThisModelPart.Nodes().size();
    const SizeType total_number_of_elements = mrThisModelPart.Elements().size();

    // First we reoder the ids
    ReorderAllIds();

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

        for (IndexType j = 0; j < number_of_layers; ++j) {
            coordinates += normal * delta_thickness;
            node_id = (j + 1) * geometry_number_of_nodes + total_number_of_nodes + i + 1;
            NodeType::Pointer p_node1 = auxiliar_model_part.CreateNewNode(node_id, coordinates[0], coordinates[1], coordinates[2]);
            // Set the DOFs in the nodes
            for (auto it_dof = dofs.begin(); it_dof != dofs.end(); ++it_dof)
                p_node1->pAddDof(*it_dof);
        }

        // We set the flag TO_ERASE for later remove the nodes
        it_node->Set(TO_ERASE, true);
    }

    const std::string& element_name = mThisParameters["element_name"].GetString();
    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    KRATOS_ERROR_IF_NOT(r_clone_element.GetGeometry().size() == 2 * TNumNodes) << "ERROR: Element " << element_name << " has a different number of nodes to " << 2 * TNumNodes << std::endl;

    // We create the new elements
    IndexType element_counter = total_number_of_elements;
    for(IndexType i = 0; i < geometry_number_of_elements; ++i) {
        auto it_elem = elements_array.begin() + i;

        auto p_prop = it_elem->pGetProperties();
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

    // Finally we remove the old nodes and elements
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // We copy the new model part to the original one
    mrThisModelPart.AddNodes( auxiliar_model_part.NodesBegin(), auxiliar_model_part.NodesEnd() );
    mrThisModelPart.AddElements( auxiliar_model_part.ElementsBegin(), auxiliar_model_part.ElementsEnd() );
    
    // Reorder again all the IDs
    ReorderAllIds();

//     // We initialize the new elements
//     InitializeElements();

    if (mThisParameters["export_to_mdpa"].GetBool()) {
        const std::string& output_name = mThisParameters["output_name"].GetString();
        std::ofstream output_file;
        ModelPartIO model_part_io(output_name, IO::WRITE);
        model_part_io.WriteModelPart(mrThisModelPart);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNumNodes>
void ShellToSolidShellProcess<TNumNodes>::ReorderAllIds()
{
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    for(SizeType i = 0; i < nodes_array.size(); ++i)
        (nodes_array.begin() + i)->SetId(i + 1);

    ElementsArrayType& element_array = mrThisModelPart.Elements();
    for(SizeType i = 0; i < element_array.size(); ++i)
        (element_array.begin() + i)->SetId(i + 1);
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

    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i)
        (nodes_array.begin() + i)->SetValue(NORMAL, ZeroVector(3));

    // Sum all the nodes normals
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

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

template class ShellToSolidShellProcess<3>;
template class ShellToSolidShellProcess<4>;
// class ShellToSolidShellProcess
} // namespace Kratos
