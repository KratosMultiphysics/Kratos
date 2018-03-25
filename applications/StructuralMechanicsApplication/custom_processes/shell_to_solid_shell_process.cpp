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
#include "utilities/mortar_utilities.h"
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
        "number_of_layers": 1
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
    ModelPart& geometry_model_part = mrThisModelPart.GetSubModelPart(model_part_name);

    // First we reoder the ids
    ReorderAllIds();

    // We copy the dof from the first node
    NodesArrayType& nodes_array = geometry_model_part.Nodes();
    NodeType::DofsContainerType dofs = nodes_array.begin()->GetDofs();

    // We iterate over the elements
    ElementsArrayType& elements_array = geometry_model_part.Elements();

    // We initialize the thickness
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->SetValue(THICKNESS, 0.0);
        it_node->SetValue(NODAL_AREA, 0.0);
    }

    // We initialize the normal
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i)
        (elements_array.begin() + i)->SetValue(NORMAL, ZeroVector(3));

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(geometry_model_part);

    // Auxiliar values
    const SizeType number_of_nodes = mrThisModelPart.Nodes().size();
    const SizeType number_of_elements = mrThisModelPart.Elements().size();

    // We compute the nodal thickness
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
        auto it_elem = elements_array.begin() + i;

        // We get the thickness
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
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        double& thickness = it_node->GetValue(THICKNESS);
        thickness /= it_node->GetValue(NODAL_AREA);
    }

    // We create the new nodes
    IndexType node_counter = number_of_nodes;
    const SizeType number_of_layers = mThisParameters["number_of_layers"].GetInt();
    for(IndexType i = 0; i < nodes_array.size(); ++i) {
        auto it_node = nodes_array.begin() + i;

        const array_1d<double, 3>& normal = it_node->FastGetSolutionStepValue(NORMAL);
        const double thickness = it_node->GetValue(THICKNESS);
        array_1d<double, 3> coordinates = it_node->Coordinates() - 0.5 * normal * thickness;
        const double delta_thickness = thickness/static_cast<double>(number_of_layers);

        mrThisModelPart.CreateNewNode(node_counter, coordinates[0], coordinates[1], coordinates[2]);
        node_counter++;

        for (IndexType j = 0; j < number_of_layers; ++j) {
            coordinates += normal * delta_thickness;
            mrThisModelPart.CreateNewNode(node_counter, coordinates[0], coordinates[1], coordinates[2]);
            node_counter++;
        }

        // We set the flag TO_ERASE for later remove the nodes
        it_node->Set(TO_ERASE, true);
    }

    const std::string& element_name = mThisParameters["element_name"].GetString();
    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    KRATOS_ERROR_IF_NOT(r_clone_element.GetGeometry().size() == TNumNodes) << "ERROR: Element " << element_name << " has a different number of nodes to " << TNumNodes << std::endl;

    // We create the new elements
    IndexType element_counter = number_of_elements;
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
        auto it_elem = elements_array.begin() + i;

        auto p_prop = it_elem->pGetProperties();
        for (IndexType j = 0; j < number_of_layers; ++j) {
            std::vector<IndexType> element_node_ids (2 * TNumNodes);
            for (IndexType k = 0; k < TNumNodes; ++k) {
                element_node_ids[k] = i * number_of_nodes + k;
                element_node_ids[k + TNumNodes] = (i + 1) * number_of_nodes + k;
            }
            mrThisModelPart.CreateNewElement(element_name, element_counter, element_node_ids,
            p_prop);
            element_counter++;
        }

        // We set the flag TO_ERASE for later remove the elements
        it_elem->Set(TO_ERASE, true);
    }

    // We initialize the new elements
    InitializeElements();

    // Reorder again all the IDs
    ReorderAllIds();

    // Finally we remove the old nodes and elements
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);
    
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

template class ShellToSolidShellProcess<3>;
template class ShellToSolidShellProcess<4>;
// class ShellToSolidShellProcess
} // namespace Kratos
