// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "utilities/timer.h"
#include "custom_io/mmg_io.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "utilities/compare_elements_and_conditions_utility.h"

// NOTE: The following contains the license of the MMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

namespace Kratos
{

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
MmgIO<TMMGLibrary>::MmgIO(
    std::string const& rFilename,
    Parameters ThisParameters,
    const Flags Options
    )
    : mFilename(rFilename)
    , mThisParameters(ThisParameters)
    , mOptions(Options)
{
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Check the mode
    if (mOptions.Is(IO::APPEND)) {
        KRATOS_ERROR << "APPEND not compatible with MmgIO" << std::endl;
    }

    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::SetOuputFile(rFilename + ".time");

    /* We restart the MMG mesh and solution */
    mMmmgUtilities.SetEchoLevel(mThisParameters["echo_level"].GetInt());
    mMmmgUtilities.InitMesh();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgIO<TMMGLibrary>::ReadModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    // Automatically read the mesh
    mMmmgUtilities.InputMesh(mFilename);

    // Automatically read the solution
    mMmmgUtilities.InputSol(mFilename);

    // Read JSON of colors
    std::unordered_map<IndexType,std::vector<std::string>> colors;  /// Where the sub model parts IDs are stored
    AssignUniqueModelPartCollectionTagUtility::ReadTagsFromJson(mFilename, colors);

    // Create the submodelparts
    const std::string& r_main_name = rModelPart.Name();
    for (auto& r_color : colors) {
        for (auto& r_name : r_color.second) {
            if (!rModelPart.HasSubModelPart(r_name) && r_main_name != r_name) {
                rModelPart.CreateSubModelPart(r_name);
            }
        }
    }

    // Some information
    MMGMeshInfo<TMMGLibrary> mmg_mesh_info;
    mMmmgUtilities.PrintAndGetMmgMeshInfo(mmg_mesh_info);

    // Creating auxiliar maps of pointers
    std::unordered_map<IndexType,Condition::Pointer> ref_condition; /// Reference condition
    std::unordered_map<IndexType,Element::Pointer> ref_element;     /// Reference element

    // Getting auxiliar properties
    auto p_auxiliar_prop = rModelPart.CreateNewProperties(0);

    // Fill the maps
    /* Conditions */
    const std::string condition_type_name = (Dimension == 2) ? "Condition2D2N" : (TMMGLibrary == MMGLibrary::MMG3D) ? "SurfaceCondition3D3N" : "Condition3D2N";
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_type_name);
    ref_condition[0] = r_clone_condition.Create(0, r_clone_condition.GetGeometry(), p_auxiliar_prop);
    for (auto& r_color : colors) {
        ref_condition[r_color.first] = r_clone_condition.Create(0, r_clone_condition.GetGeometry(), p_auxiliar_prop);
    }

    /* Elements */
    const std::string element_type_name = (Dimension == 2) ? "Element2D3N" : (TMMGLibrary == MMGLibrary::MMG3D) ? "Element3D4N" : "Element3D3N";
    Element const& r_clone_element = KratosComponents<Element>::Get(element_type_name);
    ref_element[0] = r_clone_element.Create(0, r_clone_element.GetGeometry(), p_auxiliar_prop);
    for (auto& r_color : colors) {
        ref_element[r_color.first] = r_clone_element.Create(0, r_clone_element.GetGeometry(), p_auxiliar_prop);
    }

    // Writing the new mesh data on the model part
    NodeType::DofsContainerType empty_dofs;
    mMmmgUtilities.WriteMeshDataToModelPart(rModelPart, colors, empty_dofs, mmg_mesh_info, ref_condition, ref_element);

    // Writing the new solution data on the model part
    mMmmgUtilities.WriteSolDataToModelPart(rModelPart);

    /* After that we reorder nodes, conditions and elements: */
    mMmmgUtilities.ReorderAllIds(rModelPart);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgIO<TMMGLibrary>::WriteModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    // The auxiliar color maps
    ColorsMapType aux_ref_cond, aux_ref_elem;

    // We initialize the mesh data with the given modelpart
    std::unordered_map<IndexType,std::vector<std::string>> colors;  /// Where the sub model parts IDs are stored
    mMmmgUtilities.GenerateMeshDataFromModelPart(rModelPart, colors, aux_ref_cond, aux_ref_elem);

    // Generate the maps of reference
    std::unordered_map<IndexType,Element::Pointer>   ref_element;   /// Reference element
    std::unordered_map<IndexType,Condition::Pointer> ref_condition; /// Reference condition
    mMmmgUtilities.GenerateReferenceMaps(rModelPart, aux_ref_cond, aux_ref_elem, ref_condition, ref_element);

    /* ELEMENTS */
    std::string element_name;
    Parameters elem_reference_json;
    for (auto& r_elem : ref_element) {
        CompareElementsAndConditionsUtility::GetRegisteredName(*(r_elem.second), element_name);
        const std::string name = std::to_string(r_elem.first);
        elem_reference_json.AddEmptyValue(name);
        elem_reference_json[name].SetString(element_name);
    }

    const std::string& r_elem_json_text = elem_reference_json.PrettyPrintJsonString();

    std::filebuf elem_buffer;
    elem_buffer.open(mFilename + ".elem.ref.json",std::ios::out);
    std::ostream elem_os(&elem_buffer);
    elem_os << r_elem_json_text;
    elem_buffer.close();

    /* CONDITIONS */
    std::string condition_name;
    Parameters cond_reference_json;
    for (auto& r_cond : ref_condition) {
        CompareElementsAndConditionsUtility::GetRegisteredName(*(r_cond.second), condition_name);
        const std::string name = std::to_string(r_cond.first);
        cond_reference_json.AddEmptyValue(name);
        cond_reference_json[name].SetString(condition_name);
    }

    const std::string& r_cond_json_text = cond_reference_json.PrettyPrintJsonString();

    std::filebuf cond_buffer;
    cond_buffer.open(mFilename + ".cond.ref.json",std::ios::out);
    std::ostream cond_os(&cond_buffer);
    cond_os << r_cond_json_text;
    cond_buffer.close();

    // We initialize the solution data with the given modelpart
    mMmmgUtilities.GenerateSolDataFromModelPart(rModelPart);

    // Check if the number of given entities match with mesh size
    mMmmgUtilities.CheckMeshData();

    // Automatically save the mesh
    mMmmgUtilities.OutputMesh(mFilename);

    // Automatically save the solution
    mMmmgUtilities.OutputSol(mFilename);

    // Writing the colors to a JSON
    AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson(mFilename, colors);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class MmgIO<MMGLibrary::MMG2D>;
template class MmgIO<MMGLibrary::MMG3D>;
template class MmgIO<MMGLibrary::MMGS>;

}// namespace Kratos.
