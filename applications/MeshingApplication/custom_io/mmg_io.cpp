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
    mEchoLevel = mThisParameters["echo_level"].GetInt();
    mMmmgUtilities.SetEchoLevel(mEchoLevel);
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
    std::ifstream infile(mFilename + ".json");
    KRATOS_ERROR_IF_NOT(infile.good()) << "Materials file: " << mFilename  + ".json" << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters json_text(buffer.str());
    for (auto it_param = json_text.begin(); it_param != json_text.end(); ++it_param) {
        const std::vector<std::string>& r_sub_model_part_names = it_param->GetStringArray();
        colors.insert(std::pair<IndexType,std::vector<std::string>>({std::stoi(it_param.name()), r_sub_model_part_names}));
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
    for (auto& r_color : colors)
        ref_condition[r_color.first] = r_clone_condition.Create(0, r_clone_condition.GetGeometry(), p_auxiliar_prop);

    /* Elements */
    const std::string element_type_name = (Dimension == 2) ? "Element2D3N" : (TMMGLibrary == MMGLibrary::MMG3D) ? "Element3D4N" : "Element3D3N";
    Element const& r_clone_element = KratosComponents<Element>::Get(element_type_name);
    ref_element[0] = r_clone_element.Create(0, r_clone_element.GetGeometry(), p_auxiliar_prop);
    for (auto& r_color : colors)
        ref_element[r_color.first] = r_clone_element.Create(0, r_clone_element.GetGeometry(), p_auxiliar_prop);


    // Writing the new mesh data on the model part
    NodeType::DofsContainerType empty_dofs;
    mMmmgUtilities.WriteMeshDataToModelPart(rModelPart, colors, empty_dofs, mmg_mesh_info, ref_condition, ref_element);

    // Writing the new solution data on the model part
    mMmmgUtilities.WriteSolDataToModelPart(rModelPart);

    /* After that we reorder nodes, conditions and elements: */
    mMmmgUtilities.ReorderAllIds(rModelPart);

    /* We assign flags and clear the auxiliar model parts created to reassing the flags */
    mMmmgUtilities.AssignAndClearAuxiliarSubModelPartForFlags(rModelPart);

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

    // We initialize the solution data with the given modelpart
    mMmmgUtilities.GenerateSolDataFromModelPart(rModelPart);

    // Check if the number of given entities match with mesh size
    mMmmgUtilities.CheckMeshData();

    // Retrieving current step
    const IndexType step = rModelPart.GetProcessInfo()[STEP];

    // Automatically save the mesh
    mMmmgUtilities.OutputMesh(mFilename, false, step);

    // Automatically save the solution
    mMmmgUtilities.OutputSol(mFilename, false, step);

    // Writing the colors to a JSON
    Parameters color_json;
    for (auto& r_color : colors) {
        Parameters names_array;
        for (auto& r_model_part_name : r_color.second) {
            names_array.Append(r_model_part_name);
        }
        color_json.AddValue(std::to_string(r_color.first), names_array);
    }

    const std::string& r_json_text = color_json.PrettyPrintJsonString();

    std::filebuf buffer;
    buffer.open(mFilename + ".json",std::ios::out);
    std::ostream os(&buffer);
    os << r_json_text;
    buffer.close();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class MmgIO<MMGLibrary::MMG2D>;
template class MmgIO<MMGLibrary::MMG3D>;
template class MmgIO<MMGLibrary::MMGS>;

}// namespace Kratos.
