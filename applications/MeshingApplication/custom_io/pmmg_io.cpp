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
#include "custom_io/pmmg_io.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"

// NOTE: The following contains the license of the ParMMG library
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

template<PMMGLibrary TPMMGLibrary>
ParMmgIO<TPMMGLibrary>::ParMmgIO(
        std::string const &rFilename,
        Parameters ThisParameters,
        const Flags Options
)
        : mFilename(rFilename), mThisParameters(ThisParameters), mOptions(Options)
{
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Check the mode
    if (mOptions.Is(IO::APPEND)) {
        KRATOS_ERROR << "APPEND not compatible with ParMmgIO" << std::endl;
    }

    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::SetOuputFile(rFilename + ".time");

    /* We restart the MMG mesh and solution */
    mParMmmgUtilities.SetEchoLevel(mThisParameters["echo_level"].GetInt());
    mParMmmgUtilities.InitMesh();
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgIO<TPMMGLibrary>::ReadModelPart(ModelPart &rModelPart)
{
    KRATOS_TRY ;

        // Automatically read the mesh
        mParMmmgUtilities.InputMesh(mFilename);

        // Automatically read the solution
        mParMmmgUtilities.InputSol(mFilename);

        // Read JSON of colors
        std::unordered_map<IndexType, std::vector<std::string>> colors;  /// Where the sub model parts IDs are stored
        AssignUniqueModelPartCollectionTagUtility::ReadTagsFromJson(mFilename, colors);

        // Create the submodelparts
        const std::string &r_main_name = rModelPart.Name();
        for (auto &r_color : colors) {
            for (auto &r_name : r_color.second) {
                if (!rModelPart.HasSubModelPart(r_name) && r_main_name != r_name) {
                    rModelPart.CreateSubModelPart(r_name);
                }
            }
        }

        // Some information
        PMMGMeshInfo<TPMMGLibrary> mmg_mesh_info;
        mParMmmgUtilities.PrintAndGetParMmgMeshInfo(mmg_mesh_info);

        // Creating auxiliar maps of pointers
        std::unordered_map<IndexType, Condition::Pointer> ref_condition; /// Reference condition
        std::unordered_map<IndexType, Element::Pointer> ref_element;     /// Reference element

        // Fill the maps
        mParMmmgUtilities.WriteReferenceEntitities(rModelPart, mFilename, ref_condition, ref_element);

        // Writing the new mesh data on the model part
        NodeType::DofsContainerType empty_dofs;
        mParMmmgUtilities.WriteMeshDataToModelPart(rModelPart, colors, empty_dofs, mmg_mesh_info, ref_condition,
                                                   ref_element);

        // Writing the new solution data on the model part
        mParMmmgUtilities.WriteSolDataToModelPart(rModelPart);

        /* After that we reorder nodes, conditions and elements: */
        mParMmmgUtilities.ReorderAllIds(rModelPart);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgIO<TPMMGLibrary>::WriteModelPart(ModelPart &rModelPart)
{
    KRATOS_TRY ;

        // The auxiliar color maps
        ColorsMapType aux_ref_cond, aux_ref_elem;

        // We initialize the mesh data with the given modelpart
        std::unordered_map<IndexType, std::vector<std::string>> colors;  /// Where the sub model parts IDs are stored
        mParMmmgUtilities.GenerateMeshDataFromModelPart(rModelPart, colors, aux_ref_cond, aux_ref_elem);

        // Generate the maps of reference
        std::unordered_map<IndexType, Element::Pointer> ref_element;   /// Reference element
        std::unordered_map<IndexType, Condition::Pointer> ref_condition; /// Reference condition
        mParMmmgUtilities.GenerateReferenceMaps(rModelPart, aux_ref_cond, aux_ref_elem, ref_condition, ref_element);

        // We initialize the solution data with the given modelpart
        mParMmmgUtilities.GenerateSolDataFromModelPart(rModelPart);

        // Check if the number of given entities match with mesh size
//            mParMmmgUtilities.CheckMeshData(); //not implemented

        // Automatically save the mesh
        mParMmmgUtilities.OutputMesh(mFilename);

        // Automatically save the solution
        mParMmmgUtilities.OutputSol(mFilename);

        // Output the reference files
        mParMmmgUtilities.OutputReferenceEntitities(mFilename, ref_condition, ref_element);

        // Writing the colors to a JSON
        AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson(mFilename, colors);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/


template
class ParMmgIO<PMMGLibrary::PMMG3D>;

}// namespace Kratos.
