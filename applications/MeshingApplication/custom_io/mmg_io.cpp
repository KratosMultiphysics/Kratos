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

    // TODO: Read JSON of colors

    // Some information
    MMGMeshInfo<TMMGLibrary> mmg_mesh_info;
    mMmmgUtilities.PrintAndGetMmgMeshInfo(mmg_mesh_info);

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
    mMmmgUtilities.GenerateMeshDataFromModelPart(rModelPart, mColors, aux_ref_cond, aux_ref_elem);

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

    // TODO: Write the colors to a JSON

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class MmgIO<MMGLibrary::MMG2D>;
template class MmgIO<MMGLibrary::MMG3D>;
template class MmgIO<MMGLibrary::MMGS>;

}// namespace Kratos.
