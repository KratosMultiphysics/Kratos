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
    Kratos::shared_ptr<std::fstream> pFile = Kratos::make_shared<std::fstream>();
    std::fstream::openmode OpenMode;

    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Set the mode
    if (mOptions.Is(IO::READ)) {
        OpenMode = std::fstream::in;
    } else if (mOptions.Is(IO::APPEND)) {
        OpenMode = std::fstream::in | std::fstream::app;
    } else if (mOptions.Is(IO::WRITE)) {
        OpenMode = std::fstream::out;
    } else {
        // If none of the READ, WRITE or APPEND are defined we will take READ as
        // default.
        OpenMode = std::fstream::in;
    }

    pFile->open(mFilename.c_str(), OpenMode);

    KRATOS_ERROR_IF_NOT(pFile->is_open()) << "Error opening mdpa file : " << mFilename.c_str() << std::endl;

    // Store the pointer as a regular std::iostream
    mpStream = pFile;

    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::SetOuputFile(rFilename + ".time");

    /* We restart the MMG mesh and solution */
    mEchoLevel = mThisParameters["echo_level"].GetInt();
    mMmmgUtilities.SetEchoLevel(mEchoLevel);
    mMmmgUtilities.InitMesh();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
MmgIO<TMMGLibrary>::MmgIO(
    Kratos::shared_ptr<std::iostream> Stream,
    Parameters ThisParameters,
    const Flags Options
    )
    : mThisParameters(ThisParameters)
     ,mOptions(Options)
{
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Check if the pointer is valid
    KRATOS_ERROR_IF(Stream == nullptr) << "Error: ModelPartIO Stream is invalid " << std::endl;

    // Check if the pointer was .reset() or never initialized and if its a NULL pointer)
    KRATOS_ERROR_IF(Stream == nullptr || Stream == Kratos::shared_ptr<std::iostream>(NULL)) << "Error: ModelPartIO Stream is invalid " << std::endl;

    mpStream = Stream;

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



    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgIO<TMMGLibrary>::WriteModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;



    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class MmgIO<MMGLibrary::MMG2D>;
template class MmgIO<MMGLibrary::MMG3D>;
template class MmgIO<MMGLibrary::MMGS>;

}// namespace Kratos.
