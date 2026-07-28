// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_modelers/mmg/mmg_modeler.h"
#include "custom_processes/mmg/mmg_process.h"

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

template<MMGLibrary TMMGLibrary>
MmgModeler<TMMGLibrary>::MmgModeler(
    Model& rModel,
    Parameters ModelerParameters
    ) : Modeler(rModel, ModelerParameters)
      , mpModel(&rModel)
{
    KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
        << "MmgModeler requires \"model_part_name\" in the parameters." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
Modeler::Pointer MmgModeler<TMMGLibrary>::Create(
    Model& rModel,
    const Parameters ModelParameters
    ) const
{
    return Kratos::make_shared<MmgModeler<TMMGLibrary>>(rModel, ModelParameters);
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgModeler<TMMGLibrary>::SetupModelPart()
{
    KRATOS_TRY;

    const std::string model_part_name = mParameters["model_part_name"].GetString();
    KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(model_part_name))
        << "MmgModeler: ModelPart \"" << model_part_name << "\" not found in the model." << std::endl;

    ModelPart& r_model_part = mpModel->GetModelPart(model_part_name);

    // Build the process parameters by forwarding all keys except model_part_name
    Parameters process_parameters(mParameters);
    if (process_parameters.Has("model_part_name"))
        process_parameters.RemoveValue("model_part_name");

    MmgProcess<TMMGLibrary> mmg_process(r_model_part, process_parameters);
    mmg_process.Execute();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
const Parameters MmgModeler<TMMGLibrary>::GetDefaultParameters() const
{
    // Minimal defaults — all MmgProcess parameters keep their own defaults
    // and are validated when the process is constructed in SetupModelPart.
    const Parameters default_parameters = Parameters(R"({
        "model_part_name" : ""
    })");
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class KRATOS_API(MESHING_APPLICATION) MmgModeler<MMGLibrary::MMG2D>;
template class KRATOS_API(MESHING_APPLICATION) MmgModeler<MMGLibrary::MMG3D>;
template class KRATOS_API(MESHING_APPLICATION) MmgModeler<MMGLibrary::MMGS>;

} // namespace Kratos
