// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes
// #include "med.h"

// Project includes
#include "med_model_part_io.h"


namespace Kratos {

MedModelPartIO::MedModelPartIO(const std::filesystem::path& rFileName)
    : mFileName(rFileName)
{

}

void MedModelPartIO::ReadModelPart(ModelPart& rThisModelPart)
{
    KRATOS_ERROR << "Calling base class method (ReadModelPart). Please check the definition of derived class" << std::endl;
}

} // namespace Kratos.
