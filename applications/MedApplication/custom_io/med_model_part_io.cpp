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
    KRATOS_TRY


    KRATOS_CATCH("")
}

void MedModelPartIO::ReadModelPart(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR << "Calling base class method (ReadModelPart). Please check the definition of derived class" << std::endl;

    KRATOS_CATCH("")
}

void MedModelPartIO::WriteModelPart(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    // writing should not modify the original object
    // hence using a const reference in the following
    const auto& r_model_part = rThisModelPart;

    KRATOS_ERROR << "Calling base class method (WriteModelPart). Please check the definition of derived class" << std::endl;

    KRATOS_CATCH("")
}

} // namespace Kratos.
