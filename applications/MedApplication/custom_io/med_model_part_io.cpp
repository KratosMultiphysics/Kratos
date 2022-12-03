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
#include "med.h"

// Project includes
#include "med_model_part_io.h"


namespace Kratos {

MedModelPartIO::MedModelPartIO(const std::filesystem::path& rFileName)
    : mFileName(rFileName)
{
    KRATOS_TRY

    mFileHandle = MEDfileOpen("test1.med", MED_ACC_RDWR);
    if (mFileHandle < 0) {
        KRATOS_WATCH("Erreur à la creation du fichier");
    }


    KRATOS_CATCH("")
}

MedModelPartIO::~MedModelPartIO()
{
    if (MEDfileClose(mFileHandle) < 0)
    std::cout << "Closing failed ...";
}

void MedModelPartIO::ReadModelPart(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR << "ReadModelPart is not yet implemented ..." << std::endl;

    KRATOS_CATCH("")
}

void MedModelPartIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR << "WriteModelPart is not yet implemented ..." << std::endl;

    KRATOS_CATCH("")
}

} // namespace Kratos.
