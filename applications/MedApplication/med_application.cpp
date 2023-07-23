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

// Project includes
#include "med_inc.h"
#include "med_application.h"


namespace Kratos {

KratosMedApplication::KratosMedApplication():
    KratosApplication("MedApplication")
{
    // check if the library that was used to compile is the same as the one that is loaded at runtime
    // they must match, otherwise random errors can occur
    med_int v_med_major, v_med_minor, v_med_release;
    MEDlibraryNumVersion(&v_med_major, &v_med_minor, &v_med_release);

    KRATOS_ERROR_IF(v_med_major   != MED_MAJOR_NUM ||
                    v_med_minor   != MED_NUM_MINEUR ||
                    v_med_release != MED_NUM_RELEASE) << "The MED library that was used during compilation (v" << v_med_major << "." << v_med_minor << "." << v_med_release << ") is different from the one loaded at runtime (v" << MED_MAJOR_NUM << "." << MED_NUM_MINEUR << "." << MED_NUM_RELEASE << ")!\nThis causes problems with reading/writing MED files, please check your paths for loading the library (e.g. \"PATH\" or \"LD_LIBRARY_PATH\")" << std::endl;
}

void KratosMedApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosMedApplication..." << std::endl;
}

}  // namespace Kratos.
