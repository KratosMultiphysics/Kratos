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
    // logging information about the used library for debugging purposes
    med_int v_med_major, v_med_minor, v_med_release;
    MEDlibraryNumVersion(&v_med_major, &v_med_minor, &v_med_release);

    med_int v_hdf_major, v_hdf_minor, v_hdf_release;
    MEDlibraryHdfNumVersion(&v_hdf_major, &v_hdf_minor, &v_hdf_release);

    // Note: the detail severity must be explicitly enabled, it is not shown by default
    KRATOS_DETAIL("MedApplication") << "Version of MED-library used during compilation: " << MED_MAJOR_NUM << "." << MED_NUM_MINEUR << "." << MED_NUM_RELEASE << std::endl;
    KRATOS_DETAIL("MedApplication") << "Version of MED-library loaded at runtime: " << v_med_major << "." << v_med_minor << "." << v_med_release << std::endl;
    KRATOS_DETAIL("MedApplication") << "Version of HDF-library used during compilation: " << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << std::endl;
    KRATOS_DETAIL("MedApplication") << "Version of HDF-library loaded at runtime: " << v_hdf_major << "." << v_hdf_minor << "." << v_hdf_release << std::endl;

    // check if the library that was used to compile is the same as the one that is loaded at runtime
    // they must match, otherwise random errors can occur
    KRATOS_ERROR_IF(v_med_major   != MED_MAJOR_NUM ||
                    v_med_minor   != MED_NUM_MINEUR ||
                    v_med_release != MED_NUM_RELEASE) << "The MED library that was used during compilation (v" << MED_MAJOR_NUM << "." << MED_NUM_MINEUR << "." << MED_NUM_RELEASE << ") is different from the one loaded at runtime (v" << v_med_major << "." << v_med_minor << "." << v_med_release << ")!\nThis causes problems with reading/writing MED files, please check your paths for loading the library (e.g. \"PATH\" or \"LD_LIBRARY_PATH\")" << std::endl;
}

void KratosMedApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosMedApplication..." << std::endl;
}

}  // namespace Kratos.
