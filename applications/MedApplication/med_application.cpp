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
#include "hdf5.h"

// Project includes
#include "med_application.h"


namespace Kratos {

KratosMedApplication::KratosMedApplication():
    KratosApplication("MedApplication")
{
    unsigned v_hdf_major, v_hdf_minor, v_hdf_release;
    H5get_libversion(&v_hdf_major, &v_hdf_minor, &v_hdf_release);

    // Note: the detail severity must be explicitly enabled, it is not shown by default
    KRATOS_DETAIL("MedApplication") << "Version of HDF-library used during compilation: " << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << std::endl;
    KRATOS_DETAIL("MedApplication") << "Version of HDF-library loaded at runtime: " << v_hdf_major << "." << v_hdf_minor << "." << v_hdf_release << std::endl;

    // check if the library that was used to compile is the same as the one that is loaded at runtime
    // they must match, otherwise random errors can occur
    KRATOS_ERROR_IF(v_hdf_major   != static_cast<unsigned>(H5_VERS_MAJOR) ||
                    v_hdf_minor   != static_cast<unsigned>(H5_VERS_MINOR)) << "The HDF5 library that was used during compilation (v" << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << ") is different from the one loaded at runtime (v" << v_hdf_major << "." << v_hdf_minor << "." << v_hdf_release << ")!\nThis causes problems with reading/writing MED files, please check your paths for loading the library (e.g. \"PATH\" or \"LD_LIBRARY_PATH\")" << std::endl;
}

void KratosMedApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosMedApplication..." << std::endl;
}

}  // namespace Kratos.
