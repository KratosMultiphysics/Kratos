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
#include "med_application.h"


namespace Kratos {

KratosMedApplication::KratosMedApplication():
    KratosApplication("MedApplication")
{
    // check if the library that was used to compile is the same as the one that is loaded at runtime
    med_int v_hdf_major, v_hdf_minor, v_hdf_release;
    med_int v_med_major, v_med_minor, v_med_release;

    MEDlibraryHdfNumVersion(&v_hdf_major, &v_hdf_minor, &v_hdf_release);
    MEDlibraryNumVersion(&v_med_major, &v_med_minor, &v_med_release);

    KRATOS_WATCH(v_hdf_major)
    KRATOS_WATCH(v_hdf_minor)
    KRATOS_WATCH(v_hdf_release)

    KRATOS_WATCH(v_med_major)
    KRATOS_WATCH(v_med_minor)
    KRATOS_WATCH(v_med_release)

    KRATOS_WATCH(MED_MAJOR_NUM)
    KRATOS_WATCH(MED_NUM_MINEUR)
    KRATOS_WATCH(MED_NUM_RELEASE)

    // TODO check macros vs what comes from the functions (both for HDF abd MED)
    // might indicate that versions are differen, aka different versions of the library loaded at runtime!
    // probably do in the MedApplication class, in register (actually also the filehandle type check could/should be done there ...)
}

void KratosMedApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosMedApplication..." << std::endl;
}

}  // namespace Kratos.
