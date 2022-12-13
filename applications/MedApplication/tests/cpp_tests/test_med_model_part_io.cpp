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
#include "testing/testing.h"
#include "custom_io/med_model_part_io.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(MedModelpartIO_NonExistingFile, KratosMedFastSuite)
{
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        MedModelPartIO("random_non_existing_file"),
        "File \"random_non_existing_file\" does not exist!");
}

} // Kratos::Testing
