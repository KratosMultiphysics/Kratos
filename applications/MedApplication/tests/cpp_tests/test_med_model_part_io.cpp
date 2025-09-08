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

KRATOS_TEST_CASE_IN_SUITE(MedModelpartIO_NonExistingFile_read, KratosMedFastSuite)
{
    const std::filesystem::path file_path(std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".txt");
    KRATOS_EXPECT_FALSE(std::filesystem::exists(file_path)); // make sure there are no leftovers

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MedModelPartIO dummy(file_path),
        "File \""+file_path.string()+"\" does not exist!");

    KRATOS_EXPECT_FALSE(std::filesystem::exists(file_path));
}

KRATOS_TEST_CASE_IN_SUITE(MedModelpartIO_NonExistingFile_write, KratosMedFastSuite)
{
    const std::filesystem::path file_path(std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".txt");
    KRATOS_EXPECT_FALSE(std::filesystem::exists(file_path)); // make sure there are no leftovers

    MedModelPartIO(file_path, IO::WRITE);
    KRATOS_EXPECT_TRUE(std::filesystem::exists(file_path));

    std::filesystem::remove(file_path); // clean leftover
}

KRATOS_TEST_CASE_IN_SUITE(MedModelpartIO_TextFile, KratosMedFastSuite)
{
    const std::filesystem::path file_path(std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".txt");
    std::ofstream output(file_path); // create a dummy file (that is not a hdf file)

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MedModelPartIO dummy(file_path),
        "A problem with HDF occured while trying to open file \""+file_path.string()+"\"!");

    std::filesystem::remove(file_path); // clean leftover
}

} // Kratos::Testing
