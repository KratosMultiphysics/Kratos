// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "test_utilities.h"
#include <fstream>

namespace Kratos::Testing
{

bool TestUtilities::CompareFiles(const std::filesystem::path& rPath1,
                                 const std::filesystem::path& rPath2)
{
    std::ifstream file_1(rPath1, std::ifstream::binary | std::ifstream::ate);
    std::ifstream file_2(rPath2, std::ifstream::binary | std::ifstream::ate);

    if (file_1.fail() || file_2.fail())
    {
        return false; // file problem
    }

    if (file_1.tellg() != file_2.tellg())
    {
        return false; // size mismatch
    }

    // seek back to beginning and use std::equal to compare contents
    file_1.seekg(0, std::ifstream::beg);
    file_2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(file_1.rdbuf()),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(file_2.rdbuf()));
}
} // namespace Kratos::Testing
