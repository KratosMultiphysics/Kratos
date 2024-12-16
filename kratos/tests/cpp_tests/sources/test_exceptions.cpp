//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Vicente Mataix Ferrandiz
//
//

// System includes


// External includes


// Project includes
#include "testing/testing.h"


namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ExceptionDefaultConstruction, KratosCoreFastSuite)
{
    try {
        throw Exception();
    }
    catch (Exception& e) {
        KRATOS_EXPECT_STREQ(e.what(), "Unknown Error\nin Unknown Location");
        KRATOS_EXPECT_EQ(e.where().CleanFileName(), "Unknown File");
        KRATOS_EXPECT_EQ(e.where().CleanFunctionName(), "Unknown Location");
        KRATOS_EXPECT_EQ(e.where().GetLineNumber(), 0);
    }
}

/** Checks if CPP works as exception
* Checks if CPP works as exception
*/
KRATOS_TEST_CASE_IN_SUITE(CPPExceptionTest, KratosCoreFastSuite)
{
    try {
        throw 1;
    } catch (int e) {
        KRATOS_EXPECT_EQ(e, 1);
    } catch (...)  {
        KRATOS_ERROR << std::endl;
    }
}

/** Checks if KRATOS_ERROR works as exception
 * Checks if KRATOS_ERROR works as exception
 */
KRATOS_TEST_CASE_IN_SUITE(KRATOS_ERRORTest, KratosCoreFastSuite)
{
    try {
        KRATOS_ERROR << std::endl;
    } catch(Kratos::Exception& e){
        KRATOS_EXPECT_TRUE(true);
    } catch (...)  {
        KRATOS_ERROR << "Default Exception"<< std::endl;
    }
}

/** Checks if KRATOS_ERROR_IF works as exception
 * Checks if KRATOS_ERROR_IF works as exception
 */
KRATOS_TEST_CASE_IN_SUITE(KRATOS_ERROR_IFTest, KratosCoreFastSuite)
{
    try {
        KRATOS_ERROR_IF(true) << std::endl;
    } catch(Kratos::Exception& e){
        KRATOS_EXPECT_TRUE(true);
    } catch (...)  {
        KRATOS_ERROR << "Default Exception"<< std::endl;
    }
}

/** Checks if KRATOS_ERROR_IF_NOT works as exception
 * Checks if KRATOS_ERROR_IF_NOT works as exception
 */
KRATOS_TEST_CASE_IN_SUITE(KRATOS_ERROR_IF_NOTTest, KratosCoreFastSuite)
{
    try {
        KRATOS_ERROR_IF_NOT(false) << std::endl;
    } catch(Kratos::Exception& e){
        KRATOS_EXPECT_TRUE(true);
    } catch (...)  {
        KRATOS_ERROR << "Default Exception"<< std::endl;
    }
}

}
}  // namespace Kratos.
