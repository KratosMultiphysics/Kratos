//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/exception.h"

namespace Kratos 
{
    namespace Testing 
    {
        /// Tests
       
        /** Checks if CPP works as exception
         * Checks if CPP works as exception
         */
        
        KRATOS_TEST_CASE_IN_SUITE(CPPExceptionTest, KratosCoreFastSuite)
        {
            try {
                throw 1;
            } catch (int e) {
                KRATOS_CHECK_EQUAL(e, 1);
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
                KRATOS_CHECK(true);
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
                KRATOS_CHECK(true);
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
                KRATOS_CHECK(true);
            } catch (...)  { 
                KRATOS_ERROR << "Default Exception"<< std::endl; 
            } 
        }
        
    } // namespace Testing
}  // namespace Kratos.

