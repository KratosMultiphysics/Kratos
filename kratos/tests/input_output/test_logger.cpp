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
//                   
//


// Project includes
#include "testing/testing.h"
#include "input_output/logger.h"
#include "input_output/logger_output.h"


namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(LoggerMessageStream, KratosCoreFastSuite)
		{
			LoggerMessage message("label");

			message << "Test message with number " << 12 << 'e' << "00";

			KRATOS_CHECK_C_STRING_EQUAL(message.GetLabel().c_str(), "label");
			KRATOS_CHECK_C_STRING_EQUAL(message.GetMessage().c_str(), "Test message with number 12e00");
			KRATOS_CHECK_EQUAL(message.GetSeverity(), LoggerMessage::Severity::INFO);
            KRATOS_CHECK_EQUAL(message.GetCategory(), LoggerMessage::Category::STATUS);
            KRATOS_CHECK_EQUAL(message.GetLocation().GetFileName(), "Unknown");
            KRATOS_CHECK_EQUAL(message.GetLocation().GetFunctionName(), "Unknown");
            KRATOS_CHECK_EQUAL(message.GetLocation().GetLineNumber(), -1);

			message << LoggerMessage::Severity::DETAIL 
                << LoggerMessage::Category::CRITICAL 
                << KRATOS_CODE_LOCATION << std::endl;

			KRATOS_CHECK_C_STRING_EQUAL(message.GetMessage().c_str(), "Test message with number 12e00\n");
			KRATOS_CHECK_EQUAL(message.GetSeverity(), LoggerMessage::Severity::DETAIL);
            KRATOS_CHECK_EQUAL(message.GetCategory(), LoggerMessage::Category::CRITICAL);
            KRATOS_CHECK_NOT_EQUAL(message.GetLocation().GetFileName().find("test_logger.cpp"), std::string::npos);
            KRATOS_CHECK_EQUAL(message.GetLocation().GetFunctionName(), "virtual void Kratos::Testing::TestLoggerMessageStream::TestFunction()");
            KRATOS_CHECK_EQUAL(message.GetLocation().GetLineNumber(), 40);
		}

		KRATOS_TEST_CASE_IN_SUITE(LoggerOutput, KratosCoreFastSuite)
		{
			std::stringstream buffer;
			LoggerOutput output(buffer);

			LoggerMessage message("label");
			message << "Test message with number " << 12 << 'e' << "00";

			output.WriteMessage(message);
			KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "Test message with number 12e00");
		}

		KRATOS_TEST_CASE_IN_SUITE(LoggerStream, KratosCoreFastSuite)
		{
			std::stringstream buffer;
			LoggerOutput output(buffer);
			Logger::AddOutput(output);

			Logger("TestLabel") << "Test message with number " << 12 << 'e' << "00";

			KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "Test message with number 12e00");

			Logger("TestDetail") << Logger::Severity::DETAIL << "This log has detailed severity and will not be printed in output " 
				<< Logger::Category::CRITICAL << std::endl;

			// The message has DETAIL severity and should not be written
			KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "Test message with number 12e00");
		}

		KRATOS_TEST_CASE_IN_SUITE(CheckPoint, KratosCoreFastSuite)
		{
			std::stringstream buffer;
			LoggerOutput output(buffer);

			KRATOS_CHECK_POINT("TestCheckPoint") << "The value in check point is " << 3.14;

#if defined(KRATOS_ENABLE_CHECK_POINT)
			KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "The value in check point is 3.14");
#else
			KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), ""); // should print noting
#endif
		}

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfo, KratosCoreFastSuite)
        {
            std::stringstream buffer;
			LoggerOutput output(buffer);
			Logger::AddOutput(output);

			KRATOS_INFO("TestInfo") << "Test info message";

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "Test info message");
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoIf, KratosCoreFastSuite)
        {
            std::stringstream buffer;
			LoggerOutput output(buffer);
			Logger::AddOutput(output);

            KRATOS_INFO_IF("TestInfo", true) << "Test info message";
            KRATOS_INFO_IF("TestInfo", false) << "This should not appear";

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "Test info message");
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoFirst, KratosCoreFastSuite)
        {
            std::stringstream buffer;
			LoggerOutput output(buffer);
			Logger::AddOutput(output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_FIRST("TestInfo") << "Test info message - " << i;
            }

            std::cout << buffer.str().c_str() << std::endl;

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "Test info message - 0");
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoOnce, KratosCoreFastSuite)
        {
            std::stringstream buffer;
			LoggerOutput output(buffer);
			Logger::AddOutput(output);

			for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_ONCE("TestInfo", 4) << "Test info message - " << i;
            }

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "Test info message - 4");
        }

	}   // namespace Testing
}  // namespace Kratos.


