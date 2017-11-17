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

			message << LoggerMessage::Severity::DETAIL 
				<< LoggerMessage::Category::CRITICAL << std::endl;

			KRATOS_CHECK_C_STRING_EQUAL(message.GetMessage().c_str(), "Test message with number 12e00\n");
			KRATOS_CHECK_EQUAL(message.GetSeverity(), LoggerMessage::Severity::DETAIL);
			KRATOS_CHECK_EQUAL(message.GetCategory(), LoggerMessage::Category::CRITICAL);
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

	}   // namespace Testing
}  // namespace Kratos.


