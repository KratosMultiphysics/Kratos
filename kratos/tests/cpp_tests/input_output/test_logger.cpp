//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// Project includes
#include "testing/testing.h"
#include "input_output/logger.h"
#include "input_output/logger_table_output.h"
#include "includes/data_communicator.h"


namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(LoggerMessageStream, KratosCoreFastSuite)
        {
            LoggerMessage message("label");

            message << "Test message with number " << 12 << 'e' << "00";

            KRATOS_CHECK_C_STRING_EQUAL(message.GetLabel().c_str(), "label");
            if (DataCommunicator::GetDefault().Rank() == 0) KRATOS_CHECK_C_STRING_EQUAL(message.GetMessage().c_str(), "Test message with number 12e00");
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
            KRATOS_CHECK_EQUAL(message.GetLocation().GetFunctionName(), KRATOS_CURRENT_FUNCTION);
            KRATOS_CHECK_EQUAL(message.GetLocation().GetLineNumber(), 40);
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerOutput, KratosCoreFastSuite)
        {
            std::stringstream buffer;
            LoggerOutput output(buffer);

            LoggerMessage message("label");
            message << "Test message with number " << 12 << 'e' << "00";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "label: Test message with number 12e00" : "";

            output.WriteMessage(message);
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStream, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            Logger("TestLabel") << "Test message with number " << 12 << 'e' << "00";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestLabel: Test message with number 12e00" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());

            Logger("TestDetail") << Logger::Severity::DETAIL << "This log has detailed severity and will not be printed in output "
                << Logger::Category::CRITICAL << std::endl;

            // The message has DETAIL severity and should not be written (check that nothing was added to the buffer)
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(CheckPoint, KratosCoreFastSuite)
        {
            std::stringstream buffer;
            LoggerOutput output(buffer);

            KRATOS_CHECK_POINT("TestCheckPoint") << "The value in check point is " << 3.14;

#if defined(KRATOS_ENABLE_CHECK_POINT)
            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestCheckPoint: The value in check point is 3.14" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
#else
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), ""); // should print noting
#endif
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfo, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            KRATOS_INFO("TestInfo") << "Test info message";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestInfo: Test info message" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoIf, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            KRATOS_INFO_IF("TestInfo", true) << "Test info message";
            KRATOS_INFO_IF("TestInfo", false) << "This should not appear";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestInfo: Test info message" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoOnce, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_ONCE("TestInfo") << "Test info message - " << i;
            }

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestInfo: Test info message - 0" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoFirst, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_FIRST_N("TestInfo", 4) << ".";
            }

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestInfo: .TestInfo: .TestInfo: .TestInfo: ." : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamWarning, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            KRATOS_WARNING("TestWarning") << "Test warning message";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "[WARNING] TestWarning: Test warning message" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamWarningIf, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            KRATOS_WARNING_IF("TestWarning", true) << "Test warning message";
            KRATOS_WARNING_IF("TestWarning", false) << "This should not appear";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "[WARNING] TestWarning: Test warning message" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamWarningOnce, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_WARNING_ONCE("TestWarning") << "Test warning message - " << i;
            }

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "[WARNING] TestWarning: Test warning message - 0" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamWarningFirst, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_WARNING_FIRST_N("TestWarning", 4) << ".";
            }

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "[WARNING] TestWarning: .[WARNING] TestWarning: .[WARNING] TestWarning: .[WARNING] TestWarning: ." : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamDetail, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            p_output->SetSeverity(LoggerMessage::Severity::DETAIL);
            Logger::AddOutput(p_output);

            KRATOS_DETAIL("TestDetail") << "Test detail message";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestDetail: Test detail message" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamDetailIf, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            p_output->SetSeverity(LoggerMessage::Severity::DETAIL);
            Logger::AddOutput(p_output);

            KRATOS_DETAIL_IF("TestDetail", true) << "Test detail message";
            KRATOS_DETAIL_IF("TestDetail", false) << "This should not appear";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestDetail: Test detail message" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamDetailOnce, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            p_output->SetSeverity(LoggerMessage::Severity::DETAIL);
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_DETAIL_ONCE("TestDetail") << "Test detail message - " << i;
            }

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestDetail: Test detail message - 0" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamDetailFirst, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            p_output->SetSeverity(LoggerMessage::Severity::DETAIL);
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_DETAIL_FIRST_N("TestDetail", 4) << ".";
            }

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestDetail: .TestDetail: .TestDetail: .TestDetail: ." : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerTableOutput, KratosCoreFastSuite)
        {
            int rank = DataCommunicator::GetDefault().Rank();

            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerTableOutput(buffer, {"Time Step    ", "Iteration Number        ", "Convergence        ", "Is converged"}));
            Logger::AddOutput(p_output);

            std::stringstream reference_output;
            if (rank == 0)
                reference_output << "Time Step     Iteration Number         Convergence         Is converged " << std::endl;

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), reference_output.str().c_str());

            std::size_t time_step = 1;
            Logger("Time Step") << time_step << std::endl;
            if (rank == 0) reference_output << "1             ";

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), reference_output.str().c_str());

            Logger("Label") << "This log has a label which is not in the output columns and will not be printed in output " << std::endl;

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), reference_output.str().c_str());

            double convergence = 0.00;
            for(time_step = 2;time_step < 4; time_step++){
                Logger("Time Step") << time_step << std::endl;
                for(int iteration_number = 1 ; iteration_number < 3; iteration_number++){
                    convergence = 0.3 / (iteration_number * time_step);
                    Logger("Iteration Number") << iteration_number << std::endl;
                    Logger("Convergence") << convergence << std::endl;
                }
                if(convergence < 0.06) {
                    Logger("Is converged") << "Yes" << std::endl;
                }
                else {
                    Logger("Is converged") << "No" << std::endl;
                }

            }

            if (rank == 0) {
                reference_output << std::endl << "2             1                        0.15                ";
                reference_output << std::endl << "              2                        0.075               No           ";
                reference_output << std::endl << "3             1                        0.1                 ";
                reference_output << std::endl << "              2                        0.05                Yes          ";
            }

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), reference_output.str().c_str());

            std::cout << std::endl;
            std::cout << buffer.str() << std::endl;

        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerTableDistributedOutput, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            const DataCommunicator& r_comm = DataCommunicator::GetDefault();
            int rank = r_comm.Rank();

            // table can be is defined on all ranks
            LoggerOutput::Pointer p_output(new LoggerTableOutput(buffer, {"Time Step    ", "Iteration Number        ", "Convergence        ", "Is converged"}));
            Logger::AddOutput(p_output);

            // Header is printed immediately on rank 0, but in other ranks it will be printed only if there is output.
            std::stringstream reference_output;

            if (rank == 0) reference_output << "Time Step     Iteration Number         Convergence         Is converged " << std::endl;

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), reference_output.str().c_str());

            // Print time step on all ranks
            std::size_t time_step = 1;
            Logger("Time Step") << Logger::DistributedFilter::FromAllRanks() << time_step << std::endl;

            // now the header should appear on the remaining ranks too
            if (rank != 0) reference_output << "Time Step     Iteration Number         Convergence         Is converged " << std::endl;
            // add time step to logger output on all ranks
            reference_output << "1             ";

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), reference_output.str().c_str());

            // Messages with unknown labels are ignored, regardless of the rank
            Logger("Label") << "This log has a label which is not in the output columns and will not be printed in output " << std::endl;
            Logger("Label") << Logger::DistributedFilter::FromAllRanks() << "This should not be in the output, either" << std::endl;

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), reference_output.str().c_str());

            // reported results may be different on different ranks
            Logger("Time Step") << Logger::DistributedFilter::FromAllRanks() << 2 << std::endl;
            // printing from 0-9 because I do not want to bother with spacing...
            Logger("Iteration Number") << Logger::DistributedFilter::FromAllRanks() << rank % 10 << std::endl;
            Logger("Convergence") << 0.25 << std::endl; // only on root!
            // tryign out multiple lines
            Logger("Iteration Number") << Logger::DistributedFilter::FromAllRanks() << (rank + 1) % 10 << std::endl;
            Logger("Convergence") << Logger::DistributedFilter::FromAllRanks() << 0.05 << std::endl;
            Logger("Is converged") << Logger::DistributedFilter::FromAllRanks() << (rank == 0 ? "Yes" : "No") << std::endl;

            reference_output << std::endl << "2             " << rank     % 10 << "                        ";
            if (rank == 0) {
                reference_output << "0.25                ";
            }
            reference_output << std::endl << "              " << (rank+1) % 10 << "                        " << "0.05                ";
            if (rank == 0) {
                reference_output << "Yes          ";
            }
            else {
                reference_output << "No           ";
            }

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), reference_output.str().c_str());

            std::cout << "Final table from rank " << rank << ":" << std::endl << buffer.str() << std::endl;
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoAllRanks, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            const DataCommunicator& r_comm = DataCommunicator::GetDefault();
            int rank = r_comm.Rank();

            KRATOS_INFO_ALL_RANKS("TestInfo") << "Test info message";
            std::stringstream out;
            out << "Rank " << rank << ": TestInfo: Test info message";

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), out.str().c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoIfAllRanks, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            const DataCommunicator& r_comm = DataCommunicator::GetDefault();
            int rank = r_comm.Rank();
            std::stringstream out;

            KRATOS_INFO_IF_ALL_RANKS("TestInfo", false) << "Test info if false message";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), "");

            KRATOS_INFO_IF_ALL_RANKS("TestInfo", true) << "Test info if true message";
            out << "Rank " << rank << ": TestInfo: Test info if true message";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), out.str().c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoOnceAllRanks, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            const DataCommunicator& r_comm = DataCommunicator::GetDefault();
            int rank = r_comm.Rank();

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_ONCE_ALL_RANKS("TestInfo") << "Test info message - " << i;
            }
            std::stringstream out;
            out << "Rank " << rank << ": TestInfo: Test info message - 0";

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), out.str().c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoFirstAllRanks, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            const DataCommunicator& r_comm = DataCommunicator::GetDefault();
            int rank = r_comm.Rank();

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_FIRST_N_ALL_RANKS("TestInfo", 4) << ".";
            }

            std::stringstream out;
            for(std::size_t i = 0; i < 4; i++) {
                out << "Rank " << rank << ": TestInfo: .";
            }

            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), out.str().c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerNoPrefix, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            p_output->SetSeverity(LoggerMessage::Severity::DETAIL);
            p_output->SetOption(LoggerOutput::WARNING_PREFIX, false);
            p_output->SetOption(LoggerOutput::INFO_PREFIX, false);
            p_output->SetOption(LoggerOutput::DETAIL_PREFIX, false);
            p_output->SetOption(LoggerOutput::DEBUG_PREFIX, false);
            p_output->SetOption(LoggerOutput::TRACE_PREFIX, false);
            Logger::AddOutput(p_output);

            KRATOS_WARNING("TestWarning") << "Test message\n";
            KRATOS_INFO("TestInfo") << "Test message\n";
            KRATOS_DETAIL("TestDetail") << "Test message\n";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "TestWarning: Test message\nTestInfo: Test message\nTestDetail: Test message\n" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }


        KRATOS_TEST_CASE_IN_SUITE(LoggerPrefix, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            p_output->SetSeverity(LoggerMessage::Severity::DETAIL);
            p_output->SetOption(LoggerOutput::WARNING_PREFIX, true);
            p_output->SetOption(LoggerOutput::INFO_PREFIX, true);
            p_output->SetOption(LoggerOutput::DETAIL_PREFIX, true);
            p_output->SetOption(LoggerOutput::DEBUG_PREFIX, true);
            p_output->SetOption(LoggerOutput::TRACE_PREFIX, true);
            Logger::AddOutput(p_output);

            KRATOS_WARNING("TestWarning") << "Test message\n";
            KRATOS_INFO("TestInfo") << "Test message\n";
            KRATOS_DETAIL("TestDetail") << "Test message\n";

            std::string expected_output = DataCommunicator::GetDefault().Rank() == 0 ? "[WARNING] TestWarning: Test message\n[INFO] TestInfo: Test message\n[DETAIL] TestDetail: Test message\n" : "";
            KRATOS_CHECK_C_STRING_EQUAL(buffer.str().c_str(), expected_output.c_str());
        }

    }   // namespace Testing
}  // namespace Kratos.


