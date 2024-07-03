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

#if defined(KRATOS_COLORED_OUTPUT)
#include "utilities/color_utilities.h"
static std::string TEST_KYEL=KYEL;
static std::string TEST_RST=RST;
#else
static std::string TEST_KYEL="";
static std::string TEST_RST="";
#endif

namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(LoggerMessageStream, KratosCoreFastSuite)
        {
            LoggerMessage message("label");

            message << "Test message with number " << 12 << 'e' << "00";

            KRATOS_EXPECT_STREQ(message.GetLabel().c_str(), "label");
            if (Testing::GetDefaultDataCommunicator().Rank() == 0) KRATOS_EXPECT_STREQ(message.GetMessage().c_str(), "Test message with number 12e00");
            KRATOS_EXPECT_EQ(message.GetSeverity(), LoggerMessage::Severity::INFO);
            KRATOS_EXPECT_EQ(message.GetCategory(), LoggerMessage::Category::STATUS);
            KRATOS_EXPECT_EQ(message.GetLocation().GetFileName(), "Unknown");
            KRATOS_EXPECT_EQ(message.GetLocation().GetFunctionName(), "Unknown");
            KRATOS_EXPECT_EQ(message.GetLocation().GetLineNumber(), 0);

            message << LoggerMessage::Severity::DETAIL
                << LoggerMessage::Category::CRITICAL
                << KRATOS_CODE_LOCATION << std::endl;

            KRATOS_EXPECT_STREQ(message.GetMessage().c_str(), "Test message with number 12e00\n");
            KRATOS_EXPECT_EQ(message.GetSeverity(), LoggerMessage::Severity::DETAIL);
            KRATOS_EXPECT_EQ(message.GetCategory(), LoggerMessage::Category::CRITICAL);
            KRATOS_EXPECT_NE(message.GetLocation().GetFileName().find("test_logger.cpp"), std::string::npos);
            KRATOS_EXPECT_EQ(message.GetLocation().GetFunctionName(), KRATOS_CURRENT_FUNCTION);
            KRATOS_EXPECT_EQ(message.GetLocation().GetLineNumber(), 48);
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerOutput, KratosCoreFastSuite)
        {
            std::stringstream buffer;
            LoggerOutput output(buffer);

            LoggerMessage message("label");
            message << "Test message with number " << 12 << 'e' << "00";

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "label: Test message with number 12e00" : "";

            output.WriteMessage(message);
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStream, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            Logger("TestLabel") << "Test message with number " << 12 << 'e' << "00";

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestLabel: Test message with number 12e00" : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());

            Logger("TestDetail") << Logger::Severity::DETAIL << "This log has detailed severity and will not be printed in output "
                << Logger::Category::CRITICAL << std::endl;

            // The message has DETAIL severity and should not be written (check that nothing was added to the buffer)
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(CheckPoint, KratosCoreFastSuite)
        {
            std::stringstream buffer;
            LoggerOutput output(buffer);

            KRATOS_CHECK_POINT("TestCheckPoint") << "The value in check point is " << 3.14;

#if defined(KRATOS_ENABLE_CHECK_POINT)
            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestCheckPoint: The value in check point is 3.14" : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
#else
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), ""); // should print noting
#endif
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfo, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            KRATOS_INFO("TestInfo") << "Test info message";

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestInfo: Test info message" : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoIf, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            KRATOS_INFO_IF("TestInfo", true) << "Test info message";
            KRATOS_INFO_IF("TestInfo", false) << "This should not appear";

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestInfo: Test info message" : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoOnce, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_ONCE("TestInfo") << "Test info message - " << i;
            }

#ifdef KRATOS_DEBUG
            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestInfo: Test info message - 0" : "";
#else
            std::string expected_output = "";
#endif

            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoFirst, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_FIRST_N("TestInfo", 4) << ".";
            }

#ifdef KRATOS_DEBUG
            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestInfo: .TestInfo: .TestInfo: .TestInfo: ." : "";
#else
            std::string expected_output = "";
#endif
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamWarning, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            KRATOS_WARNING("TestWarning") << "Test warning message";

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? TEST_KYEL+"[WARNING] TestWarning: Test warning message"+TEST_RST : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamWarningIf, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            KRATOS_WARNING_IF("TestWarning", true) << "Test warning message";
            KRATOS_WARNING_IF("TestWarning", false) << "This should not appear";

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? TEST_KYEL+"[WARNING] TestWarning: Test warning message"+TEST_RST : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamWarningOnce, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_WARNING_ONCE("TestWarning") << "Test warning message - " << i;
            }

#ifdef KRATOS_DEBUG
            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? TEST_KYEL+"[WARNING] TestWarning: Test warning message - 0"+TEST_RST : "";
#else
            std::string expected_output = "";
#endif

            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamWarningFirst, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_WARNING_FIRST_N("TestWarning", 4) << ".";
            }

#ifdef KRATOS_DEBUG
            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? TEST_KYEL+"[WARNING] TestWarning: ."+TEST_RST+TEST_KYEL+"[WARNING] TestWarning: ."+TEST_RST+TEST_KYEL+"[WARNING] TestWarning: ."+TEST_RST+TEST_KYEL+"[WARNING] TestWarning: ."+TEST_RST : "";
#else
            std::string expected_output = "";
#endif

            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamDetail, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            p_output->SetSeverity(LoggerMessage::Severity::DETAIL);
            Logger::AddOutput(p_output);

            KRATOS_DETAIL("TestDetail") << "Test detail message";

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestDetail: Test detail message" : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamDetailIf, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            p_output->SetSeverity(LoggerMessage::Severity::DETAIL);
            Logger::AddOutput(p_output);

            KRATOS_DETAIL_IF("TestDetail", true) << "Test detail message";
            KRATOS_DETAIL_IF("TestDetail", false) << "This should not appear";

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestDetail: Test detail message" : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
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

#ifdef KRATOS_DEBUG
            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestDetail: Test detail message - 0" : "";
#else
            std::string expected_output = "";
#endif
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
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

#ifdef KRATOS_DEBUG
            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? "TestDetail: .TestDetail: .TestDetail: .TestDetail: ." : "";
#else
            std::string expected_output = "";
#endif
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerTableOutput, KratosCoreFastSuite)
        {
            int rank = Testing::GetDefaultDataCommunicator().Rank();

            static std::stringstream buffer;
            Parameters logger_settings(R"({
                "file_header" : "My Test",
                "label"       : "TEST",
                "columns" : [
                    {
                        "column_label" : "TIME_STEP",
                        "column_header": "Time Step"
                    },
                                        {
                        "column_label" : "IT_NUMBER",
                        "column_header": "Iteration Number"
                    },
                                        {
                        "column_label" : "CONVERGENCE",
                        "column_header": "Convergence"
                    },
                                        {
                        "column_label" : "IS_CONVERGED",
                        "column_header": "Is converged"
                    }
                ]
            })");
            LoggerOutput::Pointer p_output(new LoggerTableOutput(buffer, logger_settings));
            Logger::AddOutput(p_output);

            Logger("Label") << "This log has a label which is not in the output columns and will not be printed in output " << std::endl;
            double convergence = 0.00;
            for(std::size_t time_step = 1;time_step < 4; time_step++){
                for(int iteration_number = 1 ; iteration_number < 3; iteration_number++){
                    Logger("TEST.TIME_STEP") << time_step << std::endl;
                    convergence = 0.3 / (iteration_number * time_step);
                    Logger("TEST.IT_NUMBER") << iteration_number << std::endl;
                    Logger("TEST.CONVERGENCE") << convergence << std::endl;
                    if(convergence < 0.06) {
                        Logger("TEST.IS_CONVERGED") << "Yes" << std::endl;
                    }
                    else {
                        Logger("TEST.IS_CONVERGED") << "No" << std::endl;
                    }

                }

            }

            std::stringstream reference_output;

            if (rank == 0) {
                reference_output << "My Test" << std::endl ;
                reference_output << std::endl;
                reference_output << " Time Step  Iteration Number  Convergence  Is converged " << std::endl;
                reference_output << " ---------  ----------------  -----------  ------------ " << std::endl;
                reference_output << "     1              1             0.3           No      " << std::endl;
                reference_output << "     1              2             0.15          No      " << std::endl ;
                reference_output << "     2              1             0.15          No      " << std::endl ;
                reference_output << "     2              2            0.075          No      " << std::endl ;
                reference_output << "     3              1             0.1           No      " << std::endl ;
                reference_output << "     3              2             0.05          Yes     " << std::endl ;
            }

            KRATOS_EXPECT_STREQ(buffer.str().c_str(), reference_output.str().c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerTableDistributedOutput, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
            int rank = r_comm.Rank();
            Parameters logger_settings(R"({
                "file_header" : "My Test",
                "label"       : "TEST",
                "columns" : [
                    {
                        "column_label" : "TIME_STEP",
                        "column_header": "Time Step"
                    },
                                        {
                        "column_label" : "IT_NUMBER",
                        "column_header": "Iteration Number"
                    },
                                        {
                        "column_label" : "CONVERGENCE",
                        "column_header": "Convergence"
                    },
                                        {
                        "column_label" : "IS_CONVERGED",
                        "column_header": "Is converged"
                    }
                ]
            })");
            LoggerOutput::Pointer p_output(new LoggerTableOutput(buffer, logger_settings));
            Logger::AddOutput(p_output);

            // Header is printed immediately on rank 0, but in other ranks it will be printed only if there is output.
            std::stringstream reference_output;
            if (rank == 0){
                reference_output << "My Test" << std::endl ;
                reference_output << std::endl;
                reference_output << " Time Step  Iteration Number  Convergence  Is converged " << std::endl;
                reference_output << " ---------  ----------------  -----------  ------------ " << std::endl;
            }
            // Only in rank 0 should be printed
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), reference_output.str().c_str());
            if (rank != 0){
                reference_output << "My Test" << std::endl ;
                reference_output << std::endl;
                reference_output << " Time Step  Iteration Number  Convergence  Is converged " << std::endl;
                reference_output << " ---------  ----------------  -----------  ------------ " << std::endl;
            }

            double convergence = 0.00;
            for(std::size_t time_step = 1;time_step < 4; time_step++){
                for(int iteration_number = 1 ; iteration_number < 3; iteration_number++){
                    if(time_step !=1 && iteration_number!=1)
                    {
                        Logger("TEST.TIME_STEP") << time_step << std::endl;
                        convergence = 0.3 / (iteration_number * time_step);
                        Logger("TEST.IT_NUMBER") << iteration_number << std::endl;
                        Logger("TEST.CONVERGENCE") << convergence << std::endl;
                        if(convergence < 0.06) {
                            Logger("TEST.IS_CONVERGED") << "Yes" << std::endl;
                        }
                        else {
                            Logger("TEST.IS_CONVERGED") << "No" << std::endl;
                        }
                    } else {
                        // This should be printed to all ranks + the header if not done before
                        Logger("TEST.TIME_STEP") << Logger::DistributedFilter::FromAllRanks() << time_step << std::endl;
                        convergence = 0.3 / (iteration_number * time_step);
                        Logger("TEST.IT_NUMBER") << Logger::DistributedFilter::FromAllRanks() << iteration_number << std::endl;
                        Logger("TEST.CONVERGENCE") << Logger::DistributedFilter::FromAllRanks() << convergence << std::endl;
                        if(convergence < 0.06) {
                            Logger("TEST.IS_CONVERGED") << Logger::DistributedFilter::FromAllRanks() << "Yes" << std::endl;
                        }
                        else {
                            Logger("TEST.IS_CONVERGED") << Logger::DistributedFilter::FromAllRanks() << "No" << std::endl << Logger::DistributedFilter::FromRoot();
                        }
                    }
                }

            }

            // Messages with unknown labels are ignored, regardless of the rank
            Logger("Label") << "This log has a label which is not in the output columns and will not be printed in output " << std::endl;
            Logger("Label") << Logger::DistributedFilter::FromAllRanks() << "This should not be in the output, either" << std::endl;
            reference_output << "     1              1             0.3           No      " << std::endl;

            if (rank == 0) {
                reference_output << "     1              2             0.15          No      " << std::endl ;
                reference_output << "     2              1             0.15          No      " << std::endl ;
                reference_output << "     2              2            0.075          No      " << std::endl ;
                reference_output << "     3              1             0.1           No      " << std::endl ;
                reference_output << "     3              2             0.05          Yes     " << std::endl ;
            }

            KRATOS_EXPECT_STREQ(buffer.str().c_str(), reference_output.str().c_str());

        }
        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoAllRanks, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
            int rank = r_comm.Rank();

            KRATOS_INFO_ALL_RANKS("TestInfo") << "Test info message";
            std::stringstream out;
            out << "Rank " << rank << ": TestInfo: Test info message";

            KRATOS_EXPECT_STREQ(buffer.str().c_str(), out.str().c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoIfAllRanks, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
            int rank = r_comm.Rank();
            std::stringstream out;

            KRATOS_INFO_IF_ALL_RANKS("TestInfo", false) << "Test info if false message";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), "");

            KRATOS_INFO_IF_ALL_RANKS("TestInfo", true) << "Test info if true message";
            out << "Rank " << rank << ": TestInfo: Test info if true message";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), out.str().c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoOnceAllRanks, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_ONCE_ALL_RANKS("TestInfo") << "Test info message - " << i;
            }

            std::stringstream out;
#ifdef KRATOS_DEBUG
            out << "Rank " << Testing::GetDefaultDataCommunicator().Rank() << ": TestInfo: Test info message - 0";
#else
            out << "";
#endif

            KRATOS_EXPECT_STREQ(buffer.str().c_str(), out.str().c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(LoggerStreamInfoFirstAllRanks, KratosCoreFastSuite)
        {
            static std::stringstream buffer;
            LoggerOutput::Pointer p_output(new LoggerOutput(buffer));
            Logger::AddOutput(p_output);

            for(std::size_t i = 0; i < 10; i++) {
                KRATOS_INFO_FIRST_N_ALL_RANKS("TestInfo", 4) << ".";
            }

            std::stringstream out;
#ifdef KRATOS_DEBUG
            int rank = Testing::GetDefaultDataCommunicator().Rank();
            for(std::size_t i = 0; i < 4; i++) {
                out << "Rank " << rank << ": TestInfo: .";
            }
#else
            out << "";
#endif

            KRATOS_EXPECT_STREQ(buffer.str().c_str(), out.str().c_str());
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

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? TEST_KYEL+"TestWarning: Test message\n"+TEST_RST+"TestInfo: Test message\nTestDetail: Test message\n" : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
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

            std::string expected_output = Testing::GetDefaultDataCommunicator().Rank() == 0 ? TEST_KYEL+"[WARNING] TestWarning: Test message\n"+TEST_RST+"[INFO] TestInfo: Test message\n[DETAIL] TestDetail: Test message\n" : "";
            KRATOS_EXPECT_STREQ(buffer.str().c_str(), expected_output.c_str());
        }

    }   // namespace Testing
}  // namespace Kratos.


