// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

// System includes
#include <limits>
#include <map>

/* External includes */
#include <filesystem>
#include <iostream>

/* Project includes */
#include "testing/testing.h"
#include "cpp_geomechanics_application.h"

namespace Kratos
{
    namespace Testing
    {
        void emptyProgress(double progress) {}
        void emptyLog(char *log) {}
        bool emptyCancel() {
            return false;
        }

        KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategy, KratosGeoMechanicsFastSuite)
        {
            auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30.gid";
            auto projectFile = "ProjectParameters.json";

            auto execute = KratosExecute();
            int status = execute.execute_flow_analysis(workingDirectory, projectFile, 3, 4, 0.1, "PorousDomain.Left_head", &emptyLog, &emptyProgress, &emptyLog, &emptyCancel);

            KRATOS_CHECK_EQUAL(status, 0);
        }

        KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategyTextualProgressReport, KratosGeoMechanicsFastSuite)
        {
            auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30.gid";
            auto projectFile = "ProjectParameters.json";

            auto execute = KratosExecute();
            
            bool firstMessageFound = false;
            bool finalMessageFound = false;
            int messageCount = 0;

            std::function<void(char*)> reportTextualProgress = [&firstMessageFound, &finalMessageFound, &messageCount](char* message) 
            {
                messageCount++;
                std::cout << "Captured: " << message << std::endl;

                if(strcmp(message, "Calculating head level 3m (1/12)") == 0) {
                    firstMessageFound = true;
                }

                if(strcmp(message, "Calculating head level 3.8m (9/12)") == 0) {
                    finalMessageFound = true;
                }
            };
            
            int status = execute.execute_flow_analysis(workingDirectory, projectFile, 3, 4, 0.1, "PorousDomain.Left_head", &emptyLog, &emptyProgress, reportTextualProgress, &emptyCancel);

            KRATOS_CHECK_EQUAL(status, 0);
            KRATOS_CHECK_EQUAL(firstMessageFound, true);
            KRATOS_CHECK_EQUAL(finalMessageFound, true);
            KRATOS_CHECK_EQUAL(messageCount, 9);
        }

        KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategyProgressReport, KratosGeoMechanicsFastSuite)
        {
            auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30.gid";
            auto projectFile = "ProjectParameters.json";

            auto execute = KratosExecute();
            
            bool startProgressFound = false;
            bool endProgressFound = false;
            int progressUpdates = 0;

            std::function<void(double)> reportProgress = [&startProgressFound, &endProgressFound, &progressUpdates](double progress) 
            {
                std::cout << "Progress: " << progress << std::endl;
                progressUpdates++;

                if(progress == 0.0) {
                    startProgressFound = true;
                }

                if(progress == 0.75) {
                    endProgressFound = true;
                }
            };
            
            int status = execute.execute_flow_analysis(workingDirectory, projectFile, 3, 4, 0.1, "PorousDomain.Left_head", &emptyLog, reportProgress, &emptyLog, &emptyCancel);

            KRATOS_CHECK_EQUAL(status, 0);
            KRATOS_CHECK_EQUAL(startProgressFound, true);
            KRATOS_CHECK_EQUAL(endProgressFound, true);
            KRATOS_CHECK_EQUAL(progressUpdates, 10);
        }
    }
}
