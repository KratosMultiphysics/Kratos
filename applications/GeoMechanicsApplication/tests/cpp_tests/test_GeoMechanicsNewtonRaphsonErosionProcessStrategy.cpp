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

#pragma once

// System includes
#include <limits>
#include <map>

/* External includes */
#include <filesystem>
#include <iostream>

/* Project includes */
#include "testing/testing.h"
#include "custom_workflows/dgeoflow.h"
#include "flow_stubs.h"

namespace Kratos
{
    namespace Testing
    {

        KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategy, KratosGeoMechanicsFastSuite)
        {
            auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30_fixed.gid";
            auto projectFile = "ProjectParameters.json";

            auto execute = KratosExecute();
            int status = execute.execute_flow_analysis(workingDirectory, projectFile, 10.5, 11.5, 0.1, "PorousDomain.Left_head", 
            &flow_stubs::emptyLog, &flow_stubs::emptyProgress, &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

            KRATOS_CHECK_EQUAL(status, 0);
        }

        KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategyTextualProgressReport, KratosGeoMechanicsFastSuite)
        {
            auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30_fixed.gid";
            auto projectFile = "ProjectParameters.json";

            auto execute = KratosExecute();
            
            bool firstMessageFound = false;
            bool finalMessageFound = false;
            int messageCount = 0;

            std::function<void(char*)> reportTextualProgress = [&firstMessageFound, &finalMessageFound, &messageCount](char* message) 
            {
                messageCount++;
                std::cout << "Captured: " << message << std::endl;

                if(strcmp(message, "Calculating head level 10.5m (1/12)") == 0) {
                    firstMessageFound = true;
                }

                if(strcmp(message, "Calculating head level 11.2m (8/12)") == 0) {
                    finalMessageFound = true;
                }
            };
            
            int status = execute.execute_flow_analysis(workingDirectory, projectFile, 10.5, 11.5, 0.1, "PorousDomain.Left_head", 
            &flow_stubs::emptyLog, &flow_stubs::emptyProgress, reportTextualProgress, &flow_stubs::emptyCancel);

            KRATOS_CHECK_EQUAL(status, 0);
            KRATOS_CHECK_EQUAL(firstMessageFound, true);
            KRATOS_CHECK_EQUAL(finalMessageFound, true);
            KRATOS_CHECK_EQUAL(messageCount, 8);
        }

        KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategyProgressReport, KratosGeoMechanicsFastSuite)
        {
            auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30_fixed.gid";
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

                if(progress - 0.666667 < 1e-6) {
                    endProgressFound = true;
                }
            };
            
            int status = execute.execute_flow_analysis(workingDirectory, projectFile, 10.5, 11.5, 0.1, "PorousDomain.Left_head",
            &flow_stubs::emptyLog, reportProgress, &flow_stubs::emptyLog, &flow_stubs::emptyCancel);

            KRATOS_CHECK_EQUAL(status, 0);
            KRATOS_CHECK_EQUAL(startProgressFound, true);
            KRATOS_CHECK_EQUAL(endProgressFound, true);
            KRATOS_CHECK_EQUAL(progressUpdates, 9);
        }
    }
}
