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
        void emptyLog(char *log) {}
        bool emptyCancel() {
            return false;
        }

        KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategy, KratosGeoMechanicsFastSuite)
        {
            auto workingDirectory = "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/HeightAquiferD10L30.gid";
            auto projectFile = "ProjectParameters.json";

            auto execute = KratosExecute();
            execute.execute_flow_analysis(workingDirectory, projectFile, 3, 4, 0.1, "PorousDomain.Left_head", &emptyLog, &emptyLog, &emptyCancel);
        }
    }
}
