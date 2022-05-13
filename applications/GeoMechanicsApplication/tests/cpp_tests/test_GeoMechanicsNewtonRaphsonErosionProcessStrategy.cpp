

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

/* Project includes */
#include "testing/testing.h"
#include "cpp_geomechanics_application.h"

namespace Kratos
{
	namespace Testing
    {

    	KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategy, KratosGeoMechanicsFastSuite)
        {
            auto workingDirectory = "./SteadyStatePipeElementWithEmbankment/SteadyStatePipeElementWithEmbankment.gid";
            auto projectfile = "ProjectParameters.json";

            auto execute = KratosExecute();
    		execute.cpp_geomechanics(workingDirectory, projectfile);

        }
    }
}
