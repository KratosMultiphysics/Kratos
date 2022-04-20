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
#include "cpp_geomechanics_application.hpp"

using namespace std;

namespace Kratos
{
	namespace Testing
    {

    	KRATOS_TEST_CASE_IN_SUITE(ErosionProcessStrategy, KratosGeoMechanicsFastSuite)
        {
            auto meshpath = "./SteadyStatePipeElementWithEmbankment/SteadyStatePipeElementWithEmbankment.gid/SteadyStatePipeElementWithEmbankment.mdpa";
            auto projectpath = "./SteadyStatePipeElementWithEmbankment/SteadyStatePipeElementWithEmbankment.gid/ProjectParameters.json";
            auto materialpath = "./SteadyStatePipeElementWithEmbankment/SteadyStatePipeElementWithEmbankment.gid/MaterialParameters.json";

            TCHAR pwd[MAX_PATH];
            GetCurrentDirectory(MAX_PATH, pwd);
            cout << pwd << std::endl;

            auto execute = CppExecution();
            execute.cpp_geomechanics(meshpath, projectpath, materialpath);

        }
    }
}
