//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_integration_to_python.h"

#include "integration/cad_integration_domain.h"

namespace Kratos
{

namespace Python
{

void  AddIntegrationToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def("CreateIntegrationDomain", &CadIntegrationDomain::CreateIntegrationDomain);
}

}  // namespace Python.

} // Namespace Kratos


