/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"

// Application includes
#include "custom_utilities/mpi_dem_search.h"
#include "custom_utilities/mpi_utilities.h"

namespace Kratos {
namespace Python {

using namespace pybind11;

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    typedef MPI_DEMSearch   DemSearchType;
    typedef MpiUtilities    MpiUtilitiesType;

    class_<DemSearchType, DemSearchType::Pointer, SpatialSearch>(m, "MPI_DEMSearch")
        .def(init<Communicator&>())
        ;

    class_<MpiUtilitiesType, MpiUtilitiesType::Pointer>(m, "MpiUtilities")
        .def(init<>())
        .def("Repartition",             &MpiUtilitiesType::ParallelPartitioning)
        .def("TransferModelElements",   &MpiUtilitiesType::MigrateElements)
        .def("TransferModelNodes",      &MpiUtilitiesType::MigrateNodes)
        .def("CalculateModelNewIds",    &MpiUtilitiesType::CalculateModelNewIds)
        .def("CalculateElementsNewId",  &MpiUtilitiesType::CalculateElementsNewId)
        .def("CalculateNodesNewId",     &MpiUtilitiesType::CalculateNodesNewId)
        .def("CalculateConditionsNewId",&MpiUtilitiesType::CalculateConditionsNewId)
        ;
}

}  // namespace Python.
} // Namespace Kratos

#endif // KRATOS_PYTHON defined
