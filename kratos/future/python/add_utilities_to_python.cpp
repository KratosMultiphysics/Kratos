
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "includes/model_part.h"
#include "includes/define_python.h"

// Future extensions
#include "future/python/add_utilities_to_python.h"
#include "future/utilities/csr_utilities.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddUtilitiesToPython(py::module& m)
{
    m.def_submodule("CsrUtilities")
        .def("GetEquationIdCsrIndices", &CsrUtilities::GetEquationIdCsrIndices<ModelPart::ElementsContainerType, CsrMatrix<>>)
    ;
}

}  // namespace Kratos::Future::Python.

