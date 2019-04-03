//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//



// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/matrix_market_interface.h"
#include "includes/ublas_interface.h"
#include "python/add_matrix_market_interface_to_python.h"

namespace Kratos
{

namespace Python
{
void  AddMatrixMarketInterfaceToPython(pybind11::module& m)
{

    namespace py = pybind11;

    m.def("ReadMatrixMarketMatrix", ReadMatrixMarketMatrix <Kratos::CompressedMatrix>);
    m.def("WriteMatrixMarketMatrix", WriteMatrixMarketMatrix <Kratos::CompressedMatrix>);

    m.def("ReadMatrixMarketVector", ReadMatrixMarketVector <Kratos::Vector>);
    m.def("WriteMatrixMarketVector", WriteMatrixMarketVector <Kratos::Vector>);

}

}  // namespace Python.

} // Namespace Kratos

