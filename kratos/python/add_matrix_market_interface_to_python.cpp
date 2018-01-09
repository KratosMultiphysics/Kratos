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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/matrix_market_interface.h"
#include "includes/ublas_interface.h"
#include "python/add_matrix_market_interface_to_python.h"

namespace Kratos
{

namespace Python
{
void  AddMatrixMarketInterfaceToPython()
{

    using namespace boost::python;

    def("ReadMatrixMarketMatrix", ReadMatrixMarketMatrix <Kratos::CompressedMatrix>);
    def("WriteMatrixMarketMatrix", WriteMatrixMarketMatrix <Kratos::CompressedMatrix>);

    def("ReadMatrixMarketVector", ReadMatrixMarketVector <Kratos::Vector>);
    def("WriteMatrixMarketVector", WriteMatrixMarketVector <Kratos::Vector>);

}

}  // namespace Python.

} // Namespace Kratos

