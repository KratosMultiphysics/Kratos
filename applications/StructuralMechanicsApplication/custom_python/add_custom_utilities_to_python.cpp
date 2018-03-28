// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"


//Utilities
#include "custom_utilities/sprism_neighbours.hpp"
#include "custom_utilities/formfinding_IO_utility.h"


namespace Kratos
{
namespace Python
{

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

//     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//     typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_<FormfindingIOUtility>("FormfindingIOUtility", init<const ModelPart&, const Parameters>())
    .def("PrintModelPart",&FormfindingIOUtility::PrintModelPart)
    .def("ReadPrestressData",&FormfindingIOUtility::ReadPrestressData )
    .def("PrintPrestressData",&FormfindingIOUtility::PrintPrestressData )
    ;

}

}  // namespace Python.  

} // Namespace Kratos

