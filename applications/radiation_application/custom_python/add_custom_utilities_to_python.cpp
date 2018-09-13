//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//

// External includes 
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rad_face_utilities.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
	
namespace Python
{
void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;


    class_<RadFaceUtilities>(m,"RadFaceUtilities").def(init<>())
    .def("ConditionModelPart",&RadFaceUtilities::ConditionModelPart)
    ;



    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


}
}  // namespace Python.

} // Namespace Kratos

