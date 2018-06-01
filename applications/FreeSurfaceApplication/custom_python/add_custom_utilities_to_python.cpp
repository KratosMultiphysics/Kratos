//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/edge_data.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{


  void  AddCustomUtilitiesToPython(pybind11::module& pymodule)
  {
	using namespace pybind11;


		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
  		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    	typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
  
    	class_< MatrixContainer < 2, SparseSpaceType> > (pymodule,"MatrixContainer2D")
        .def(init< >())
    	.def("ConstructCSRVector", &MatrixContainer < 2, SparseSpaceType >::ConstructCSRVector)
    	.def("BuildCSRData", &MatrixContainer < 2, SparseSpaceType >::BuildCSRData)
    	.def("Clear", &MatrixContainer < 2, SparseSpaceType >::Clear)
    	;

    	class_< MatrixContainer < 3, SparseSpaceType> > (pymodule,"MatrixContainer3D")
        .def(init< >())
    	.def("ConstructCSRVector", &MatrixContainer < 3, SparseSpaceType >::ConstructCSRVector)
    	.def("BuildCSRData", &MatrixContainer < 3, SparseSpaceType >::BuildCSRData)
    	.def("Clear", &MatrixContainer < 3, SparseSpaceType >::Clear)
    	;


  }





}  // namespace Python.

} // Namespace Kratos
