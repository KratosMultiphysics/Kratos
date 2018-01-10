//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/NurbsShapeFunctionModeler.h"

namespace Kratos
{

namespace Python
{


  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;


		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

		class_<NurbsShapeFunctionModeler, boost::noncopyable>("NurbsShapeFunctionModeler", init<>())
			.def("SetUp", &NurbsShapeFunctionModeler::SetUp)
			.def("EvaluateShapeFunction", &NurbsShapeFunctionModeler::EvaluateShapeFunction)
			.def("EvaluateShapeFunctionSecondOrder", &NurbsShapeFunctionModeler::EvaluateShapeFunctionSecondOrder)
			;
  }





}  // namespace Python.

} // Namespace Kratos
