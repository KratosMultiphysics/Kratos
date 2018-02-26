//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
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

#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"



namespace Kratos
{

namespace Python
{


  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;


// 		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
// 		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
// 		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( UPPER_SURFACE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOWER_SURFACE )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( WAKE_NORMAL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( PROJECTION_MATRIX )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( UPPER_PROJECTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOWER_PROJECTION )

  }





}  // namespace Python.

} // Namespace Kratos
