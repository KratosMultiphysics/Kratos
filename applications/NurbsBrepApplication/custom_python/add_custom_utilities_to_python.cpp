//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/NurbBrepApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
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

#include "custom_utilities/NurbsBrepModeler.h"

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

    class_<NurbsBrepModeler, boost::noncopyable>("NurbsBrepModeler", init<Parameters&, ModelPart&>());
      //.def("SetUp", &NurbsBrepModeler::SetUp)
      //;
    //  .def("EvaluateShapeFunction", &NurbsShapeFunctionModeler::EvaluateShapeFunction)
    //  .def("EvaluateShapeFunctionSecondOrder", &NurbsShapeFunctionModeler::EvaluateShapeFunctionSecondOrder)
    //  ;
  }





}  // namespace Python.

} // Namespace Kratos
