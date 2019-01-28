//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_linear_solvers_to_python.h"

//linear solvers
#ifdef INCLUDE_SUPERLU_MT
  #include "linear_system/linear_solvers/superlu_mt_direct_solver.hpp"
#else
  #include "linear_system/linear_solvers/superlu_direct_solver.hpp"
#endif

#ifdef INCLUDE_FEAST
  #include "linear_system/linear_solvers/feast_solver.hpp"
#endif

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void AddCustomLinearSolversToPython(pybind11::module& m)
{
  //base types
  typedef Kratos::Vector                                                     DenseVectorType;
  typedef Kratos::Matrix                                                     DenseMatrixType;
  typedef boost::numeric::ublas::vector<double>                             SparseVectorType;
  typedef UblasSpace<double, CompressedMatrix, SparseVectorType>             SparseSpaceType;
  typedef UblasSpace<double, DenseMatrixType, DenseVectorType>                LocalSpaceType;
  typedef DirectSolver<SparseSpaceType, LocalSpaceType>                     DirectSolverType;

#ifdef INCLUDE_SUPERLU_MT
  typedef SuperLUmtDirectSolver<SparseSpaceType, LocalSpaceType>   SuperLUmtDirectSolverType;

  py::class_<SuperLUmtDirectSolverType, typename SuperLUmtDirectSolverType::Pointer, DirectSolverType>
      (m, "superlu_direct")
      .def(py::init<>() )
      .def(py::init<Parameters>());
#else
  typedef SuperLUDirectSolver<SparseSpaceType, LocalSpaceType>       SuperLUDirectSolverType;

  py::class_<SuperLUDirectSolverType, typename SuperLUDirectSolverType::Pointer, DirectSolverType>
      (m, "superlu_direct")
      .def(py::init<>() )
      .def(py::init<Parameters>());

  //typedef SuperLUIterativeSolver<SparseSpaceType, LocalSpaceType> SuperLUIterativeSolverType;

  // py::class_<SuperLUIterativeSolverType, typename SuperLUIterativeSolverType::Pointer, SuperLUDirectSolverType>
  //     (m, "superlu_iterative")
  //     .def(py::init<>() )
  //     .def(py::init<Parameters>());
#endif

#ifdef INCLUDE_FEAST
  typedef FEASTEigenValueSolver<SparseSpaceType, LocalSpaceType>                    FEASTEigenValueSolverType;
  typedef LinearSolver<SparseSpaceType, LocalSpaceType>                                      LinearSolverType;
  typedef DenseVector<std::complex<double> >                                                ComplexVectorType;
  typedef DenseMatrix<std::complex<double> >                                                ComplexMatrixType;
  typedef UblasSpace<std::complex<double>, ComplexCompressedMatrix, ComplexVectorType> ComplexSparseSpaceType;
  typedef UblasSpace<std::complex<double>, ComplexMatrixType, ComplexVectorType>        ComplexLocalSpaceType;
  typedef LinearSolver<ComplexSparseSpaceType, ComplexLocalSpaceType>                 ComplexLinearSolverType;

  py::class_<FEASTEigenValueSolverType, FEASTEigenValueSolverType::Pointer, LinearSolverType>
      (m, "feast_eigen")
      .def(py::init<Parameters>() )
      .def(py::init<Parameters, ComplexLinearSolverType::Pointer>() )
      ;
#endif


}

}  // namespace Python.

} // Namespace Kratos
