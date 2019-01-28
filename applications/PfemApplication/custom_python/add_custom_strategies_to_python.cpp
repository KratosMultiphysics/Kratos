//
//   Project Name:        KratosPfemApplication           $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "utilities/openmp_utils.h"
#include "custom_python/add_custom_strategies_to_python.h"

// Solution schemes
#include "custom_solvers/solution_schemes/ale_solution_scheme.hpp"


namespace Kratos
{

namespace Python
{

//base types
typedef Kratos::Vector                                                                        DenseVectorType;
typedef Kratos::Matrix                                                                        DenseMatrixType;
typedef boost::numeric::ublas::vector<double>                                                     SparseVectorType;
typedef boost::numeric::ublas::matrix<double>                                                     SparseMatrixType;
typedef UblasSpace<double, CompressedMatrix, SparseVectorType>                                     SparseSpaceType;
typedef UblasSpace<double, DenseMatrixType, DenseVectorType>                                        LocalSpaceType;

void  AddCustomStrategiesToPython(pybind11::module& m)
{

  namespace py = pybind11;

  // Solution scheme types
  typedef DynamicScheme<SparseSpaceType, LocalSpaceType>                                DynamicSchemeType;
  typedef AleSolutionScheme<SparseSpaceType, LocalSpaceType>                            AleSolutionSchemeType;

  // Time integration methods for vectors
  typedef array_1d<double, 3>                                                                      VectorType;
  typedef Variable<VectorType>                                                             VariableVectorType;
  typedef TimeIntegrationMethod<VariableVectorType, VectorType>               TimeIntegrationMethodVectorType;

  typedef std::vector<TimeIntegrationMethodVectorType::Pointer>                  TimeVectorIntegrationMethods;

  // Time integration methods for scalars
  typedef Variable<double>                                                                 VariableScalarType;
  typedef TimeIntegrationMethod<VariableScalarType, double>                   TimeIntegrationMethodScalarType;

  typedef std::vector<TimeIntegrationMethodScalarType::Pointer>                  TimeScalarIntegrationMethods;

  //*************************SHCHEME CLASSES****************************

  // Dynamic Scheme Type
  py::class_<AleSolutionSchemeType, typename AleSolutionSchemeType::Pointer, DynamicSchemeType>(m,"AleDynamicScheme")
      .def(py::init<TimeVectorIntegrationMethods&, TimeScalarIntegrationMethods&>())
      .def(py::init<TimeVectorIntegrationMethods&, TimeScalarIntegrationMethods&, Flags&>())
      ;

}

}  // namespace Python.

} // Namespace Kratos
