//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/default_spaces.h"

// Schemes
#include "custom_strategies/schemes/residual_based_bossak_scheme.hpp"

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
  typedef TDefaultSparseSpace<double>                         SparseSpaceType;
  typedef TDefaultDenseSpace<double>                                    LocalSpaceType;
  typedef Scheme< SparseSpaceType, LocalSpaceType >                                 SchemeType;

  //custom scheme types
  typedef ResidualBasedBossakScheme< SparseSpaceType, LocalSpaceType >              ResidualBasedBossakSchemeType;


  //*************************SHCHEME CLASSES****************************


  // Residual Based Bossak Scheme Type
  py::class_<ResidualBasedBossakSchemeType, typename ResidualBasedBossakSchemeType::Pointer, SchemeType>
      (m,"ResidualBasedBossakScheme")
      .def(py::init< double , double >())
      .def("Initialize", &ResidualBasedBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
      ;


}

}  // namespace Python.

} // Namespace Kratos

