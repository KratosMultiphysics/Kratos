//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                 October 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

// Utilities
#include "custom_utilities/properties_layout.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{

  namespace py = pybind11;

  typedef TableKeyVariables<double,double>  TableKeyScalarVariablesType;

  py::class_<PropertiesLayout, PropertiesLayout::Pointer>(m,"PropertiesLayout")
      .def(py::init<>())
      .def("Clone", &PropertiesLayout::Clone)
      .def("RegisterTable", &PropertiesLayout::RegisterTable)
      .def("__repr__", &PropertiesLayout::Info)
      DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON(PropertiesLayout)
      DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON(PropertiesLayout)
      ;

  py::class_<Variable<PropertiesLayout>, VariableData>(m,"PropertiesLayoutVariable")
      ;

}

}  // namespace Python.

} // Namespace Kratos
