//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                     2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "add_cross_sections_to_python.h"
#include "custom_utilities/shell_cross_section.hpp"

#include "custom_elements/shell_elements/shell_thick_element_3D4N.hpp"
#include "custom_elements/shell_elements/shell_thin_element_3D3N.hpp"


namespace Kratos
{
namespace Python
{

namespace py = pybind11;

void Helper_SetCrossSectionsOnIntegrationPoints_Thin(ShellThinElement3D3N& el, const pybind11::list& seclist)
{
  int n = len(seclist);
  std::vector<ShellCrossSection::Pointer> shell_sec_list;
  for(int i = 0; i < n; i++) {
    auto p = pybind11::cast<ShellCrossSection::Pointer >( seclist[i] );
    shell_sec_list.push_back(p);
  }
  el.SetCrossSectionsOnIntegrationPoints(shell_sec_list);
}
void Helper_SetCrossSectionsOnIntegrationPoints_Thick(ShellThickElement3D4N& el, const pybind11::list& seclist)
{
  int n = len(seclist);
  std::vector<ShellCrossSection::Pointer> shell_sec_list;
    for(int i = 0; i < n; i++)
    {
        auto p = pybind11::cast<ShellCrossSection::Pointer>( seclist[i] );
        shell_sec_list.push_back(p);
    }
    el.SetCrossSectionsOnIntegrationPoints(shell_sec_list);
}

void AddCrossSectionsToPython(pybind11::module& m)
{

  py::class_<ShellCrossSection, ShellCrossSection::Pointer >(m,"ShellCrossSection")
      .def(py::init<>())
      .def("BeginStack", &ShellCrossSection::BeginStack)
      .def("AddPly", &ShellCrossSection::AddPly)
      .def("EndStack", &ShellCrossSection::EndStack)
      .def("SetOffset", &ShellCrossSection::SetOffset)
      .def("Clone", &ShellCrossSection::Clone)
      .def("NumberOfPlies", &ShellCrossSection::NumberOfPlies)
      .def("NumberOfIntegrationPointsAt", &ShellCrossSection::NumberOfIntegrationPointsAt)
      .def("SetConstitutiveLawAt", &ShellCrossSection::SetConstitutiveLawAt)
      .def("__repr__", &ShellCrossSection::Info )
      DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON_AS_POINTER(ShellCrossSection)
      DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON_AS_POINTER(ShellCrossSection)
      ;

  py::class_<Variable<ShellCrossSection::Pointer>, VariableData>(m,"ShellCrossSectionVariable")
      ;

  py::class_<ShellThinElement3D3N, ShellThinElement3D3N::Pointer, Element >(m,"ShellThinElement3D3N")
      .def("SetCrossSectionsOnIntegrationPoints", &Helper_SetCrossSectionsOnIntegrationPoints_Thin)
      ;

  py::class_<ShellThickElement3D4N, ShellThickElement3D4N::Pointer, Element >(m,"ShellThickElement3D4N")
      .def("SetCrossSectionsOnIntegrationPoints", &Helper_SetCrossSectionsOnIntegrationPoints_Thick)
      ;


}

}

}

