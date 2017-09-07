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
#include <boost/python.hpp>
#include "includes/constitutive_law.h"
#include "includes/properties.h"

// Project includes
#include "add_cross_sections_to_python.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "custom_elements/shell_elements/shell_thick_element_3D4N.hpp"
#include "custom_elements/shell_elements/shell_thin_element_3D3N.hpp"


namespace Kratos
{
  namespace Python
  {

    using namespace boost::python;

    void Helper_SetCrossSectionsOnIntegrationPoints_Thin(ShellThinElement3D3N& el, const boost::python::list& seclist)
    {
      int n = len(seclist);
      std::vector<ShellCrossSection::Pointer> shell_sec_list;
      for(int i = 0; i < n; i++) {
	shell_sec_list.push_back(boost::python::extract<ShellCrossSection::Pointer>(seclist[i]));
      }
      el.SetCrossSectionsOnIntegrationPoints(shell_sec_list);
    }
    void Helper_SetCrossSectionsOnIntegrationPoints_Thick(ShellThickElement3D4N& el, const boost::python::list& seclist)
    {
      int n = len(seclist);
      std::vector<ShellCrossSection::Pointer> shell_sec_list;
      for(int i = 0; i < n; i++) {
	shell_sec_list.push_back(boost::python::extract<ShellCrossSection::Pointer>(seclist[i]));
      }
      el.SetCrossSectionsOnIntegrationPoints(shell_sec_list);
    }

    void AddCrossSectionsToPython()
    {

      class_<ShellCrossSection, ShellCrossSection::Pointer, boost::noncopyable >(
										 "ShellCrossSection", 
										 init<>())
	.def("BeginStack", &ShellCrossSection::BeginStack)
	.def("AddPly", &ShellCrossSection::AddPly)
	.def("EndStack", &ShellCrossSection::EndStack)
	.def("SetOffset", &ShellCrossSection::SetOffset)
	.def("Clone", &ShellCrossSection::Clone)
	.def("NumberOfPlies", &ShellCrossSection::NumberOfPlies)
	.def("NumberOfIntegrationPointsAt", &ShellCrossSection::NumberOfIntegrationPointsAt)
	.def("SetConstitutiveLawAt", &ShellCrossSection::SetConstitutiveLawAt)
	.def(self_ns::str(self))
	DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON_AS_POINTER(ShellCrossSection)
	DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON_AS_POINTER(ShellCrossSection)
	;

      class_<Variable<ShellCrossSection::Pointer> , bases<VariableData>, boost::noncopyable >("ShellCrossSectionVariable", no_init)
	;

      class_<ShellThinElement3D3N, ShellThinElement3D3N::Pointer, bases<Element>, boost::noncopyable >(
												       "ShellThinElement3D3N",
												       no_init)
	.def("SetCrossSectionsOnIntegrationPoints", &Helper_SetCrossSectionsOnIntegrationPoints_Thin)
	;

      class_<ShellThickElement3D4N, ShellThickElement3D4N::Pointer, bases<Element>, boost::noncopyable >(
													 "ShellThickElement3D4N",
													 no_init)
	.def("SetCrossSectionsOnIntegrationPoints", &Helper_SetCrossSectionsOnIntegrationPoints_Thick)
	;


    }

  }

}
