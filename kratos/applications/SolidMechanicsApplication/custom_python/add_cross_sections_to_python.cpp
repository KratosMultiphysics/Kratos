/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
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

#include "custom_elements/shell_thick_element_3D4N.hpp"
#include "custom_elements/shell_thin_element_3D3N.hpp"

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