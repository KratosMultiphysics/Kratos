// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes
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

using namespace pybind11;

void Helper_SetCrossSectionsOnIntegrationPoints_Thin(ShellThinElement3D3N& el, const pybind11::list& seclist)
{
    int n = len(seclist);
    std::vector<ShellCrossSection::Pointer> shell_sec_list;
    for(int i = 0; i < n; i++)
    {
        auto p = pybind11::cast<ShellCrossSection::Pointer >( seclist[i] );
        shell_sec_list.push_back(p);
//         shell_sec_list.push_back(boost::python::extract<ShellCrossSection::Pointer>(seclist[i]));
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
//         shell_sec_list.push_back(boost::python::extract<ShellCrossSection::Pointer>(seclist[i]));
    }
    el.SetCrossSectionsOnIntegrationPoints(shell_sec_list);
}

void AddCrossSectionsToPython(pybind11::module& m)
{

    class_<ShellCrossSection, ShellCrossSection::Pointer >(m,"ShellCrossSection")
    .def(init<>())
    .def("BeginStack", &ShellCrossSection::BeginStack)
    .def("AddPly", &ShellCrossSection::AddPly)
    .def("EndStack", &ShellCrossSection::EndStack)
    // .def("SetOffset", &ShellCrossSection::SetOffset)
    .def("Clone", &ShellCrossSection::Clone)
    .def("NumberOfPlies", &ShellCrossSection::NumberOfPlies)
    .def("NumberOfIntegrationPointsAt", &ShellCrossSection::NumberOfIntegrationPointsAt)
    .def("SetConstitutiveLawAt", &ShellCrossSection::SetConstitutiveLawAt)
    // .def("__repr__", &ShellCrossSection::Info )
    DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON_AS_POINTER(ShellCrossSection)
    DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON_AS_POINTER(ShellCrossSection)
    ;

    class_<Variable<ShellCrossSection::Pointer>,VariableData >(m,"ShellCrossSectionVariable")
    ;

    class_<ShellThinElement3D3N, ShellThinElement3D3N::Pointer, Element >(m,"ShellThinElement3D3N")
    .def("SetCrossSectionsOnIntegrationPoints", &Helper_SetCrossSectionsOnIntegrationPoints_Thin)
    ;

    class_<ShellThickElement3D4N, ShellThickElement3D4N::Pointer, Element >(m,"ShellThickElement3D4N")
    .def("SetCrossSectionsOnIntegrationPoints", &Helper_SetCrossSectionsOnIntegrationPoints_Thick)
    ;


}

}

}
