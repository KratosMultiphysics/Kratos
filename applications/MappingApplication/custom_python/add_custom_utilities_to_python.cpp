
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/mapper_utilities.h"
#include "custom_utilities/mapping_intersection_utilities.h"
#include "custom_utilities/iga_mapping_intersection_utilities.h"


namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto m_mapper_utilities = m.def_submodule("MapperUtilities");
    auto m_mapping_intersection_utils = m.def_submodule("MappingIntersectionUtilities");
    auto m_iga_mapping_intersection_utils = m.def_submodule("IgaMappingIntersectionUtilities");

    m_mapper_utilities.def("SaveCurrentConfiguration", &MapperUtilities::SaveCurrentConfiguration);
    m_mapper_utilities.def("RestoreCurrentConfiguration", &MapperUtilities::RestoreCurrentConfiguration);
    m_mapping_intersection_utils.def("FindIntersection1DGeometries2D", &MappingIntersectionUtilities::FindIntersection1DGeometries2D);
    m_mapping_intersection_utils.def("CreateQuadraturePointsCoupling1DGeometries2D", &MappingIntersectionUtilities::CreateQuadraturePointsCoupling1DGeometries2D);
}

}  // namespace Python.
} // Namespace Kratos
