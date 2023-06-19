//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Internal includes
#include "add_custom_utilities_to_python.h"

// Project includes
#include "custom_utilities/model_subdivision_utilities.h"


namespace Kratos
{
namespace Python
{


void AddCustomUtilitiesToPython(pybind11::module& rModule)
{
    pybind11::class_<Wind::ModelSubdivisionUtilities>(rModule, "ModelSubdivisionUtilities")
        .def_static("SortNodesBySlabs", &Wind::ModelSubdivisionUtilities::SortNodesBySlabs)
        ;
}


} // namespace Python
} // namespace Kratos