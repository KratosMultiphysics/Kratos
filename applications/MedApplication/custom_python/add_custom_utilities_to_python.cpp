// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/med_testing_utilities.h"


namespace Kratos::Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MedTestingUtilities>(m,"MedTestingUtilities")
        .def_static("CheckModelPartsAreEqual", &MedTestingUtilities::CheckModelPartsAreEqual,
            py::arg("model_part_1"),
            py::arg("model_part_2"),
            py::arg("check_sub_model_parts")=true)
        .def_static("AddGeometriesFromElements", &MedTestingUtilities::AddGeometriesFromElements)
        .def_static("ComputeLength", &MedTestingUtilities::ComputeLength)
        .def_static("ComputeArea", &MedTestingUtilities::ComputeArea)
        .def_static("ComputeVolume", &MedTestingUtilities::ComputeVolume)
        .def_static("ComputeDomainSize", &MedTestingUtilities::ComputeDomainSize)
        ;
}

} // namespace Kratos::Python
