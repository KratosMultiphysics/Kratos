//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#if defined(USE_EIGEN_MKL)
#include "custom_utilities/mkl_utilities.h"
#endif

namespace Kratos::Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

#if defined(USE_EIGEN_MKL)
    // Expose enum MKLThreadSetting
    py::enum_<MKLUtilities::MKLThreadSetting>(m, "MKLThreadSetting")
        .value("Minimal", MKLUtilities::MKLThreadSetting::Minimal)
        .value("Consistent", MKLUtilities::MKLThreadSetting::Consistent)
        .value("Do_nothing", MKLUtilities::MKLThreadSetting::Do_nothing)
        .value("Manual", MKLUtilities::MKLThreadSetting::Manual)
        .export_values();

    // Expose class MKLUtilities
    py::class_<MKLUtilities>(m, "MKLUtilities")
        .def_static("GetNumThreads", &MKLUtilities::GetNumThreads)
        .def_static("SetNumThreads", &MKLUtilities::SetNumThreads, py::arg("NumThreads"))
        .def_static("CheckThreadNumber", &MKLUtilities::CheckThreadNumber, py::arg("NumberOfMKLThreads"))
        .def_static("CheckThreadNumber", [](const MKLUtilities::MKLThreadSetting Setting) {
            return MKLUtilities::CheckThreadNumber(static_cast<int>(Setting));
        }, py::arg("ThreadSetting"))
        .def_static("ComputeMKLThreadCount", [](Kratos::Parameters Settings) {
            return MKLUtilities::ComputeMKLThreadCount(Settings);
        }, py::arg("Settings"))
        .def_static("ComputeMKLThreadCount", [](const int NumberOfMKLThreads) {
            return MKLUtilities::ComputeMKLThreadCount(NumberOfMKLThreads);
        }, py::arg("NumberOfMKLThreads"))
        .def_static("ComputeMKLThreadCount", [](const MKLUtilities::MKLThreadSetting Setting) {
            return MKLUtilities::ComputeMKLThreadCount(static_cast<int>(Setting));
        }, py::arg("ThreadSetting"))
        ;
#endif

}

} // namespace Kratos::Python
