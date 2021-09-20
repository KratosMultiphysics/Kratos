//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

// Schemes
#include "custom_strategies/velocity_height_residual_based_bdf_scheme.h"
#include "custom_strategies/shallow_water_residual_based_bdf_scheme.h"
#include "custom_strategies/flux_corrected_shallow_water_scheme.h"

namespace Kratos
{

namespace Python
{

  void AddCustomStrategiesToPython(pybind11::module& m)
  {
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // Schemes
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    typedef VelocityHeightResidualBasedBDFScheme<SparseSpaceType, LocalSpaceType> VelocityHeightResidualBasedBDFSchemeType;
    py::class_<VelocityHeightResidualBasedBDFSchemeType, VelocityHeightResidualBasedBDFSchemeType::Pointer, BaseSchemeType>
    (m, "VelocityHeightResidualBasedBDFScheme")
    .def(py::init<std::size_t>())
    ;

    typedef ShallowWaterResidualBasedBDFScheme<SparseSpaceType, LocalSpaceType> ShallowWaterResidualBasedBDFSchemeType;
    py::class_<ShallowWaterResidualBasedBDFSchemeType, ShallowWaterResidualBasedBDFSchemeType::Pointer, BaseSchemeType>
    (m, "ShallowWaterResidualBasedBDFScheme")
    .def(py::init<std::size_t>())
    .def(py::init<std::size_t, bool>())
    ;

    typedef FluxCorrectedShallowWaterScheme<SparseSpaceType, LocalSpaceType> FluxCorrectedShallowWaterSchemeType;
    py::class_<FluxCorrectedShallowWaterSchemeType, FluxCorrectedShallowWaterSchemeType::Pointer, BaseSchemeType>
    (m, "FluxCorrectedShallowWaterScheme")
    .def(py::init<std::size_t>())
    .def(py::init<std::size_t, bool>())
    .def(py::init<Parameters>())
    ;

  }

}  // namespace Python.

} // Namespace Kratos
