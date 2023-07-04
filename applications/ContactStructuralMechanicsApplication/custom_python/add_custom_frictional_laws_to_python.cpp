// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_frictional_laws_to_python.h"

/* Utilities */
#include "custom_frictional_laws/frictional_law.h"
#include "custom_frictional_laws/frictional_law_with_derivative.h"
#include "custom_frictional_laws/tresca_frictional_law.h"
#include "custom_frictional_laws/coulomb_frictional_law.h"

namespace Kratos::Python
{
namespace py = pybind11;

void  AddCustomFrictionalLawsToPython(pybind11::module& m)
{
    // TODO: Reduce code duplication with templates

    /// Frictional laws
    /* Base class */
    using FrictionalLaw2D2N = FrictionalLawWithDerivative<2, 2, false, 2>;
    using FrictionalLaw3D3N = FrictionalLawWithDerivative<3, 3, false, 3>;
    using FrictionalLaw3D4N = FrictionalLawWithDerivative<3, 4, false, 4>;
    using FrictionalLaw3D3N4N = FrictionalLawWithDerivative<3, 3, false, 4>;
    using FrictionalLaw3D4N3N = FrictionalLawWithDerivative<3, 4, false, 3>;
    using FrictionalLaw2D2NNV = FrictionalLawWithDerivative<2, 2, true, 2>;
    using FrictionalLaw3D3NNV = FrictionalLawWithDerivative<3, 3, true, 3>;
    using FrictionalLaw3D4NNV = FrictionalLawWithDerivative<3, 4, true, 4>;
    using FrictionalLaw3D3N4NNV = FrictionalLawWithDerivative<3, 3, true, 4>;
    using FrictionalLaw3D4N3NNV = FrictionalLawWithDerivative<3, 4, true, 3>;

    /* Tresca */
    using TrescaFrictionalLaw2D2N = TrescaFrictionalLaw<2, 2, false, 2>;
    using TrescaFrictionalLaw3D3N = TrescaFrictionalLaw<3, 3, false, 3>;
    using TrescaFrictionalLaw3D4N = TrescaFrictionalLaw<3, 4, false, 4>;
    using TrescaFrictionalLaw3D3N4N = TrescaFrictionalLaw<3, 3, false, 4>;
    using TrescaFrictionalLaw3D4N3N = TrescaFrictionalLaw<3, 4, false, 3>;
    using TrescaFrictionalLaw2D2NNV = TrescaFrictionalLaw<2, 2, true, 2>;
    using TrescaFrictionalLaw3D3NNV = TrescaFrictionalLaw<3, 3, true, 3>;
    using TrescaFrictionalLaw3D4NNV = TrescaFrictionalLaw<3, 4, true, 4>;
    using TrescaFrictionalLaw3D3N4NNV = TrescaFrictionalLaw<3, 3, true, 4>;
    using TrescaFrictionalLaw3D4N3NNV = TrescaFrictionalLaw<3, 4, true, 3>;

    /* Coulomb */
    using CoulombFrictionalLaw2D2N = CoulombFrictionalLaw<2, 2, false, 2>;
    using CoulombFrictionalLaw3D3N = CoulombFrictionalLaw<3, 3, false, 3>;
    using CoulombFrictionalLaw3D4N = CoulombFrictionalLaw<3, 4, false, 4>;
    using CoulombFrictionalLaw3D3N4N = CoulombFrictionalLaw<3, 3, false, 4>;
    using CoulombFrictionalLaw3D4N3N = CoulombFrictionalLaw<3, 4, false, 3>;
    using CoulombFrictionalLaw2D2NNV = CoulombFrictionalLaw<2, 2, true, 2>;
    using CoulombFrictionalLaw3D3NNV = CoulombFrictionalLaw<3, 3, true, 3>;
    using CoulombFrictionalLaw3D4NNV = CoulombFrictionalLaw<3, 4, true, 4>;
    using CoulombFrictionalLaw3D3N4NNV = CoulombFrictionalLaw<3, 3, true, 4>;
    using CoulombFrictionalLaw3D4N3NNV = CoulombFrictionalLaw<3, 4, true, 3>;

    // Base class
    py::class_<FrictionalLaw, typename FrictionalLaw::Pointer>(m, "FrictionalLaw")
    .def(py::init<>())
    .def("GetFrictionCoefficient",&FrictionalLaw::GetFrictionCoefficient)
    .def("GetThresholdValue",&FrictionalLaw::GetThresholdValue)
    ;

    /* 2D 2N */
    // Base class
    py::class_<FrictionalLaw2D2N, typename FrictionalLaw2D2N::Pointer,FrictionalLaw>(m, "FrictionalLaw2D2N")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw2D2N, typename TrescaFrictionalLaw2D2N::Pointer, FrictionalLaw2D2N>(m, "TrescaFrictionalLaw2D2N")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw2D2N, typename CoulombFrictionalLaw2D2N::Pointer, FrictionalLaw2D2N>(m, "CoulombFrictionalLaw2D2N")
    .def(py::init<>())
    ;

    /* 3D 3N */
    // Base class
    py::class_<FrictionalLaw3D3N, typename FrictionalLaw3D3N::Pointer, FrictionalLaw>(m, "FrictionalLaw3D3N")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw3D3N, typename TrescaFrictionalLaw3D3N::Pointer, FrictionalLaw3D3N>(m, "TrescaFrictionalLaw3D3N")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw3D3N, typename CoulombFrictionalLaw3D3N::Pointer, FrictionalLaw3D3N>(m, "CoulombFrictionalLaw3D3N")
    .def(py::init<>())
    ;

    /* 3D 4N */
    // Base class
    py::class_<FrictionalLaw3D4N, typename FrictionalLaw3D4N::Pointer, FrictionalLaw>(m, "FrictionalLaw3D4N")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw3D4N, typename TrescaFrictionalLaw3D4N::Pointer, FrictionalLaw3D4N>(m, "TrescaFrictionalLaw3D4N")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw3D4N, typename CoulombFrictionalLaw3D4N::Pointer, FrictionalLaw3D4N>(m, "CoulombFrictionalLaw3D4N")
    .def(py::init<>())
    ;

    /* 3D 3N-4N */
    // Base class
    py::class_<FrictionalLaw3D3N4N, typename FrictionalLaw3D3N4N::Pointer,FrictionalLaw>(m, "FrictionalLaw3D3N4N")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw3D3N4N, typename TrescaFrictionalLaw3D3N4N::Pointer, FrictionalLaw3D3N4N>(m, "TrescaFrictionalLaw3D3N4N")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw3D3N4N, typename CoulombFrictionalLaw3D3N4N::Pointer, FrictionalLaw3D3N4N>(m, "CoulombFrictionalLaw3D3N4N")
    .def(py::init<>())
    ;

    /* 3D 4N-3N */
    // Base class
    py::class_<FrictionalLaw3D4N3N, typename FrictionalLaw3D4N3N::Pointer, FrictionalLaw>(m, "FrictionalLaw3D4N3N")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw3D4N3N, typename TrescaFrictionalLaw3D4N3N::Pointer, FrictionalLaw3D4N3N>(m, "TrescaFrictionalLaw3D4N3N")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw3D4N3N, typename CoulombFrictionalLaw3D4N3N::Pointer, FrictionalLaw3D4N3N>(m, "CoulombFrictionalLaw3D4N3N")
    .def(py::init<>())
    ;

    /* 2D 2N NV */
    // Base class
    py::class_<FrictionalLaw2D2NNV, typename FrictionalLaw2D2NNV::Pointer, FrictionalLaw>(m, "FrictionalLaw2D2NNV")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw2D2NNV, typename TrescaFrictionalLaw2D2NNV::Pointer, FrictionalLaw2D2NNV>(m, "TrescaFrictionalLaw2D2NNV")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw2D2NNV, typename CoulombFrictionalLaw2D2NNV::Pointer, FrictionalLaw2D2NNV>(m, "CoulombFrictionalLaw2D2NNV")
    .def(py::init<>())
    ;

    /* 3D 3N NV */
    // Base class
    py::class_<FrictionalLaw3D3NNV, typename FrictionalLaw3D3NNV::Pointer, FrictionalLaw>(m, "FrictionalLaw3D3NNV")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw3D3NNV, typename TrescaFrictionalLaw3D3NNV::Pointer, FrictionalLaw3D3NNV>(m, "TrescaFrictionalLaw3D3NNV")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw3D3NNV, typename CoulombFrictionalLaw3D3NNV::Pointer, FrictionalLaw3D3NNV>(m, "CoulombFrictionalLaw3D3NNV")
    .def(py::init<>())
    ;

    /* 3D 4N NV */
    // Base class
    py::class_<FrictionalLaw3D4NNV, typename FrictionalLaw3D4NNV::Pointer, FrictionalLaw>(m, "FrictionalLaw3D4NNV")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw3D4NNV, typename TrescaFrictionalLaw3D4NNV::Pointer, FrictionalLaw3D4NNV>(m, "TrescaFrictionalLaw3D4NNV")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw3D4NNV, typename CoulombFrictionalLaw3D4NNV::Pointer, FrictionalLaw3D4NNV>(m, "CoulombFrictionalLaw3D4NNV")
    .def(py::init<>())
    ;

    /* 3D 3N-4N NV */
    // Base class
    py::class_<FrictionalLaw3D3N4NNV, typename FrictionalLaw3D3N4NNV::Pointer, FrictionalLaw>(m, "FrictionalLaw3D3N4NNV")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw3D3N4NNV, typename TrescaFrictionalLaw3D3N4NNV::Pointer, FrictionalLaw3D3N4NNV>(m, "TrescaFrictionalLaw3D3N4NNV")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw3D3N4NNV, typename CoulombFrictionalLaw3D3N4NNV::Pointer, FrictionalLaw3D3N4NNV>(m, "CoulombFrictionalLaw3D3N4NNV")
    .def(py::init<>())
    ;

    /* 3D 4N-3N NV */
    // Base class
    py::class_<FrictionalLaw3D4N3NNV, typename FrictionalLaw3D4N3NNV::Pointer, FrictionalLaw>(m, "FrictionalLaw3D4N3NNV")
    .def(py::init<>())
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw3D4N3NNV, typename TrescaFrictionalLaw3D4N3NNV::Pointer, FrictionalLaw3D4N3NNV>(m, "TrescaFrictionalLaw3D4N3NNV")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw3D4N3NNV, typename CoulombFrictionalLaw3D4N3NNV::Pointer, FrictionalLaw3D4N3NNV>(m, "CoulombFrictionalLaw3D4N3NNV")
    .def(py::init<>())
    ;
}

}  // namespace Kratos::Python.
