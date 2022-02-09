// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_frictional_laws_to_python.h"

// Utilities
#include "custom_frictional_laws/frictional_law.h"
#include "custom_frictional_laws/frictional_law_with_derivative.h"
#include "custom_frictional_laws/tresca_frictional_law.h"
#include "custom_frictional_laws/coulomb_frictional_law.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void  AddCustomFrictionalLawsToPython(pybind11::module& m)
{
    /// Frictional laws
    /* Base class */
    typedef FrictionalLawWithDerivative<2,2,false,2> FrictionalLaw2D2N;
    typedef FrictionalLawWithDerivative<3,3,false,3> FrictionalLaw3D3N;
    typedef FrictionalLawWithDerivative<3,4,false,4> FrictionalLaw3D4N;
    typedef FrictionalLawWithDerivative<3,3,false,4> FrictionalLaw3D3N4N;
    typedef FrictionalLawWithDerivative<3,4,false,3> FrictionalLaw3D4N3N;
    typedef FrictionalLawWithDerivative<2,2,true,2> FrictionalLaw2D2NNV;
    typedef FrictionalLawWithDerivative<3,3,true,3> FrictionalLaw3D3NNV;
    typedef FrictionalLawWithDerivative<3,4,true,4> FrictionalLaw3D4NNV;
    typedef FrictionalLawWithDerivative<3,3,true,4> FrictionalLaw3D3N4NNV;
    typedef FrictionalLawWithDerivative<3,4,true,3> FrictionalLaw3D4N3NNV;

    /* Tresca */
    typedef TrescaFrictionalLaw<2,2,false,2> TrescaFrictionalLaw2D2N;
    typedef TrescaFrictionalLaw<3,3,false,3> TrescaFrictionalLaw3D3N;
    typedef TrescaFrictionalLaw<3,4,false,4> TrescaFrictionalLaw3D4N;
    typedef TrescaFrictionalLaw<3,3,false,4> TrescaFrictionalLaw3D3N4N;
    typedef TrescaFrictionalLaw<3,4,false,3> TrescaFrictionalLaw3D4N3N;
    typedef TrescaFrictionalLaw<2,2,true,2> TrescaFrictionalLaw2D2NNV;
    typedef TrescaFrictionalLaw<3,3,true,3> TrescaFrictionalLaw3D3NNV;
    typedef TrescaFrictionalLaw<3,4,true,4> TrescaFrictionalLaw3D4NNV;
    typedef TrescaFrictionalLaw<3,3,true,4> TrescaFrictionalLaw3D3N4NNV;
    typedef TrescaFrictionalLaw<3,4,true,3> TrescaFrictionalLaw3D4N3NNV;

    /* Coulomb */
    typedef CoulombFrictionalLaw<2,2,false,2> CoulombFrictionalLaw2D2N;
    typedef CoulombFrictionalLaw<3,3,false,3> CoulombFrictionalLaw3D3N;
    typedef CoulombFrictionalLaw<3,4,false,4> CoulombFrictionalLaw3D4N;
    typedef CoulombFrictionalLaw<3,3,false,4> CoulombFrictionalLaw3D3N4N;
    typedef CoulombFrictionalLaw<3,4,false,3> CoulombFrictionalLaw3D4N3N;
    typedef CoulombFrictionalLaw<2,2,true,2> CoulombFrictionalLaw2D2NNV;
    typedef CoulombFrictionalLaw<3,3,true,3> CoulombFrictionalLaw3D3NNV;
    typedef CoulombFrictionalLaw<3,4,true,4> CoulombFrictionalLaw3D4NNV;
    typedef CoulombFrictionalLaw<3,3,true,4> CoulombFrictionalLaw3D3N4NNV;
    typedef CoulombFrictionalLaw<3,4,true,3> CoulombFrictionalLaw3D4N3NNV;

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

}  // namespace Python.

} // Namespace Kratos

