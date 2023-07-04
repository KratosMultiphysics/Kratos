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

/**
 * @brief Registers frictional laws in the given module.
 * @param m the module to register the laws in
 * @param rEndName the name to append to the law names
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster>
void RegisterFrictionalLaws(pybind11::module& m, const std::string& rEndName)
{
    // Base class
    using FrictionalLawWithDerivativeType = Kratos::FrictionalLawWithDerivative<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>;
    std::string name = "FrictionalLaw" + rEndName;
    pybind11::class_<FrictionalLawWithDerivativeType, typename FrictionalLawWithDerivativeType::Pointer, Kratos::FrictionalLaw>(m, name.c_str())
    .def(pybind11::init<>())
    ;

    // Tresca frictional law
    using TrescaFrictionalLawType = Kratos::TrescaFrictionalLaw<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>;
    name = "TrescaFrictionalLaw" + rEndName;
    pybind11::class_<TrescaFrictionalLawType, typename TrescaFrictionalLawType::Pointer, FrictionalLawWithDerivativeType>(m, name.c_str())
    .def(pybind11::init<>())
    ;

    // Coulomb frictional law
    using CoulombFrictionalLawType = Kratos::CoulombFrictionalLaw<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>;
    name = "CoulombFrictionalLaw" + rEndName;
    pybind11::class_<CoulombFrictionalLawType, typename CoulombFrictionalLawType::Pointer, FrictionalLawWithDerivativeType>(m, name.c_str())
    .def(pybind11::init<>())
    ;
}

namespace Kratos::Python
{
namespace py = pybind11;

void  AddCustomFrictionalLawsToPython(pybind11::module& m)
{
    RegisterFrictionalLaws<2, 2, false, 2>(m, "2D2N");
    RegisterFrictionalLaws<3, 3, false, 3>(m, "3D3N");
    RegisterFrictionalLaws<3, 4, false, 4>(m, "3D4N");
    RegisterFrictionalLaws<3, 3, false, 4>(m, "3D3N4N");
    RegisterFrictionalLaws<3, 4, false, 3>(m, "3D4N3N");
    RegisterFrictionalLaws<2, 2, true, 2>(m, "2D2NNV");
    RegisterFrictionalLaws<3, 3, true, 3>(m, "3D3NNV");
    RegisterFrictionalLaws<3, 4, true, 4>(m, "3D4NNV");
    RegisterFrictionalLaws<3, 3, true, 4>(m, "3D3N4NNV");
    RegisterFrictionalLaws<3, 4, true, 3>(m, "3D4N3NNV");
}

}  // namespace Kratos::Python.
