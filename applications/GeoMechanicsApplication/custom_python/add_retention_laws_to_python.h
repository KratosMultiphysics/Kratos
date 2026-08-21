// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra,
//                   Richard Faasse 
//

#pragma once

#include <pybind11/pybind11.h>

namespace Kratos::Python
{

void AddRetentionLawsToPython(const pybind11::module& rModule);

} // namespace Kratos::Python.
