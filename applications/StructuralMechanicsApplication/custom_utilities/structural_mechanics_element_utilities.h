//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined( KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos {
namespace StructuralMechanicsElementUtilities {

bool ComputeLumpedMassMatrix(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

bool HasRayleighDamping(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

double GetRayleighAlpha(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

double GetRayleighBeta(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

void CalculateRayleighDampingMatrix(
    Element& rElement,
    Element::MatrixType& rDampingMatrix,
    /*const*/ ProcessInfo& rCurrentProcessInfo,
    const std::size_t MatrixSize);

} // namespace StructuralMechanicsElementUtilities.
}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED  defined
