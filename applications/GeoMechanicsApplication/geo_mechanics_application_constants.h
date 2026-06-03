// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "includes/define.h"

namespace Kratos
{
/* Constants used in GeoMechanicsApplication */
// The size type definition
using SizeType = std::size_t;

// Static definition of the dimension
constexpr SizeType N_DIM_3D = 3;
constexpr SizeType N_DIM_2D = 2;
constexpr SizeType N_DIM_1D = 1;

// Limits
constexpr double TINY  = 1.0e-60;
constexpr double LARGE = 1.0e10;

// factor for pore pressure calculations: 1: mechanical sign convention, -1: soil mechanics sign convention
constexpr double PORE_PRESSURE_SIGN_FACTOR = 1.0;

// Static definition of the size of stress tensor (n x n)
constexpr SizeType STRESS_TENSOR_SIZE_2D = 3;
constexpr SizeType STRESS_TENSOR_SIZE_3D = 3;

// Static definition of the VoigtSize
constexpr SizeType VOIGT_SIZE_3D              = 6;
constexpr SizeType VOIGT_SIZE_2D_PLANE_STRESS = 3;
constexpr SizeType VOIGT_SIZE_2D_PLANE_STRAIN = 4;
constexpr SizeType VOIGT_SIZE_2D_AXISYMMETRIC = 4;
constexpr SizeType VOIGT_SIZE_2D_INTERFACE    = 2;
constexpr SizeType VOIGT_SIZE_3D_INTERFACE    = 3;

// stress/strain vector indices (3D interface):
enum class indexStress3DInterface : int {
    INDEX_3D_INTERFACE_XZ,
    INDEX_3D_INTERFACE_YZ,
    INDEX_3D_INTERFACE_ZZ
};

enum class IsDiffOrderElement : int { Yes, No };

enum class PlasticityStatus : int {
    ELASTIC,
    TENSION_APEX,
    TENSION_CUT_OFF,
    TENSION_MOHR_COULOMB_CORNER,
    MOHR_COULOMB_FAILURE
};

enum class DrainageType : int { DRAINED, FULLY_COUPLED, UNDRAINED, CONSTANT_WATER_PRESSURE };

} // namespace Kratos