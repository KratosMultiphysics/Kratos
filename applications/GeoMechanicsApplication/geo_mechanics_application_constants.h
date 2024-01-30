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
    /* Contants used in GeoMechanicsApplication */
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
    constexpr SizeType VOIGT_SIZE_3D = 6;
    constexpr SizeType VOIGT_SIZE_2D_PLANE_STRESS = 3;
    constexpr SizeType VOIGT_SIZE_2D_PLANE_STRAIN = 4;
    constexpr SizeType VOIGT_SIZE_2D_AXISYMMETRIC = 4;
    constexpr SizeType VOIGT_SIZE_2D_INTERFACE    = 2;
    constexpr SizeType VOIGT_SIZE_3D_INTERFACE    = 3;

    // DOF indices (3D):
    enum indexDOF3D: int{ INDEX_X,
                          INDEX_Y,
                          INDEX_Z };

    // stress/strain vector indices (3D):
    enum indexStress3D: int{ INDEX_3D_XX,
                             INDEX_3D_YY,
                             INDEX_3D_ZZ,
                             INDEX_3D_XY,
                             INDEX_3D_YZ,
                             INDEX_3D_XZ };

    // stress/strain vector indices (2D plane strain):
    enum indexStress2DPlaneStrain: int{ INDEX_2D_PLANE_STRAIN_XX,
                                        INDEX_2D_PLANE_STRAIN_YY,
                                        INDEX_2D_PLANE_STRAIN_ZZ,
                                        INDEX_2D_PLANE_STRAIN_XY };

    // stress/strain vector indices (2D plane stress):
    enum indexStress2DPlaneStress: int{ INDEX_2D_PLANE_STRESS_XX,
                                        INDEX_2D_PLANE_STRESS_YY,
                                        INDEX_2D_PLANE_STRESS_XY };

    // stress/strain vector indices (2D axisymmetric):
    enum indexStress2DAxisymmetric: int{ INDEX_2D_AXI_SYMMETRIC_XX,
                                         INDEX_2D_AXI_SYMMETRIC_YY,
                                         INDEX_2D_AXI_SYMMETRIC_RR,
                                         INDEX_2D_AXI_SYMMETRIC_XY};

    // stress/strain vector indices (2D interface):
    enum indexStress2DInterface: int{ INDEX_2D_INTERFACE_XZ,
                                      INDEX_2D_INTERFACE_ZZ };

    // stress/strain vector indices (3D interface):
    enum indexStress3DInterface: int{ INDEX_3D_INTERFACE_XZ,
                                      INDEX_3D_INTERFACE_YZ,
                                      INDEX_3D_INTERFACE_ZZ };

    // stress/strain vector indices 2D beam:
    enum indexStress2DBeam: int{ INDEX_2D_BEAM_XX,
                                 INDEX_2D_BEAM_YY,
                                 INDEX_2D_BEAM_XY };
    // DOF indices 2D beam:
    enum indexDOF2DBeam: int{ INDEX_2D_BEAM_X,
                              INDEX_2D_BEAM_Y,
                              INDEX_2D_BEAM_T };

    // Heat vector indices:
    enum class indexThermalFlux : int {
        X, Y, Z
    };

}