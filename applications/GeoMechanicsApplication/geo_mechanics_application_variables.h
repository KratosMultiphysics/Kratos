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

#if !defined(KRATOS_GEO_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_GEO_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geo_mechanics_application_constants.h"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include <string>

namespace Kratos
{
/* Variables from GeoMechanicsApplication */
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, VELOCITY_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, DT_PRESSURE_COEFFICIENT)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, DT_WATER_PRESSURE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NORMAL_FLUID_FLUX)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, HYDRAULIC_HEAD)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, HYDRAULIC_DISCHARGE)

// K0 values
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, int, K0_MAIN_DIRECTION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, K0_VALUE_XX)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, K0_VALUE_YY)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, K0_VALUE_ZZ)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, K0_NC)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, OCR)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, POISSON_UNLOADING_RELOADING)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, POP)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, DENSITY_SOLID)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, BULK_MODULUS_SOLID)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, BULK_MODULUS_FLUID)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PERMEABILITY_XX)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PERMEABILITY_YY)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PERMEABILITY_ZZ)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PERMEABILITY_XY)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PERMEABILITY_YZ)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PERMEABILITY_ZX)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PERMEABILITY_CHANGE_INVERSE_FACTOR)

// Mohr-Coulomb
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, GEO_COHESION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, GEO_FRICTION_ANGLE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, GEO_DILATANCY_ANGLE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, GEO_TENSILE_STRENGTH)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, SPECIFIC_HEAT_CAPACITY_WATER)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, SPECIFIC_HEAT_CAPACITY_SOLID)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, THERMAL_CONDUCTIVITY_WATER)
KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, THERMAL_CONDUCTIVITY_SOLID)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, SOLID_COMPRESSIBILITY)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, DT_TEMPERATURE_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, DT_TEMPERATURE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NORMAL_HEAT_FLUX)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, std::string, THERMAL_LAW_NAME)

// Variables for Micro-Climate boundary
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, AIR_TEMPERATURE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, SOLAR_RADIATION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, AIR_HUMIDITY)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PRECIPITATION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, WIND_SPEED)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, A1_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, A2_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, A3_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, ALPHA_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, QF_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, SMIN_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, SMAX_COEFFICIENT)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, MINIMUM_JOINT_WIDTH)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, TRANSVERSAL_PERMEABILITY)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, FLUID_FLUX_VECTOR)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, LOCAL_FLUID_FLUX_VECTOR)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, LOCAL_STRESS_VECTOR)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, LOCAL_RELATIVE_DISPLACEMENT_VECTOR)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Matrix, PERMEABILITY_MATRIX)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Matrix, LOCAL_PERMEABILITY_MATRIX)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, ACCUMULATED_STRAIN)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, CRITICAL_DISPLACEMENT)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, TOTAL_DISPLACEMENT)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, INCREMENTAL_DISPLACEMENT)

KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, TOTAL_ROTATION)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(GEO_MECHANICS_APPLICATION, INCREMENTAL_ROTATION)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, IS_CONVERGED)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Matrix, TOTAL_STRESS_TENSOR)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Vector, TOTAL_STRESS_VECTOR)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Matrix, CAUCHY_STRAIN_TENSOR)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Vector, CAUCHY_STRAIN_VECTOR)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, TIME_UNIT_CONVERTER)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, LOCAL_EQUIVALENT_STRAIN)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NONLOCAL_EQUIVALENT_STRAIN)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, JOINT_WIDTH)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, NODAL_SMOOTHING)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Matrix, NODAL_CAUCHY_STRESS_TENSOR)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Matrix, ENGINEERING_STRAIN_TENSOR)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Vector, ENGINEERING_STRAIN_VECTOR)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NODAL_DAMAGE_VARIABLE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NODAL_JOINT_AREA)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NODAL_JOINT_WIDTH)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NODAL_JOINT_DAMAGE)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, BIOT_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PLATE_SHAPE_CORRECTION_FACTOR)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, MEAN_EFFECTIVE_STRESS)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, MEAN_STRESS)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, ENGINEERING_VOLUMETRIC_STRAIN)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, ENGINEERING_VON_MISES_STRAIN)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, GREEN_LAGRANGE_VOLUMETRIC_STRAIN)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, GREEN_LAGRANGE_VON_MISES_STRAIN)

// Nodal load variables
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(STRUCTURAL_MECHANICS_APPLICATION, POINT_LOAD)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(STRUCTURAL_MECHANICS_APPLICATION, LINE_LOAD)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(STRUCTURAL_MECHANICS_APPLICATION, SURFACE_LOAD)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, RESET_DISPLACEMENTS)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, CONSIDER_GEOMETRIC_STIFFNESS)

// Gap/closure
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, CONSIDER_GAP_CLOSURE)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, USE_CONSISTENT_MASS_MATRIX)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, IGNORE_UNDRAINED)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, USE_HENCKY_STRAIN)

// Retention models
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, SATURATED_SATURATION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, RESIDUAL_SATURATION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, VAN_GENUCHTEN_AIR_ENTRY_PRESSURE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, VAN_GENUCHTEN_GN)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, VAN_GENUCHTEN_GL)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, MINIMUM_RELATIVE_PERMEABILITY)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, std::string, RETENTION_LAW)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, DEGREE_OF_SATURATION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, EFFECTIVE_SATURATION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, BISHOP_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, DERIVATIVE_OF_SATURATION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, RELATIVE_PERMEABILITY)

// absorbing boundary
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Vector, ABSORBING_FACTORS)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, VIRTUAL_THICKNESS)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, CONFINED_STIFFNESS)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, SHEAR_STIFFNESS)

// pipe elements
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, IS_PIPING_CONVERGED)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PIPE_ETA)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PIPE_THETA)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PIPE_D_70)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, int, PIPE_START_ELEMENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PIPE_ELEMENT_LENGTH)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, PIPE_IN_EQUILIBRIUM)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, PIPE_MODIFIED_D)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PIPE_MODEL_FACTOR)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PIPE_WIDTH_FACTOR)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PIPE_HEIGHT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, PREV_PIPE_HEIGHT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, DIFF_PIPE_HEIGHT)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, PIPE_EROSION)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, PIPE_ACTIVE)

// Filter
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, FILTER_LENGTH)

// UDSM
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, std::string, UDSM_NAME) // also for umat
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, int, UDSM_NUMBER)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, bool, IS_FORTRAN_UDSM) // also for umat

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Vector, UMAT_PARAMETERS)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, int, INDEX_OF_UMAT_C_PARAMETER)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, int, INDEX_OF_UMAT_PHI_PARAMETER)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Vector, STATE_VARIABLES)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_1)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_2)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_3)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_4)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_5)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_6)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_7)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_8)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_9)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_10)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_11)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_12)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_13)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_14)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_15)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_16)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_17)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_18)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_19)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_20)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_21)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_22)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_23)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_24)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_25)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_26)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_27)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_28)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_29)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_30)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_31)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_32)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_33)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_34)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_35)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_36)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_37)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_38)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_39)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_40)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_41)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_42)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_43)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_44)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_45)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_46)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_47)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_48)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_49)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, STATE_VARIABLE_50)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Vector, STRAINS_OF_PIECEWISE_LINEAR_LAW)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, Vector, STRESSES_OF_PIECEWISE_LINEAR_LAW)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, INTERFACE_NORMAL_STIFFNESS)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, INTERFACE_SHEAR_STIFFNESS)

KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, GEO_SHEAR_CAPACITY)

} // namespace Kratos

#endif /* KRATOS_GEO_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
