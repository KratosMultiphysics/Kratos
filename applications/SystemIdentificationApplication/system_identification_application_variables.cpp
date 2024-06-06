//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Ihar Antonau
//                   Fabian Meister
//

// Include base h
#include "system_identification_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_DISPLACEMENT)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_ROTATION)

    KRATOS_CREATE_VARIABLE(double, PERTURBATION_SIZE)
    KRATOS_CREATE_VARIABLE(bool, ADAPT_PERTURBATION_SIZE)
    KRATOS_CREATE_VARIABLE(bool, HAS_ROTATION_DOFS)

    KRATOS_CREATE_VARIABLE(std::string, SENSOR_NAME)
    KRATOS_CREATE_VARIABLE(std::string, TEST_ANALYSIS_NAME)
    KRATOS_CREATE_VARIABLE(double, SENSOR_MEASURED_VALUE)
    KRATOS_CREATE_VARIABLE(Vector, SENSOR_MEASURED_VALUE_VECTOR)
    KRATOS_CREATE_VARIABLE(double, SENSOR_ERROR)
    KRATOS_CREATE_VARIABLE(int, SENSOR_ELEMENT_ID)

    KRATOS_CREATE_VARIABLE(Vector, EIGENVALUE_VECTOR)
    KRATOS_CREATE_VARIABLE(Matrix, EIGENVECTOR_MATRIX)

    KRATOS_CREATE_VARIABLE( std::string, DESIGN_VARIABLE_NAME)
    KRATOS_CREATE_VARIABLE( Matrix, LHS_DESIGN_DERIVATIVE )
    KRATOS_CREATE_VARIABLE( Matrix, MASS_DESIGN_DERIVATIVE )
    KRATOS_CREATE_VARIABLE( Matrix, LHS_PRIMAL_ELEMENT )
    KRATOS_CREATE_VARIABLE( Matrix, MASS_MATRIX_PRIMAL_ELEMENT )


} // namespace Kratos
