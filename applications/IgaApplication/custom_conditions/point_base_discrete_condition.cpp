//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// System includes

// External includes

// Project includes
#include "point_base_discrete_condition.h"


namespace Kratos
{
    void PointBaseDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_ERROR << "You have called to the CalculateAll from the base class for discrete point conditions" << std::endl;
    }

    void PointBaseDiscreteCondition::CalculateBaseVectorOnSurface(
        Vector& rBaseVector)
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        KRATOS_ERROR_IF(DN_De.size2() != 2) << "PointBaseDiscreteCondition::CalculateBaseVectorOnSurface: Cannot compute base vector on surface, because dimension of underlying topology is not size 2." << std::endl;

        if (rBaseVector.size() != 3)
            rBaseVector.resize(3);
        rBaseVector = ZeroVector(3);

        Matrix Jacobian;
        CalculateJacobian(DN_De, Jacobian, 3, 2);

        Vector g1 = ZeroVector(3);
        Vector g2 = ZeroVector(3);

        g1[0] = Jacobian(0, 0);
        g2[0] = Jacobian(0, 1);
        g1[1] = Jacobian(1, 0);
        g2[1] = Jacobian(1, 1);
        g1[2] = Jacobian(2, 0);
        g2[2] = Jacobian(2, 1);

        MathUtils<double>::CrossProduct(rBaseVector, g1, g2);
    }
} // Namespace Kratos


