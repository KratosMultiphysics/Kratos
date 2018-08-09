//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "custom_conditions/point_base_discrete_condition.h"

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"

namespace Kratos
{
    void PointBaseDiscreteCondition::Initialize()
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }


    //************************************************************************************
    //************************************************************************************
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

    //***********************************************************************************
    //***********************************************************************************
    void PointBaseDiscreteCondition::CalculateBaseVector(
        Vector& rBaseVector,
        const Matrix& rDN_De)
    {
        if (rBaseVector.size() != 3)
            rBaseVector.resize(3);
        rBaseVector = ZeroVector(3);

        Matrix Jacobian;
        CalculateJacobian(rDN_De, Jacobian, 3, 2);

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


