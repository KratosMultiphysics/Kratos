//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/curve_base_discrete_condition.h"

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"

namespace Kratos
{
    void CurveBaseDiscreteCondition::Initialize()
    {
        KRATOS_TRY

        Vector base_vector = ZeroVector(3);
        CalculateBaseVector(base_vector, GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES));
        mBaseVector0 = base_vector;

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************
    void CurveBaseDiscreteCondition::CalculateBaseVector(
        Vector& rBaseVector,
        const Matrix& rDN_De
    )
    {
        if (rBaseVector.size() != 3)
            rBaseVector.resize(3);
        rBaseVector = ZeroVector(3);

        // this is valid for all parameter edges
        if (Has(TANGENTS))
        {
            GetBoundaryEdgeBaseVector(rDN_De, GetValue(TANGENTS), rBaseVector);
        }
        else
        {
            int number_of_control_points = GetGeometry().size();

            for (int i = 0; i < number_of_control_points; i++)
            {
                rBaseVector[0] += rDN_De(0, i) * GetGeometry()[i].X();
                rBaseVector[1] += rDN_De(0, i) * GetGeometry()[i].Y();
                rBaseVector[2] += rDN_De(0, i) * GetGeometry()[i].Z();
            }
        }
    }

    //************************************************************************************
    //************************************************************************************
    void CurveBaseDiscreteCondition::CalculateNormalVector(
        Vector& rNormalVector,
        const Matrix& rDN_De
    )
    {
        if (rNormalVector.size() != 3)
            rNormalVector.resize(3);
        rNormalVector = ZeroVector(3);

        Matrix J;
        CalculateJacobian(rDN_De, J);

        //basis vectors g1 and g2
        array_1d<double, 3> g1;
        array_1d<double, 3> g2;

        g1[0] = J(0, 0);
        g2[0] = J(0, 1);
        g1[1] = J(1, 0);
        g2[1] = J(1, 1);
        g1[2] = J(2, 0);
        g2[2] = J(2, 1);

        MathUtils<double>::CrossProduct(rNormalVector, g1, g2);
    }

    //************************************************************************************
    //************************************************************************************
    void CurveBaseDiscreteCondition::GetBoundaryEdgeBaseVector(
        const Matrix& DN_De,
        const array_1d<double, 2>& Tangents,
        Vector& rBaseVector)
    {
        Matrix J;
        CalculateJacobian(DN_De, J);

        //basis vectors g1 and g2
        array_1d<double, 3> g1;
        array_1d<double, 3> g2;

        g1[0] = J(0, 0);
        g2[0] = J(0, 1);
        g1[1] = J(1, 0);
        g2[1] = J(1, 1);
        g1[2] = J(2, 0);
        g2[2] = J(2, 1);

        rBaseVector = g1 * Tangents[0] + g2 * Tangents[1];
    }
} // Namespace Kratos


