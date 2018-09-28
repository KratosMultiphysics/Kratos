/*
//  KRATOS .___  ________    _____
//         |   |/  _____/   /  _  \
//         |   /   \  ___  /  /_\  \
//         |   \    \_\  \/    |    \
//         |___|\______  /\____|__  /
//                     \/         \/  Application
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//
//  Authors: Tobias Teschemacher
*/

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/curve_base_discrete_element.h"

#include "iga_application.h"
#include "iga_application_variables.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"

namespace Kratos
{
    void CurveBaseDiscreteElement::Initialize()
    {
        KRATOS_TRY

        //Constitutive Law initialisation
        BaseDiscreteElement::InitializeMaterial();

        Vector base_vector = ZeroVector(3);
        GetBaseVector(base_vector, GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES));
        mBaseVector0 = base_vector;

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************
    void CurveBaseDiscreteElement::GetBaseVector(
        Vector& rBaseVector,
        const Matrix& rDN_De
    )
    {
        if (rBaseVector.size() != 3)
            rBaseVector.resize(3);
        noalias(rBaseVector) = ZeroVector(3);

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
    void CurveBaseDiscreteElement::GetBoundaryEdgeBaseVector(
        const Matrix& DN_De,
        const array_1d<double, 2>& Tangents,
        Vector& rBaseVector)
    {
        Matrix J;
        Jacobian(DN_De, J);

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

    //************************************************************************************
    //************************************************************************************
    void CurveBaseDiscreteElement::Get1stVariationsAxialStrain(
        Vector& rEpsilon1stVariationDoF,
        const Vector& rBaseVector,
        const int& rNumberOfDoFs,
        const Matrix& rDN_De)
    {
        int mat_size = rDN_De.size1()*rNumberOfDoFs;

        if (rEpsilon1stVariationDoF.size() != mat_size)
            rEpsilon1stVariationDoF.resize(mat_size, false);
        rEpsilon1stVariationDoF = ZeroVector(mat_size);

        for (std::size_t r = 0; r < mat_size; r++)
        {
            int xyz_r = r % rNumberOfDoFs; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
            int i = r / rNumberOfDoFs;     // index for the shape functions
            if (xyz_r>2)
                rEpsilon1stVariationDoF[r] = 0.0;
            else
                rEpsilon1stVariationDoF[r] = rBaseVector[xyz_r] * rDN_De(i, 0);
        }
    }

    //************************************************************************************
    //************************************************************************************
    void CurveBaseDiscreteElement::Get2ndVariationsAxialStrain(
        Matrix& rEpsilon2ndVariationDoF,
        const int& rNumberOfDoFs,
        const Matrix& rDN_De)
    {
        int mat_size = rDN_De.size1()*rNumberOfDoFs;

        if ((rEpsilon2ndVariationDoF.size1() != mat_size) || (rEpsilon2ndVariationDoF.size2() != mat_size))
            rEpsilon2ndVariationDoF.resize(mat_size, mat_size, false);
        rEpsilon2ndVariationDoF = ZeroMatrix(mat_size, mat_size);

        for (std::size_t r = 0; r<mat_size; r++) //in the case
        {
            int xyz_r = r % rNumberOfDoFs; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
            int i = r / rNumberOfDoFs;     // index for the shape functions
            if (xyz_r>2)
                for (std::size_t s = 0; s<mat_size; s++)
                    rEpsilon2ndVariationDoF(r, s) = 0.0;
            else
            {
                for (std::size_t s = 0; s<mat_size; s++)
                {
                    int xyz_s = s % rNumberOfDoFs; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int j = s / rNumberOfDoFs;     // index for the shape functions
                    if (xyz_s>2)
                        rEpsilon2ndVariationDoF(r, s) = 0;
                    else
                        if (xyz_r == xyz_s)
                            rEpsilon2ndVariationDoF(r, s) = rDN_De(i, 0) * rDN_De(j, 0);
                        else
                            rEpsilon2ndVariationDoF(r, s) = 0;
                }
            }
        }

    }
} // Namespace Kratos


