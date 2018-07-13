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
//                   Riccardo Rossi
//


// System includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"

// External includes


// Project includes
#include "custom_elements/meshless_base_element.h"
#include "custom_elements/truss_discrete_element.h"

// Application includes
#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    void TrussDiscreteElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        // definition of problem size
        const unsigned int number_of_control_points = GetGeometry().size();
        unsigned int number_of_dofs = number_of_control_points * 3;

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != number_of_dofs)
                rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_dofs, number_of_dofs); //resetting LHS
        }
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (rRightHandSideVector.size() != number_of_dofs)
                rRightHandSideVector.resize(number_of_dofs);
            rRightHandSideVector = ZeroVector(number_of_dofs); //resetting RHS
        }

        //reading in of integration weight, shape function values and shape function derivatives
        const double& integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        const Vector&   N     = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix&  DN_De  = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        //Bending stabilization
        const double E = GetProperties()[YOUNG_MODULUS];
        const double area = GetProperties()[CROSS_AREA];
        const double prestress = GetProperties()[PRESTRESS_CAUCHY];

        Vector base_vector = ZeroVector(3);
        GetBaseVector(base_vector, DN_De);

        double A0 = norm_2(mBaseVector0);
        double a = norm_2(base_vector);

        //stresses
        //Green Lagrange formulation (strain)
        double E11_membrane = 0.5 * (pow(a, 2) - pow(A0,2)); 

        //normal force
        double S11_membrane = prestress * area + E11_membrane * area * E / pow(A0,2);

        // 1st variation of the axial strain 
        Vector epsilon_var_dof = ZeroVector(number_of_dofs);
        Get1stVariationsAxialStrain(epsilon_var_dof, base_vector, 3, DN_De);
        epsilon_var_dof = epsilon_var_dof / pow(A0, 2);

        // 2nd variation of the axial strain 
        Matrix epsilon_var_2_dof = ZeroMatrix(number_of_dofs, number_of_dofs);
        Get2ndVariationsAxialStrain(epsilon_var_2_dof, 3, DN_De);
        epsilon_var_2_dof = epsilon_var_2_dof / pow(A0, 2);

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            for (int r = 0; r < number_of_dofs; r++) {
                for (int s = 0; s < number_of_dofs; s++) {
                    rLeftHandSideMatrix(r, s) = E * area * epsilon_var_dof[r] * epsilon_var_dof[s] + S11_membrane * epsilon_var_2_dof(r, s);
                }
            }
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            rRightHandSideVector = -S11_membrane * epsilon_var_dof;
        }

        rLeftHandSideMatrix = rLeftHandSideMatrix * integration_weight;
        rRightHandSideVector = rRightHandSideVector * integration_weight;

        KRATOS_CATCH("");
    }
} // Namespace Kratos


