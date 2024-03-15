//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "custom_conditions/grid_based_conditions/mpm_grid_navier_slip_friction_condition.h"
#include "includes/checks.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::InitializeSolutionStep(
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;

        bool is_active = true;
        const GeometryType& r_geometry = this->GetGeometry();

        for (const auto& r_node : r_geometry) {
            if (r_node.IsNot(ACTIVE)) {
                is_active = false;
                break;
            }
        }

        this->Set(ACTIVE, is_active);

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        KRATOS_TRY;

        if (rResult.size() != LocalSize) {
            rResult.resize(LocalSize);
        }

        const GeometryType& r_geometry = GetGeometry();
        const unsigned int pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

        unsigned int index = 0;

        for (const auto& r_node : r_geometry) {
            if (r_node.Is(ACTIVE)) {
                rResult[index++] = r_node.GetDof(DISPLACEMENT_X,pos  ).EquationId();
                rResult[index++] = r_node.GetDof(DISPLACEMENT_Y,pos+1).EquationId();
                if constexpr (TDim == 3) {
                    rResult[index++] = r_node.GetDof(DISPLACEMENT_Z,pos+2).EquationId();
                }
            } else {
                for (unsigned int i=0; i<TDim; ++i) rResult[i] = 0;
            }
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::GetDofList(
        DofsVectorType& rConditionDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        KRATOS_TRY;

        rConditionDofList.resize(0);
        rConditionDofList.reserve(LocalSize);

        const GeometryType& r_geometry = GetGeometry();
        const unsigned int pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

        for (const auto& r_node : r_geometry) {
            if (r_node.Is(ACTIVE)) {
                rConditionDofList.push_back(r_node.pGetDof(DISPLACEMENT_X,pos));
                rConditionDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y,pos+1));
                if constexpr (TDim == 3) {
                    rConditionDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z,pos+2));
                }
            }
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::GetFirstDerivativesVector(
        Vector& rValues,
        int Step
        ) const
    {
        KRATOS_TRY;

        if (rValues.size() != LocalSize) {
            rValues.resize(LocalSize, false);
        }

        noalias(rValues) = ZeroVector(LocalSize);

        const GeometryType& r_geometry = GetGeometry();

        for (unsigned int i = 0; i < TNumNodes; i++) {
            const array_1d<double,3> & r_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * TDim;
            for(unsigned int k = 0; k<TDim; ++k) {
                rValues[index + k] = r_velocity[k];
            }
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        )
    {

        if ( rDampingMatrix.size1() != LocalSize ) {
            rDampingMatrix.resize(LocalSize, LocalSize, false);
        }

        noalias(rDampingMatrix) = ZeroMatrix(LocalSize,LocalSize);

        const auto& r_geom = this->GetGeometry();

        // Compute condition unit normal vector
        // We are assuming a constant normal over the entire condition
        array_1d<double,3> Normal = r_geom.UnitNormal(0, GeometryData::IntegrationMethod::GI_GAUSS_1);

        // Set the tangential projection matrix
        BoundedMatrix<double,TDim,TDim> tang_proj_mat;
        noalias(tang_proj_mat) = IdentityMatrix(TDim,TDim);
        for (IndexType d1 = 0; d1 < TDim; ++d1) {
            for (IndexType d2 = 0; d2 < TDim; ++d2) {
                tang_proj_mat(d1,d2) -= Normal[d1]*Normal[d2];
            }
        }

        // Calculate the required Gauss points integration data
        const IntegrationMethod& r_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& r_integration_pts = r_geom.IntegrationPoints(r_integration_method);
        const SizeType n_gauss = r_integration_pts.size();
        Vector GaussPtsWeights;
        r_geom.DeterminantOfJacobian(GaussPtsWeights, r_integration_method);
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            GaussPtsWeights[i_gauss] = GaussPtsWeights[i_gauss] * r_integration_pts[i_gauss].Weight();
        }
        const Matrix ShapeFunctionsContainer = r_geom.ShapeFunctionsValues(r_integration_method);

        const double& friction_coefficient = this->GetValue(FRICTION_COEFFICIENT);

        // Calculate the Navier-slip contribution
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Get current Gauss point data
            const double w_gauss = GaussPtsWeights[i_gauss];
            const auto& N_gauss = row(ShapeFunctionsContainer, i_gauss);
            const double aux_val = w_gauss * friction_coefficient;
            // Assemble RHS and LHS contributions
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                    for (IndexType d1 = 0; d1 < TDim; ++d1) {
                        for (IndexType d2 = 0; d2 < TDim; ++d2) {
                            rDampingMatrix(i_node*BlockSize + d1, j_node*BlockSize + d2) += aux_val * N_gauss[i_node] * N_gauss[j_node] * tang_proj_mat(d1,d2);
                        }
                    }
                }
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = false;
        VectorType temp = Vector();

        CalculateAll(
            rLeftHandSideMatrix,
            temp,
            rCurrentProcessInfo,
            CalculateStiffnessMatrixFlag,
            CalculateResidualVectorFlag);
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll(
            temp,
            rRightHandSideVector,
            rCurrentProcessInfo,
            CalculateStiffnessMatrixFlag,
            CalculateResidualVectorFlag);
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll(
            rLeftHandSideMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            CalculateStiffnessMatrixFlag,
            CalculateResidualVectorFlag);
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    void MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag
        )
    {
        if (CalculateStiffnessMatrixFlag && rLeftHandSideMatrix.size1() != LocalSize) {
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
        }
        if (CalculateResidualVectorFlag && rRightHandSideVector.size() != LocalSize) {
            rRightHandSideVector.resize(LocalSize, false);
        }
    }

    //************************************************************************************
    //************************************************************************************

    template<unsigned int TDim, unsigned int TNumNodes>
    int MPMGridNavierSlipFrictionCondition<TDim,TNumNodes>::Check(
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        KRATOS_TRY;

        // Base check
        int check = Condition::Check(rCurrentProcessInfo);

        if (check!=0) {
            return check;
        } else {
            // Check if friction coefficient is defined
            KRATOS_CHECK(this->Has(FRICTION_COEFFICIENT));
            KRATOS_CHECK_GREATER(this->GetValue(FRICTION_COEFFICIENT),0);
            const GeometryType& r_geometry = this->GetGeometry();
            for (const auto& r_node : r_geometry) {
                // Check nodal variables
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL, r_node);
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(IS_STRUCTURE, r_node);
                // Check DOFs
                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node);
            }
            return 0;
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    template class MPMGridNavierSlipFrictionCondition<2,2>;
    template class MPMGridNavierSlipFrictionCondition<3,3>;

} // Namespace Kratos
