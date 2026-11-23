//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2026-10-27 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/monolithic_2fluid_2d_aux.h"
#include "pfem_2_application_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

namespace
{
    template <class TMatrixType>
    void FillConsistentMassMatrix(TMatrixType& rMassMatrix, const double Area)
    {
        const double factor = Area / 12.0;
        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)
            {
                rMassMatrix(i, j) = factor;
            }
            rMassMatrix(i, i) += factor;
        }
    }
}

    //************************************************************************************
    //************************************************************************************
    MonolithicPFEM22DAux::MonolithicPFEM22DAux(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {}

    //************************************************************************************
    //************************************************************************************
    MonolithicPFEM22DAux::MonolithicPFEM22DAux(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {}

    Element::Pointer MonolithicPFEM22DAux::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new MonolithicPFEM22DAux(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer MonolithicPFEM22DAux::Create(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new MonolithicPFEM22DAux(NewId, pGeometry, pProperties));
    }

    MonolithicPFEM22DAux::~MonolithicPFEM22DAux()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DAux::InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo)
    {
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DAux::GetNodalSource(array_1d<double, 3>& rSource, const IndexType StepIndex) const
    {
        const GeometryType& rGeom = this->GetGeometry();
        rSource = ZeroVector(3);

        bool found_nodal = false;
        for (unsigned int i = 0; i < 3; ++i)
        {
            const IndexType step = (rGeom[i].GetBufferSize() > StepIndex) ? StepIndex : 0;
            if (rGeom[i].SolutionStepsDataHas(HEAT_FLUX))
            {
                rSource[i] = rGeom[i].FastGetSolutionStepValue(HEAT_FLUX, step);
                found_nodal = true;
            }
            else if (rGeom[i].SolutionStepsDataHas(BODY_FORCE))
            {
                rSource[i] = rGeom[i].FastGetSolutionStepValue(BODY_FORCE, step)[0];
                found_nodal = true;
            }
        }

        // If no nodal source variable is present, fall back to property-level HEAT_FLUX
        if (!found_nodal)
        {
            double prop_source = 0.0;
            if (this->GetProperties().Has(HEAT_FLUX))
                prop_source = this->GetProperties()[HEAT_FLUX];

            for (unsigned int i = 0; i < 3; ++i)
                rSource[i] = prop_source;
        }
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DAux::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int NumNodes = 3;
        if (rMassMatrix.size1() != NumNodes || rMassMatrix.size2() != NumNodes)
            rMassMatrix.resize(NumNodes, NumNodes, false);
        noalias(rMassMatrix) = ZeroMatrix(NumNodes, NumNodes);

        double Area;
        BoundedMatrix<double, 3, 2> DN_DX;
        array_1d<double, 3> N;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Consistent mass matrix for linear 2D triangle:
        // M_ij = (Area / 12) * (1 + delta_ij)
        FillConsistentMassMatrix(rMassMatrix, Area);

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DAux::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int NumNodes = 3;

        if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
            rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes, NumNodes);

        if (rRightHandSideVector.size() != NumNodes)
            rRightHandSideVector.resize(NumNodes, false);
        noalias(rRightHandSideVector) = ZeroVector(NumNodes);

        KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(DELTA_TIME))
            << "DELTA_TIME is not set in ProcessInfo." << std::endl;
        const double dt = rCurrentProcessInfo[DELTA_TIME];
        KRATOS_ERROR_IF(dt <= 0.0)
            << "DELTA_TIME must be positive, got " << dt << std::endl;
        const double inv_dt = 1.0 / dt;

        // theta = 1.0 -> BDF1 (Backward Euler), 0.5 -> Crank-Nicolson, 0.0 -> Forward Euler
        const double theta = rCurrentProcessInfo.Has(TIMEDIS_THETA)
            ? rCurrentProcessInfo[TIMEDIS_THETA]
            : 1.0;

        const double diffusivity = this->GetProperties().Has(CONDUCTIVITY)
            ? this->GetProperties()[CONDUCTIVITY]
            : 1.0;

        double Area;
        BoundedMatrix<double, 3, 2> DN_DX;
        array_1d<double, 3> N;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Consistent mass matrix M_ij = (Area / 12) * (1 + delta_ij)
        BoundedMatrix<double, 3, 3> M;
        FillConsistentMassMatrix(M, Area);

        // Diffusion (stiffness) matrix K_ij = Area * k * (grad N_i . grad N_j)
        const BoundedMatrix<double, 3, 3> K = Area * diffusivity * prod(DN_DX, trans(DN_DX));

        array_1d<double, 3> T_curr, T_prev;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            const auto& rNode = this->GetGeometry()[i];
            T_curr[i] = rNode.FastGetSolutionStepValue(TEMPERATURE, 0);
            T_prev[i] = (rNode.GetBufferSize() > 1)
                ? rNode.FastGetSolutionStepValue(TEMPERATURE, 1)
                : T_curr[i];
        }

        array_1d<double, 3> f_curr, f_prev;
        this->GetNodalSource(f_curr, 0);
        this->GetNodalSource(f_prev, 1);

        const array_1d<double, 3> f_theta = theta * f_curr + (1.0 - theta) * f_prev;
        const array_1d<double, 3> T_theta = theta * T_curr + (1.0 - theta) * T_prev;

        // Theta-method residual form: (1/dt * M + theta * K) * dT = M * f_theta - (1/dt) * M * (T_curr - T_prev) - K * T_theta
        noalias(rLeftHandSideMatrix) = inv_dt * M + theta * K;
        noalias(rRightHandSideVector) = prod(M, f_theta)
            - inv_dt * prod(M, T_curr - T_prev)
            - prod(K, T_theta);

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DAux::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        MatrixType LHS;
        this->CalculateLocalSystem(LHS, rRightHandSideVector, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DAux::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
    {
        const SizeType NumNodes = 3;
        const GeometryType& rGeom = this->GetGeometry();

        if (rResult.size() != NumNodes)
            rResult.resize(NumNodes, false);

        for (SizeType i = 0; i < NumNodes; ++i)
        {
            rResult[i] = rGeom[i].GetDof(TEMPERATURE).EquationId();
        }
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DAux::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const
    {
        const SizeType NumNodes = 3;
        const GeometryType& rGeom = this->GetGeometry();

        if (ElementalDofList.size() != NumNodes)
            ElementalDofList.resize(NumNodes);

        for (SizeType i = 0; i < NumNodes; ++i)
        {
            ElementalDofList[i] = rGeom[i].pGetDof(TEMPERATURE);
        }
    }

} // Namespace Kratos
