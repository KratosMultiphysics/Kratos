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

    void MonolithicPFEM22DAux::GetNodalSource(array_1d<double, 3>& rSource) const
    {
        const GeometryType& rGeom = this->GetGeometry();
        rSource = ZeroVector(3);

        bool found_nodal = false;
        for (unsigned int i = 0; i < 3; ++i)
        {
            if (rGeom[i].SolutionStepsDataHas(HEAT_FLUX))
            {
                rSource[i] = rGeom[i].FastGetSolutionStepValue(HEAT_FLUX);
                found_nodal = true;
            }
            else if (rGeom[i].SolutionStepsDataHas(BODY_FORCE))
            {
                rSource[i] = rGeom[i].FastGetSolutionStepValue(BODY_FORCE)[0];
                found_nodal = true;
            }
        }

        // If no non-zero nodal source was found, check property-level HEAT_FLUX / HEAT_SOURCE
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
        const double factor = Area / 12.0;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                rMassMatrix(i, j) = factor;
            }
            rMassMatrix(i, i) += factor;
        }

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void MonolithicPFEM22DAux::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int NumNodes = 3;

        // Check sizes
        if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
            rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes, NumNodes);

        if (rRightHandSideVector.size() != NumNodes)
            rRightHandSideVector.resize(NumNodes, false);
        noalias(rRightHandSideVector) = ZeroVector(NumNodes);

        // Compute consistent mass matrix M_ij = \int N_i N_j d\Omega
        this->CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);

        // Get nodal source values f_j
        array_1d<double, 3> f_nodal;
        this->GetNodalSource(f_nodal);

        // (w, f)_i = \int N_i f(\mathbf{x}) d\Omega = \sum_j M_ij * f_j
        noalias(rRightHandSideVector) = prod(rLeftHandSideMatrix, f_nodal);

        // Subtract current solution contribution to obtain standard residual form:
        // Residual = (w, f) - M * T_current
        array_1d<double, 3> current_temperature;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            current_temperature[i] = this->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
        }
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, current_temperature);

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
