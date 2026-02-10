#include "compressible_sbm_slip_condition.h"

namespace Kratos
{


    // EquationIdVector

    template <>
    void CompressibleSBMSlipCondition<2, 2>::EquationIdVector(
        EquationIdVectorType &rResult,
        const ProcessInfo &) const
    {
        KRATOS_WATCH("SBM::EquationIdVector<2,2>");
        rResult.resize(4, false);

        rResult[0] = GetGeometry()[0].GetDof(VELOCITY_X).EquationId();
        rResult[1] = GetGeometry()[0].GetDof(VELOCITY_Y).EquationId();
        rResult[2] = GetGeometry()[1].GetDof(VELOCITY_X).EquationId();
        rResult[3] = GetGeometry()[1].GetDof(VELOCITY_Y).EquationId();
    }

    template <>
    void CompressibleSBMSlipCondition<3, 3>::EquationIdVector(
        EquationIdVectorType &rResult,
        const ProcessInfo &) const
    {
        KRATOS_WATCH("SBM::EquationIdVector<3,3>");
        rResult.resize(9, false);

        unsigned int k = 0;
        for (unsigned int i = 0; i < 3; ++i)
        {
            rResult[k++] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
            rResult[k++] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
            rResult[k++] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
        }
    }

    // GetDofList 

    template <>
    void CompressibleSBMSlipCondition<2, 2>::GetDofList(
        DofsVectorType &rDofs,
        const ProcessInfo &) const
    {
        KRATOS_WATCH("SBM::GetDofList<2,2>");
        rDofs.resize(4);

        rDofs[0] = GetGeometry()[0].pGetDof(VELOCITY_X);
        rDofs[1] = GetGeometry()[0].pGetDof(VELOCITY_Y);
        rDofs[2] = GetGeometry()[1].pGetDof(VELOCITY_X);
        rDofs[3] = GetGeometry()[1].pGetDof(VELOCITY_Y);
    }

    template <>
    void CompressibleSBMSlipCondition<3, 3>::GetDofList(
        DofsVectorType &rDofs,
        const ProcessInfo &) const
    {
        KRATOS_WATCH("SBM::GetDofList<3,3>");
        rDofs.resize(9);

        unsigned int k = 0;
        for (unsigned int i = 0; i < 3; ++i)
        {
            rDofs[k++] = GetGeometry()[i].pGetDof(VELOCITY_X);
            rDofs[k++] = GetGeometry()[i].pGetDof(VELOCITY_Y);
            rDofs[k++] = GetGeometry()[i].pGetDof(VELOCITY_Z);
        }
    }


    // CalculateRightHandSide — PROOF OF EXECUTION 

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleSBMSlipCondition<TDim, TNumNodes>::CalculateRightHandSide(
        VectorType &rRHS,
        const ProcessInfo &)
    {
        KRATOS_WATCH("SBM CONDITION RHS CALLED ");
        rRHS.resize(TDim * TNumNodes, false);
        rRHS.clear();
    }


    /* Explicit instantiation */

    template class CompressibleSBMSlipCondition<2, 2>;
    template class CompressibleSBMSlipCondition<3, 3>;

} // namespace Kratos
