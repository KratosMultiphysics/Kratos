//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "conservative_condition.h"
#include "includes/checks.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"

namespace Kratos
{

template<std::size_t TNumNodes>
int ConservativeCondition<TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int err = Condition::Check(rCurrentProcessInfo);
    if (err != 0) return err;

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry())
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, r_node)

        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, r_node)
    }

    return err;

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void ConservativeCondition<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rResult.size() != mLocalSize)
        rResult.resize(mLocalSize, false); // False says not to preserve existing storage!!

    const GeometryType& r_geom = this->GetGeometry();
    int counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = r_geom[i].GetDof(MOMENTUM_X).EquationId();
        rResult[counter++] = r_geom[i].GetDof(MOMENTUM_Y).EquationId();
        rResult[counter++] = r_geom[i].GetDof(HEIGHT).EquationId();
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void ConservativeCondition<TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rConditionDofList.size() != mLocalSize)
        rConditionDofList.resize(mLocalSize);

    const GeometryType& r_geom = this->GetGeometry();
    int counter=0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rConditionDofList[counter++] = r_geom[i].pGetDof(MOMENTUM_X);
        rConditionDofList[counter++] = r_geom[i].pGetDof(MOMENTUM_Y);
        rConditionDofList[counter++] = r_geom[i].pGetDof(HEIGHT);
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void ConservativeCondition<TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(MOMENTUM_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(MOMENTUM_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(HEIGHT, Step);
    }
}

template<std::size_t TNumNodes>
void ConservativeCondition<TNumNodes>::AddWaveTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const double Weight)
{
    const bool integrate_by_parts = true;
    const auto z = rData.topography;
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];
    const double c2 = rData.gravity * rData.height;
    const auto n = rData.normal;

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const IndexType i_block = 3 * i;
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const IndexType j_block = 3 * j;

            double n_ij;
            if (integrate_by_parts) {
                n_ij = rN[i] * rN[j];
            } else {
                n_ij = 0.0;
            }

            /* First component
             * A_1 = {{2 * u_1   0   -u_1^2 + gh},
             *        {  u_2    u_1   -u_1 * u_2},
             *        {   1      0        0     }}
             * b_1 = {gh, 0, 0}^T
            */
            rMatrix(i_block,     j_block)     -= Weight * n_ij * n[0] * 2*u_1;
            rMatrix(i_block,     j_block + 2) -= Weight * n_ij * n[0] * (-u_1*u_1 + c2);
            rMatrix(i_block + 1, j_block)     -= Weight * n_ij * n[0] * u_2;
            rMatrix(i_block + 1, j_block + 1) -= Weight * n_ij * n[0] * u_1;
            rMatrix(i_block + 1, j_block + 2) -= Weight * n_ij * n[0] * (-u_1*u_2);
            rMatrix(i_block + 2, j_block)     -= Weight * n_ij * n[0];
            rVector[i_block]                  += Weight * n_ij * n[0] * c2 * z[j];

            /* Second component
             * A_2 = {{u_2    u_1      -u_1 * u_2},
             *        { 0   2 * u_2   -u_2^2 + gh},
             *        { 0      1            0    }}
             * b_2 = {0, gh, 0}^T
            */
            rMatrix(i_block,     j_block)     -= Weight * n_ij * n[1] * u_2;
            rMatrix(i_block,     j_block + 1) -= Weight * n_ij * n[1] * u_1;
            rMatrix(i_block,     j_block + 2) -= Weight * n_ij * n[1] * (-u_1*u_2);
            rMatrix(i_block + 1, j_block + 1) -= Weight * n_ij * n[1] * 2*u_1;
            rMatrix(i_block + 1, j_block + 2) -= Weight * n_ij * n[1] * (-u_1*u_1 + c2);
            rMatrix(i_block + 2, j_block + 1) -= Weight * n_ij * n[1];
            rVector[i_block + 1]              += Weight * n_ij * n[1] * c2 * z[j];
        }
    }
}

template<std::size_t TNumNodes>
void ConservativeCondition<TNumNodes>::AddFluxTerms(
    LocalVectorType& rVector,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const double Weight)
{
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];
    const double c2 = rData.gravity * rData.height;
    const double l = this->StabilizationParameter(rData);

}

template<std::size_t TNumNodes>
double ConservativeCondition<TNumNodes>::StabilizationParameter(const ConditionData& rData) const
{
    const double lambda = std::sqrt(rData.gravity * rData.height) + norm_2(rData.velocity);
    const double epsilon = 1e-6;
    const double threshold = rData.relative_dry_height * rData.length;
    const double w = ShallowWaterUtilities().WetFraction(rData.height, threshold);
    return w * rData.length * rData.stab_factor / (lambda + epsilon);
}

template class ConservativeCondition<2>;

} // namespace Kratos
