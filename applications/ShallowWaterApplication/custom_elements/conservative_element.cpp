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
#include "includes/checks.h"
#include "conservative_element.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"

namespace Kratos
{

template<std::size_t TNumNodes>
int ConservativeElement<TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int err = Element::Check(rCurrentProcessInfo);
    if (err != 0) return err;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry())
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MANNING, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VERTICAL_VELOCITY, r_node)

        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, r_node)
    }

    return err;

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
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
void ConservativeElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rElementalDofList.size() != mLocalSize)
        rElementalDofList.resize(mLocalSize);

    const GeometryType& r_geom = this->GetGeometry();
    int counter=0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = r_geom[i].pGetDof(MOMENTUM_X);
        rElementalDofList[counter++] = r_geom[i].pGetDof(MOMENTUM_Y);
        rElementalDofList[counter++] = r_geom[i].pGetDof(HEIGHT);
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
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
void ConservativeElement<TNumNodes>::AddWaveTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    const bool integrate_by_parts = true;
    const auto z = rData.topography;
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];
    const double c2 = rData.gravity * rData.height;
    const double l = this->StabilizationParameter(rData);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const IndexType i_block = 3 * i;
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const IndexType j_block = 3 * j;

            double g1_ij;
            double g2_ij;
            double d_ij;
            if (integrate_by_parts) {
                g1_ij = -rDN_DX(i,0) * rN[j];
                g2_ij = -rDN_DX(i,1) * rN[j];
            } else {
                g1_ij = rN[i] * rDN_DX(j,0);
                g2_ij = rN[i] * rDN_DX(j,1);
            }

            /* First component
             * A_1 = {{2 * u_1   0   -u_1^2 + gh},
             *        {  u_2    u_1   -u_1 * u_2},
             *        {   1      0        0     }}
             * b_1 = {gh, 0, 0}^T
             */
            rMatrix(i_block,     j_block)     += Weight * g1_ij * 2*u_1;
            rMatrix(i_block,     j_block + 2) += Weight * g1_ij * (-u_1*u_1 + c2);
            rMatrix(i_block + 1, j_block)     += Weight * g1_ij * u_2;
            rMatrix(i_block + 1, j_block + 1) += Weight * g1_ij * u_1;
            rMatrix(i_block + 1, j_block + 2) -= Weight * g1_ij * u_1*u_2;
            rMatrix(i_block + 2, j_block)     += Weight * g1_ij;
            rVector[i_block]                  -= Weight * g1_ij * c2 * z[j];

            /* Second component
             * A_2 = {{u_2    u_1      -u_1 * u_2},
             *        { 0   2 * u_2   -u_2^2 + gh},
             *        { 0      1            0    }}
             * b_2 = {0, gh, 0}^T
             */
            rMatrix(i_block,     j_block)     += Weight * g2_ij * u_2;
            rMatrix(i_block,     j_block + 1) += Weight * g2_ij * u_1;
            rMatrix(i_block,     j_block + 2) -= Weight * g2_ij * u_1*u_2;
            rMatrix(i_block + 1, j_block + 1) += Weight * g2_ij * 2*u_2;
            rMatrix(i_block + 1, j_block + 2) += Weight * g2_ij * (-u_2*u_2 + c2);
            rMatrix(i_block + 2, j_block + 1) += Weight * g2_ij;
            rVector[i_block + 1]              -= Weight * g2_ij * c2 * z[j];

            /* Stabilization x-x
             * A1*A1
             * A1*b1
             */
            d_ij = rDN_DX(i,0) * rDN_DX(j,0);
            rMatrix(i_block,     j_block)     += Weight * l * d_ij * (3*pow(u_1,2) + c2);
            rMatrix(i_block,     j_block + 2) += Weight * l * d_ij * (-2*pow(u_1,3) + 2*u_1*c2);
            rMatrix(i_block + 1, j_block)     += Weight * l * d_ij * 2*u_1*u_2;
            rMatrix(i_block + 1, j_block + 1) += Weight * l * d_ij * pow(u_1,2);
            rMatrix(i_block + 1, j_block + 2) += Weight * l * d_ij * (-2*pow(u_1,2)*u_2 + u_2*c2);
            rMatrix(i_block + 2, j_block)     += Weight * l * d_ij * 2*u_1;
            rMatrix(i_block + 2, j_block + 2) += Weight * l * d_ij * (-pow(u_1,2) + c2);
            rVector[i_block]                  -= Weight * l * d_ij * 2*u_1*c2 * z[j];
            rVector[i_block + 1]              -= Weight * l * d_ij * u_2*c2 * z[j];
            rVector[i_block + 2]              -= Weight * l * d_ij * c2 * z[j];

            /* Stabilization y-y
             * A2*A2
             * A2*b2
             */
            d_ij = rDN_DX(i,1) * rDN_DX(j,1);
            rMatrix(i_block,     j_block)     += Weight * l * d_ij * pow(u_2,2);
            rMatrix(i_block,     j_block + 1) += Weight * l * d_ij * 2*u_1*u_2;
            rMatrix(i_block,     j_block + 2) += Weight * l * d_ij * (-2*u_1*pow(u_2,2) + u_1*c2);
            rMatrix(i_block + 1, j_block + 1) += Weight * l * d_ij * (3*(pow(u_2,2) + c2));
            rMatrix(i_block + 1, j_block + 2) += Weight * l * d_ij * (-2*pow(u_2,3) + 2*u_2*c2);
            rMatrix(i_block + 2, j_block + 1) += Weight * l * d_ij * 2*u_2;
            rMatrix(i_block + 2, j_block + 2) += Weight * l * d_ij * (-pow(u_2,2) + c2);
            rVector[i_block]                  -= Weight * l * d_ij * u_1*c2 * z[j];
            rVector[i_block + 1]              -= Weight * l * d_ij * 2*u_2*c2 * z[j];
            rVector[i_block + 2]              -= Weight * l * d_ij * c2 * z[j];

            /* Stabilization x-y
             * A1*A2
             * A1*b1
             */
            d_ij = rDN_DX(i,0) * rDN_DX(j,1);
            rMatrix(i_block,     j_block)     += Weight * l * d_ij * 2*u_1*u_2;
            rMatrix(i_block,     j_block + 1) += Weight * l * d_ij * (pow(u_1,2)+c2);
            rMatrix(i_block,     j_block + 2) -= Weight * l * d_ij * 2*pow(u_1,2)*u_2;
            rMatrix(i_block + 1, j_block)     += Weight * l * d_ij * pow(u_2,2);
            rMatrix(i_block + 1, j_block + 1) += Weight * l * d_ij * 2*u_1*u_2;
            rMatrix(i_block + 1, j_block + 2) += Weight * l * d_ij * (-2*u_1*pow(u_2,2) + u_1*c2);
            rMatrix(i_block + 2, j_block)     += Weight * l * d_ij * u_2;
            rMatrix(i_block + 2, j_block + 1) += Weight * l * d_ij * u_1;
            rMatrix(i_block + 2, j_block + 2) -= Weight * l * d_ij * u_1*u_2;
            rVector[i_block + 1]              -= Weight * l * d_ij * u_1*c2 * z[j];

            /* Stabilization y-x
             * A2*A1
             * A2*b1
             */
            d_ij = rDN_DX(j,1) * rDN_DX(i,0);
            rMatrix(i_block,     j_block)     += Weight * l * d_ij * 2*u_1*u_2;
            rMatrix(i_block,     j_block + 1) += Weight * l * d_ij * pow(u_1,2);
            rMatrix(i_block,     j_block + 2) += Weight * l * d_ij * (-2*pow(u_1,2)*u_2 + u_2*c2);
            rMatrix(i_block + 1, j_block)     += Weight * l * d_ij * (pow(u_2,2) + c2);
            rMatrix(i_block + 1, j_block + 1) += Weight * l * d_ij * 2*u_1*u_2;
            rMatrix(i_block + 1, j_block + 2) -= Weight * l * d_ij * 2*u_1*pow(u_2,2);
            rMatrix(i_block + 2, j_block)     += Weight * l * d_ij * u_2;
            rMatrix(i_block + 2, j_block + 1) += Weight * l * d_ij * u_1;
            rMatrix(i_block + 2, j_block + 2) -= Weight * l * d_ij * u_1*u_2;
            rVector[i_block]                  -= Weight * l * d_ij * u_2*c2 * z[j];
        }
    }
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::AddFrictionTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    const auto u = rData.velocity;
    const double u_1 = u[0];
    const double u_2 = u[1];
    const double g = rData.gravity;
    const double s = rData.p_bottom_friction->CalculateLHS(rData.height, u);
    const double lumping_factor = 1.0 / TNumNodes;
    const double l = StabilizationParameter(rData);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const IndexType i_block = 3 * i;

        rMatrix(i_block,     i_block)     += Weight * lumping_factor * g*s;
        rMatrix(i_block + 1, i_block + 1) += Weight * lumping_factor * g*s;

        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const IndexType j_block = 3 * j;

            /* Stabilization x
             * l*A1*Sf
             */
            const double g1_ij = -rDN_DX(i,0) * rN[j];
            rMatrix(i_block,     j_block)     += Weight * l * g1_ij * 2*g*s*u_1;
            rMatrix(i_block + 1, j_block)     += Weight * l * g1_ij * g*s*u_2;
            rMatrix(i_block + 1, j_block + 1) += Weight * l * g1_ij * g*s*u_2;
            rMatrix(i_block + 2, j_block)     += Weight * l * g1_ij * g*s;

            /* Stabilization y
             * l*A2*Sf
             */
            const double g2_ij = -rDN_DX(i,1) * rN[j];
            rMatrix(i_block,     j_block)     += Weight * l * g2_ij * g*s*u_2;
            rMatrix(i_block,     j_block + 1) += Weight * l * g2_ij * g*s*u_1;
            rMatrix(i_block + 1, j_block + 1) += Weight * l * g2_ij * 2*g*s*u_2;
            rMatrix(i_block + 2, j_block + 1) += Weight * l * g2_ij * g*s;
        }
    }
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::AddArtificialViscosityTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];
    const double c2 = rData.gravity * rData.height;
    const double l = this->StabilizationParameter(rData);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const IndexType i_block = 3 * i;
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const IndexType j_block = 3 * j;

            /* Inertia terms
             */
            const double n_ij = rN[i] * rN[j];
            rMatrix(i_block,     j_block)     += Weight * n_ij;
            rMatrix(i_block + 1, j_block + 1) += Weight * n_ij;
            rMatrix(i_block + 2, j_block + 2) += Weight * n_ij;

            /* Stabilization x
             * l*A1*N
             */
            const double g1_ij = -rDN_DX(i,0) * rN[j];
            rMatrix(i_block,     j_block)     += Weight * l * g1_ij * 2*u_1;
            rMatrix(i_block,     j_block + 2) += Weight * l * g1_ij * (-u_1*u_1 + c2);
            rMatrix(i_block + 1, j_block)     += Weight * l * g1_ij * u_2;
            rMatrix(i_block + 1, j_block + 1) += Weight * l * g1_ij * u_1;
            rMatrix(i_block + 1, j_block + 2) -= Weight * l * g1_ij * u_1*u_2;
            rMatrix(i_block + 2, j_block)     += Weight * l * g1_ij;

            /* Stabilization y
             * l*A2*N
             */
            const double g2_ij = -rDN_DX(i,1) * rN[j];
            rMatrix(i_block,     j_block)     += Weight * l * g2_ij * u_2;
            rMatrix(i_block,     j_block + 1) += Weight * l * g2_ij * u_1;
            rMatrix(i_block,     j_block + 2) -= Weight * l * g2_ij * u_1*u_2;
            rMatrix(i_block + 1, j_block + 1) += Weight * l * g2_ij * 2*u_2;
            rMatrix(i_block + 1, j_block + 2) += Weight * l * g2_ij * (-u_2*u_2 + c2);
            rMatrix(i_block + 2, j_block + 1) += Weight * l * g2_ij;
        }
    }
}

template<std::size_t TNumNodes>
double ConservativeElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    const double lambda = std::sqrt(rData.gravity * rData.height) + norm_2(rData.velocity);
    const double epsilon = 1e-6;
    const double threshold = rData.relative_dry_height * rData.length;
    const double w = ShallowWaterUtilities().WetFraction(rData.height, threshold);
    return w * rData.length * rData.stab_factor / (lambda + epsilon);
}


template class ConservativeElement<3>;

} // namespace Kratos
