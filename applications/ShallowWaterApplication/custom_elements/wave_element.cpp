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
#include "utilities/geometry_utilities.h"
#include "custom_friction_laws/manning_law.h"
#include "custom_utilities/shallow_water_utilities.h"
#include "shallow_water_application_variables.h"
#include "wave_element.h"

namespace Kratos
{

template<std::size_t TNumNodes>
int WaveElement<TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int err = Element::Check(rCurrentProcessInfo);
    if (err != 0) return err;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry())
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MANNING, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VERTICAL_VELOCITY, r_node)

        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, r_node)
    }

    return err;

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rResult.size() != mLocalSize)
        rResult.resize(mLocalSize, false); // False says not to preserve existing storage!!

    const GeometryType& r_geom = this->GetGeometry();
    int counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(1)).EquationId();
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(2)).EquationId();
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(3)).EquationId();
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rElementalDofList.size() != mLocalSize)
        rElementalDofList.resize(mLocalSize);

    const GeometryType& r_geom = this->GetGeometry();
    int counter=0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(1));
        rElementalDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(2));
        rElementalDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(3));
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(1), Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(2), Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(3), Step);
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(VERTICAL_VELOCITY, Step);
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_ERROR << "WaveElement::GetSecondDerivativesVector This method is not supported by the formulation" << std::endl;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType PointNumber = 0; PointNumber < 1; PointNumber++)
        rValues[PointNumber] = double(this->GetValue(rVariable));
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
}

template<std::size_t TNumNodes>
const Variable<double>& WaveElement<TNumNodes>::GetUnknownComponent(int Index) const
{
    switch (Index) {
        case 1: return VELOCITY_X;
        case 2: return VELOCITY_Y;
        case 3: return HEIGHT;
        default: KRATOS_ERROR << "WaveElement::GetUnknownComponent index out of bounds." << std::endl;
    }
}

template<std::size_t TNumNodes>
typename WaveElement<TNumNodes>::LocalVectorType WaveElement<TNumNodes>::GetUnknownVector(ElementData& rData)
{
    std::size_t index = 0;
    array_1d<double,mLocalSize> unknown;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        unknown[index++] = rData.nodal_v[i][0];
        unknown[index++] = rData.nodal_v[i][1];
        unknown[index++] = rData.nodal_h[i];
    }
    return unknown;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::InitializeData(ElementData& rData, const ProcessInfo& rCurrentProcessInfo)
{
    rData.stab_factor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    rData.shock_stab_factor = rCurrentProcessInfo[SHOCK_STABILIZATION_FACTOR];
    rData.relative_dry_height = rCurrentProcessInfo[RELATIVE_DRY_HEIGHT];
    rData.gravity = rCurrentProcessInfo[GRAVITY_Z];
    rData.p_bottom_friction = Kratos::make_shared<ManningLaw>();
    rData.p_bottom_friction->Initialize(this->GetGeometry(), rCurrentProcessInfo);
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetNodalData(ElementData& rData, const GeometryType& rGeometry, int Step)
{
    rData.length = rGeometry.Length();

    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rData.nodal_h[i] = rGeometry[i].FastGetSolutionStepValue(HEIGHT, Step);
        rData.nodal_z[i] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY, Step);
        rData.nodal_v[i] = rGeometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        rData.nodal_q[i] = rGeometry[i].FastGetSolutionStepValue(MOMENTUM, Step);

        // const IndexType block = 3 * i;

        // const auto h = rGeometry[i].FastGetSolutionStepValue(HEIGHT, Step);
        // const array_1d<double,3> v = rGeometry[i].FastGetSolutionStepValue(VELOCITY, Step);

        // rData.topography[i] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY, Step);

        // rData.unknown[block]     = v[0];
        // rData.unknown[block + 1] = v[1];
        // rData.unknown[block + 2] = h;
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateGaussPointData(
    ElementData& rData,
    const array_1d<double,TNumNodes>& rN)
{
    const double h = inner_prod(rData.nodal_h, rN);
    const array_1d<double,3> v = VectorProduct(rData.nodal_v, rN);

    rData.height = h;
    rData.velocity = v;

    /**
     * A_1 = {{ 0   0   g },
     *        { 0   0   0 },
     *        { h   0   0 }}
     */
    rData.A1 = ZeroMatrix(3, 3);
    rData.A1(0,2) = rData.gravity;
    rData.A1(2,0) = h;

    /*
     * A_2 = {{ 0   0   0 },
     *        { 0   0   g },
     *        { 0   h   0 }}
     */
    rData.A2 = ZeroMatrix(3, 3);
    rData.A2(1,2) = rData.gravity;
    rData.A2(2,1) = h;

    rData.b1 = ZeroVector(3);
    rData.b1[0] = rData.gravity;

    rData.b2 = ZeroVector(3);
    rData.b2[1] = rData.gravity;
    // rData.height = 0.0;
    // rData.velocity = ZeroVector(3);

    // for (IndexType i = 0; i < TNumNodes; i++)
    // {
    //     const IndexType block = 3 * i;

    //     const auto h = rData.unknown[block + 2];
    //     array_1d<double,3> v;
    //     v[0] = rData.unknown[block];
    //     v[1] = rData.unknown[block + 1];
    //     v[2] = 0.0;

    //     rData.height += rN[i] * h;
    //     rData.velocity += rN[i] * v;
    // }
}

template<std::size_t TNumNodes>
virtual void CalculateArtificialViscosityData(
    ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX)
{
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateGeometryData(
    Vector &rGaussWeights,
    Matrix &rNContainer,
    ShapeFunctionsGradientsType &rDN_DXContainer) const
{
    Vector det_j_vector;
    const auto& r_geom = this->GetGeometry();
    const auto integration_method = r_geom.GetDefaultIntegrationMethod();
    r_geom.ShapeFunctionsIntegrationPointsGradients(rDN_DXContainer, det_j_vector, integration_method, rNContainer);

    const unsigned int number_of_gauss_points = r_geom.IntegrationPointsNumber(integration_method);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);

    if (rGaussWeights.size() != number_of_gauss_points)
        rGaussWeights.resize(number_of_gauss_points, false);

    for (unsigned int g = 0; g < number_of_gauss_points; ++g)
        rGaussWeights[g] = det_j_vector[g] * integration_points[g].Weight();
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::AddWaveTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    const bool integrate_by_parts = false;
    // const auto h = rData.height;
    // const auto g = rData.gravity;
    const auto z = rData.nodal_z;
    // const double c = std::sqrt(g*h);
    const double l = StabilizationParameter(rData);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const range i_range (3 * i, 3 * i + 3);
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const range j_range (3 * j, 3 * j + 3);

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

            /// First component
            project(rMatrix, i_range, j_range) += Weight * g1_ij * rData.A1;
            project(rVector, i_range)          -= Weight * g1_ij * rData.b1 * z[j];

            /// Second component
            project(rMatrix, i_range, j_range) += Weight * g2_ij * rData.A2;
            project(rVector, i_range)          -= Weight * g2_ij * rData.b2 * z[j];

            /// Stabilization x-x
            d_ij = rDN_DX(i,0) * rDN_DX(j,0);
            project(rMatrix, i_range, j_range) += Weight * l * d_ij * prod(rData.A1, rData.A1);
            project(rVector, i_range)          -= Weight * l * d_ij * prod(rData.A1, rData.b1) * z[j];

            /// Stabilization y-y
            d_ij = rDN_DX(i,1) * rDN_DX(j,1);
            project(rMatrix, i_range, j_range) += Weight * l * d_ij * prod(rData.A2, rData.A2);
            project(rVector, i_range)          -= Weight * l * d_ij * prod(rData.A2, rData.b2) * z[j];

            /// Stabilization x-y
            d_ij = rDN_DX(i,0) * rDN_DX(j,1);
            project(rMatrix, i_range, j_range) += Weight * l * d_ij * prod(rData.A1, rData.A2);
            project(rVector, i_range)          -= Weight * l * d_ij * prod(rData.A1, rData.b2) * z[j];

            /// Stabilization y-x
            d_ij = rDN_DX(i,1) * rDN_DX(j,0);
            project(rMatrix, i_range, j_range) += Weight * l * d_ij * prod(rData.A2, rData.A1);
            project(rVector, i_range)          -= Weight * l * d_ij * prod(rData.A2, rData.b1) * z[j];

            // /* First component
            //  * A_1 = {{ 0   0   g },
            //  *        { 0   0   0 },
            //  *        { h   0   0 }}
            //  */
            // rMatrix(i_block,     j_block + 2) += Weight * g1_ij * g;
            // rMatrix(i_block + 2, j_block)     += Weight * g1_ij * h;
            // rVector[i_block]                  -= Weight * g1_ij * g * z[j];

            // /* Second component
            //  * A_2 = {{ 0   0   0 },
            //  *        { 0   0   g },
            //  *        { 0   h   0 }}
            //  */
            // rMatrix(i_block + 1, j_block + 2) += Weight * g2_ij * g;
            // rMatrix(i_block + 2, j_block + 1) += Weight * g2_ij * h;
            // rVector[i_block + 1]              -= Weight * g2_ij * g * z[j];

            // /* Stabilization x-x
            //  * l / sqrt(gh) * A1 * A1
            //  */
            // d_ij = -rDN_DX(i,0) * rDN_DX(j,0);
            // rMatrix(i_block,     j_block)     -= Weight * l * c * d_ij;
            // rMatrix(i_block + 2, j_block + 2) -= Weight * l * c * d_ij;
            // rVector[i_block + 2]              += Weight * l * c * d_ij * z[j];

            // /* Stabilization y-y
            //  * l / sqrt(gh) * A2 * A2
            //  */
            // d_ij = -rDN_DX(i,1) * rDN_DX(j,1);
            // rMatrix(i_block + 1, j_block + 1) -= Weight * l * c * d_ij;
            // rMatrix(i_block + 2, j_block + 2) -= Weight * l * c * d_ij;
            // rVector[i_block + 2]              += Weight * l * c * d_ij * z[j];

            // /* Stabilization x-y
            //  * l / sqrt(gh) * A1 * A2
            //  */
            // d_ij = -rDN_DX(i,1) * rDN_DX(j,0);
            // rMatrix(i_block,     j_block + 1) -= Weight * l * c * d_ij;

            // /* Stabilization y-x
            //  * l / sqrt(gh) * A1 * A2
            //  */
            // d_ij = -rDN_DX(i,0) * rDN_DX(j,1);
            // rMatrix(i_block + 1, j_block)     -= Weight * l * c * d_ij;
        }
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::AddFrictionTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    // const auto v = rData.velocity;
    // const double g = rData.gravity;
    // const double h = rData.height;
    const double s = rData.p_bottom_friction->CalculateLHS(rData.height, rData.velocity);
    const double lumping_factor = 1.0 / TNumNodes;
    // const double inv_c = std::sqrt(InverseHeight(rData) / g);
    const double l = StabilizationParameter(rData);
    BoundedMatrix<double,3,3> Sf = ZeroMatrix(3,3);
    Sf(0,0) = rData.gravity*s;
    Sf(1,1) = rData.gravity*s;

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        // const IndexType i_block = 3 * i;
        const range i_range (3*i, 3*i + 3);

        project(rMatrix, i_range, i_range) += Weight * lumping_factor * Sf;

        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const range j_range (3*j, 3*j + 3);

            /// Stabilization x
            const double g1_ij = rDN_DX(i,0) * rN[j];
            project(rMatrix, i_range, j_range) += Weight * l * g1_ij * prod(rData.A1, Sf);

            /// Stabilization y
            const double g2_ij = rDN_DX(i,1) * rN[j];
            project(rMatrix, i_range, j_range) += Weight * l * g2_ij * prod(rData.A2, Sf);
        }
    }

    //     rMatrix(i_block,     i_block)     += lumping_factor * g * s * Weight;
    //     rMatrix(i_block + 1, i_block + 1) += lumping_factor * g * s * Weight;

    //     for (IndexType j = 0; j < TNumNodes; ++j)
    //     {
    //         const IndexType j_block = 3 * j;

    //         /* Stabilization x
    //         * l / sqrt(gh) * A1 * Sf
    //         */
    //         const double g1_ij = -rN[j] * rDN_DX(i,0);
    //         rMatrix(i_block + 2, j_block)     -= l * g1_ij * h * inv_c * g * s * Weight;

    //         /* Stabilization y
    //         * l / sqrt(gh) * A2 * Sf
    //         */
    //         const double g2_ij = -rN[j] * rDN_DX(i,1);
    //         rMatrix(i_block + 2, j_block + 1) -= l * g2_ij * h * inv_c * g * s * Weight;
    //     }
    // }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::AddArtificialViscosityTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{}

template<>
void WaveElement<3>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX,
    const double Weight)
{
    // const double h = rData.height;
    // const double g = rData.gravity;
    // const double inv_c = std::sqrt(InverseHeight(rData) / g);
    const double l = StabilizationParameter(rData);
    BoundedMatrix<double, 3, 3> M = IdentityMatrix(3);

    // Algebraic factor
    const double one_twelfth = 1.0 / 12.0;
    const double one_sixth = 1.0 / 6.0;

    for (IndexType i = 0; i < 3; ++i)
    {
        const range i_range (3*i, 3*i + 3);
        for (IndexType j = 0; j < 3; ++j)
        {
            const range j_range (3*j, 3*j + 3);

            // Consistent mass matrix
            const double n = (i == j) ? one_sixth : one_twelfth;
            project(rMatrix, i_range, j_range) += Weight * n * M;

            /// Stabilization x
            const double g1_ij = rDN_DX(i,0) * rN[j];
            project(rMatrix, i_range, j_range) += Weight * l * g1_ij * rData.A1;

            /// Stabilization y
            const double g2_ij = rDN_DX(i,1) * rN[j];
            project(rMatrix, i_range, j_range) += Weight * l * g2_ij * rData.A2;
        }
    }

    // for (IndexType i = 0; i < 3; ++i)
    // {
    //     const IndexType i_block = 3 * i;
    //     for (IndexType j = 0; j < 3; ++j)
    //     {
    //         const IndexType j_block = 3 * j;

    //         // Algebraic mass matrix
    //         const double n = (i == j) ? 2*one_twelve : one_twelve;
    //         rMatrix(i_block,     j_block)     += n;
    //         rMatrix(i_block + 1, j_block + 1) += n;
    //         rMatrix(i_block + 2, j_block + 2) += n;

    //         /* Stabilization x
    //          * l / sqrt(gh) * A1 * N
    //          */
    //         const double g1_ij = -rN[j] * rDN_DX(i,0);
    //         rMatrix(i_block,     j_block + 2) -= l * g1_ij * g * inv_c;
    //         rMatrix(i_block + 2, j_block)     -= l * g1_ij * h * inv_c;

    //         /* Stabilization y
    //          * l / sqrt(gh) * A2 * N
    //          */
    //         const double g2_ij = -rN[j] * rDN_DX(i,1);
    //         rMatrix(i_block + 1, j_block + 2) -= l * g2_ij * g * inv_c;
    //         rMatrix(i_block + 2, j_block + 1) -= l * g2_ij * h * inv_c;
    //     }
    // }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    // const double h = rData.height;
    // const double g = rData.gravity;
    // const double inv_c = std::sqrt(InverseHeight(rData) / g);
    const double l = StabilizationParameter(rData);
    BoundedMatrix<double, 3, 3> M = IdentityMatrix(3);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const range i_range (3*i, 3*i + 3);
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            const range j_range (3*j, 3*j + 3);

            // Consistent mass matrix
            const double n_ij = rN[i] * rN[j];
            project(rMatrix, i_range, j_range) += Weight * n_ij * M;

            /// Stabilization x
            const double g1_ij = rDN_DX(i,0) * rN[j];
            project(rMatrix, i_range, j_range) += Weight * l * g1_ij * rData.A1;

            /// Stabilization y
            const double g2_ij = rDN_DX(i,1) * rN[j];
            project(rMatrix, i_range, j_range) += Weight * l * g2_ij * rData.A2;
        }
    }

    // for (IndexType i = 0; i < TNumNodes; ++i)
    // {
    //     const IndexType i_block = 3 * i;
    //     for (IndexType j = 0; j < TNumNodes; ++j)
    //     {
    //         const IndexType j_block = 3 * j;

    //         /* Inertia terms
    //          */
    //         const double n_ij = rN[i] * rN[j];
    //         rMatrix(i_block,     j_block)     += Weight * n_ij;
    //         rMatrix(i_block + 1, j_block + 1) += Weight * n_ij;
    //         rMatrix(i_block + 2, j_block + 2) += Weight * n_ij;

    //         /* Stabilization x
    //          * l / sqrt(gh) * A1 * N
    //          */
    //         const double g1_ij = -rN[j] * rDN_DX(i,0);
    //         rMatrix(i_block,     j_block + 2) -= Weight * l * g1_ij * g * inv_c;
    //         rMatrix(i_block + 2, j_block)     -= Weight * l * g1_ij * h * inv_c;

    //         /* Stabilization y
    //          * l / sqrt(gh) * A2 * N
    //          */
    //         const double g2_ij = -rN[j] * rDN_DX(i,1);
    //         rMatrix(i_block + 1, j_block + 2) -= Weight * l * g2_ij * g * inv_c;
    //         rMatrix(i_block + 2, j_block + 1) -= Weight * l * g2_ij * h * inv_c;
    //     }
    // }
}

template<std::size_t TNumNodes>
double WaveElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    const double inv_c = std::sqrt(InverseHeight(rData) / rData.gravity);
    return rData.length * rData.stab_factor * inv_c;
}

template<std::size_t TNumNodes>
double WaveElement<TNumNodes>::InverseHeight(const ElementData& rData) const
{
    const double height = rData.height;
    const double epsilon = rData.relative_dry_height * rData.length;
    return ShallowWaterUtilities().InverseHeight(height, epsilon);
}

template<std::size_t TNumNodes>
const array_1d<double,3> WaveElement<TNumNodes>::VectorProduct(
    const array_1d<array_1d<double,3>,TNumNodes>& rV,
    const array_1d<double,TNumNodes>& rN) const
{
    array_1d<double,3> result = ZeroVector(3);
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        result += rV[i] * rN[i];
    }
    return result;
}

template<>
void WaveElement<3>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != mLocalSize)
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);

    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    LocalMatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType rhs = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, this->GetGeometry());

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, area);

    CalculateGaussPointData(data, N);

    AddWaveTerms(lhs, rhs, data, N, DN_DX);
    AddFrictionTerms(lhs, rhs, data, N, DN_DX);
    AddArtificialViscosityTerms(lhs, data, DN_DX);

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rhs) -= prod(lhs, this->GetUnknownVector(data));

    noalias(rLeftHandSideMatrix) = area * lhs;
    noalias(rRightHandSideVector) = area * rhs;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != mLocalSize)
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);

    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    LocalMatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType rhs = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, this->GetGeometry());

    Vector weights;
    Matrix N_container;
    ShapeFunctionsGradientsType DN_DX_container;
    CalculateGeometryData(weights, N_container, DN_DX_container);
    const IndexType num_gauss_points = weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);
        const BoundedMatrix<double,TNumNodes,2> DN_DX = DN_DX_container[g];

        CalculateGaussPointData(data, N);

        AddWaveTerms(lhs, rhs, data, N, DN_DX, weight);
        AddFrictionTerms(lhs, rhs, data, N, DN_DX, weight);
        AddArtificialViscosityTerms(lhs, data, DN_DX, weight);
    }

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rhs) -= prod(lhs, this->GetUnknownVector(data));

    noalias(rLeftHandSideMatrix) = lhs;
    noalias(rRightHandSideVector) = rhs;
}

template<>
void WaveElement<3>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if(rMassMatrix.size1() != mLocalSize)
        rMassMatrix.resize(mLocalSize, mLocalSize, false);

    LocalMatrixType mass_matrix = ZeroMatrix(mLocalSize, mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, this->GetGeometry());

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, area);

    CalculateGaussPointData(data, N);

    AddMassTerms(mass_matrix, data, N, DN_DX);

    noalias(rMassMatrix) = area * mass_matrix;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if(rMassMatrix.size1() != mLocalSize)
        rMassMatrix.resize(mLocalSize, mLocalSize, false);

    LocalMatrixType mass_matrix = ZeroMatrix(mLocalSize, mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, this->GetGeometry());

    Vector weights;
    Matrix N_container;
    ShapeFunctionsGradientsType DN_DX_container;
    CalculateGeometryData(weights, N_container, DN_DX_container);
    const IndexType num_gauss_points = weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);
        const BoundedMatrix<double,TNumNodes,2> DN_DX = DN_DX_container[g];

        CalculateGaussPointData(data, N);

        AddMassTerms(mass_matrix, data, N, DN_DX, weight);
    }
    noalias(rMassMatrix) = mass_matrix;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

template class WaveElement<3>;
template class WaveElement<4>;
template class WaveElement<6>;
template class WaveElement<8>;
template class WaveElement<9>;

} // namespace Kratos
