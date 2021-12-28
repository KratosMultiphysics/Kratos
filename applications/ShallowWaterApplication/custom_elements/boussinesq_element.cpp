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
#include "boussinesq_element.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/phase_function.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
const Variable<double>& BoussinesqElement<TNumNodes>::GetUnknownComponent(int Index) const
{
    switch (Index) {
        case 0: return VELOCITY_X;
        case 1: return VELOCITY_Y;
        case 2: return FREE_SURFACE_ELEVATION;
        default: KRATOS_ERROR << "BoussinesqElement::GetUnknownComponent index out of bounds." << std::endl;
    }
}

template<std::size_t TNumNodes>
typename BoussinesqElement<TNumNodes>::LocalVectorType BoussinesqElement<TNumNodes>::GetUnknownVector(const ElementData& rData) const
{
    std::size_t index = 0;
    array_1d<double,mLocalSize> unknown;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        unknown[index++] = rData.nodal_v[i][0];
        unknown[index++] = rData.nodal_v[i][1];
        unknown[index++] = rData.nodal_f[i];
    }
    return unknown;
}

template<std::size_t TNumNodes>
typename BoussinesqElement<TNumNodes>::LocalVectorType BoussinesqElement<TNumNodes>::ConservativeVector(
    const LocalVectorType& rVector,
    const ElementData& rData) const
{
    array_1d<double,mLocalSize> conservative_vector;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        conservative_vector[3*i    ] = -rVector[3*i    ] * rData.nodal_z[i];
        conservative_vector[3*i + 1] = -rVector[3*i + 1] * rData.nodal_z[i];
        conservative_vector[3*i + 2] =  rVector[3*i + 2];
    }
    return conservative_vector;
}

template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::GetNodalData(ElementData& rData, const GeometryType& rGeometry, int Step)
{
    rData.length = rGeometry.Length();

    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rData.nodal_f[i] = rGeometry[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION, Step);
        rData.nodal_h[i] = rGeometry[i].FastGetSolutionStepValue(HEIGHT, Step);
        rData.nodal_z[i] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY, Step);
        rData.nodal_v[i] = rGeometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        rData.nodal_a[i] = rGeometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
        rData.nodal_v_lap[i] = rGeometry[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN, Step);
    }
}

template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::CalculateGaussPointData(ElementData& rData, const array_1d<double,TNumNodes>& rN)
{
    const double eta = inner_prod(rData.nodal_f, rN);
    const double H = -inner_prod(rData.nodal_z, rN);
    const double g = rData.gravity;
    const array_1d<double,3> v = WaveElementType::VectorProduct(rData.nodal_v, rN);

    double e; // the non linearity ratio
    if (rData.amplitude) {
        e = rData.amplitude / H;
    } else {
        e = std::abs(eta) / H;
    }

    rData.depth = H;
    rData.height = H + e * eta;
    rData.velocity = v;

    /**
     * A_1 = {{ e*u_1      0     g  },
     *        {   0      e*u_1   0  },
     *        {H + e*eta   0   e*u_1}}
     */
    rData.A1(0,0) = e*v[0];
    rData.A1(0,1) = 0;
    rData.A1(0,2) = g;
    rData.A1(1,0) = 0;
    rData.A1(1,1) = e*v[0];
    rData.A1(1,2) = 0;
    rData.A1(2,0) = H + e*eta;
    rData.A1(2,1) = 0;
    rData.A1(2,2) = e*v[0];

    /*
     * A_2 = {{ e*u_2      0      0  },
     *        {   0      e*u_2    g  },
     *        {   0   H + e*eta e*u_2}}
     */
    rData.A2(0,0) = e*v[1];
    rData.A2(0,1) = 0;
    rData.A2(0,2) = 0;
    rData.A2(1,0) = 0;
    rData.A2(1,1) = e*v[1];
    rData.A2(1,2) = g;
    rData.A2(2,0) = 0;
    rData.A2(2,1) = H + e*eta;
    rData.A2(2,2) = e*v[1];

    /// b_1
    rData.b1[0] = 0;
    rData.b1[1] = 0;
    rData.b1[2] = -v[0];

    /// b_2
    rData.b2[0] = 0;
    rData.b2[1] = 0;
    rData.b2[2] = -v[1];
}

template<std::size_t TNumNodes>
double BoussinesqElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    const double lambda = std::sqrt(rData.gravity * std::abs(rData.height)) + norm_2(rData.velocity);
    const double epsilon = 1e-6;
    const double threshold = rData.relative_dry_height * rData.length;
    const double w = PhaseFunction::WetFraction(rData.height, threshold);
    return w * rData.length * rData.stab_factor / (lambda + epsilon);
}

template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddDispersiveTerms(
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    // Constants
    const double dispersive_ratio = rData.depth / rData.wavelength;
    const double m2 = std::pow(dispersive_ratio, 2);
    const double beta = -0.531;
    const double C1 = 0.5 * std::pow(beta, 2) - 0.166666666666;
    const double C3 = beta + 0.5;
    const double H = rData.depth;

    // Projected auxiliary field
    array_1d<array_1d<double,3>,TNumNodes> nodal_Jf;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        nodal_Jf[i] = (C1 + C3) * std::pow(H, 3) * rData.nodal_v_lap[i];
    }

    // Stabilization constants
    const double l = StabilizationParameter(rData);
    const array_1d<double,3> A1i3 = column(rData.A1, 2);
    const array_1d<double,3> A2i3 = column(rData.A2, 2);

    // Adding the contribution of the dispersive fields
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            double g1_ij;
            double g2_ij;
            double d_ij;
            if (rData.integrate_by_parts) {
                g1_ij = -rDN_DX(i,0) * rN[j];
                g2_ij = -rDN_DX(i,1) * rN[j];
            } else {
                g1_ij = rN[i] * rDN_DX(j,0);
                g2_ij = rN[i] * rDN_DX(j,1);
            }

            /// Gradient contribution
            rVector[3*i + 2] -= Weight*m2*g1_ij*nodal_Jf[j][0];
            rVector[3*i + 2] -= Weight*m2*g2_ij*nodal_Jf[j][1];

            /// Stabilization x-x
            d_ij = rDN_DX(i,0) * rDN_DX(j,0);
            MathUtils<double>::AddVector(rVector, -Weight*m2*l*d_ij*A1i3*nodal_Jf[j][0], 3*i);

            /// Stabilization y-y
            d_ij = rDN_DX(i,1) * rDN_DX(j,1);
            MathUtils<double>::AddVector(rVector, -Weight*m2*l*d_ij*A2i3*nodal_Jf[j][1], 3*i);

            /// Stabilization x-y
            d_ij = rDN_DX(i,0) * rDN_DX(j,1);
            MathUtils<double>::AddVector(rVector, -Weight*m2*l*d_ij*A1i3*nodal_Jf[j][1], 3*i);

            /// Stabilization y-x
            d_ij = rDN_DX(i,1) * rDN_DX(j,0);
            MathUtils<double>::AddVector(rVector, -Weight*m2*l*d_ij*A2i3*nodal_Jf[j][0], 3*i);
        }
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    // Constants
    const double dispersive_ratio = rData.depth / rData.wavelength;
    const double m2 = std::pow(dispersive_ratio, 2);
    const double beta = -0.531;
    const double C2 = 0.5 * std::pow(beta, 2);
    const double C4 = beta;
    const double H = rData.depth;
    const double l = StabilizationParameter(rData);
    BoundedMatrix<double,3,3> M = IdentityMatrix(3);
    BoundedMatrix<double,2,2> K;
    array_1d<double,2> derivatives_i;
    array_1d<double,2> derivatives_j;

    // Adding the contribution of the dispersive fields
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        derivatives_i[0] = rDN_DX(i,0);
        derivatives_i[1] = rDN_DX(i,1);
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            /// Consistent mass matrix
            const double n_ij = WaveElementType::ShapeFunctionProduct(rN, i, j);//rN[i] * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*n_ij*M, 3*i, 3*j);

            /// Stabilization x
            const double g1_ij = rDN_DX(i,0) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g1_ij*rData.A1, 3*i, 3*j);

            /// Stabilization y
            const double g2_ij = rDN_DX(i,1) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g2_ij*rData.A2, 3*i, 3*j);

            /// Disperive term
            derivatives_j[0] = rDN_DX(j,0);
            derivatives_j[1] = rDN_DX(j,1);
            noalias(K) = -outer_prod(derivatives_i, derivatives_j);
            MathUtils<double>::AddMatrix(rMatrix, Weight*m2*(C2+C4)*std::pow(H,2)*K, 3*i, 3*j);
        }
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddAuxiliaryLaplacian(
    LocalMatrixType& rLaplacian,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    array_1d<double,3> gradients_vector_i;
    array_1d<double,3> gradients_vector_j;

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        gradients_vector_i[0] = rDN_DX(i,0);
        gradients_vector_i[1] = rDN_DX(i,1);
        gradients_vector_i[2] = 0.0;

        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            gradients_vector_j[0] = rDN_DX(j,0);
            gradients_vector_j[1] = rDN_DX(j,1);
            gradients_vector_j[2] = 0.0;

            // Projector for the auxiliary field
            MathUtils<double>::AddMatrix(rLaplacian, -Weight*outer_prod(gradients_vector_i, gradients_vector_j), 3*i, 3*j);
        }
    }
}


template<>
void BoussinesqElement<3>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    auto& r_geom = this->GetGeometry();

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, r_geom);

    // Geometrical data
    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, area);

    // Auxiliary fields
    LocalMatrixType laplacian = ZeroMatrix(mLocalSize, mLocalSize);

    // Gauss point contribution
    CalculateGaussPointData(data, N);
    AddAuxiliaryLaplacian(laplacian, data, N, DN_DX, area);

    // Calculate the lapacian vector
    const LocalVectorType& vel_vector = this->GetUnknownVector(data);
    LocalVectorType vel_laplacian_vector = prod(laplacian, vel_vector);

    // Add the elemental contribution to the nodes
    array_1d<double,3> vel_laplacian = ZeroVector(3);
    for (std::size_t i = 0; i < 3; ++i)
    {
        std::size_t block = 3 * i;
        vel_laplacian[0] = vel_laplacian_vector[block];
        vel_laplacian[1] = vel_laplacian_vector[block + 1];
        r_geom[i].SetLock();
        r_geom[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN) += vel_laplacian;
        r_geom[i].UnSetLock();
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddRightHandSide(
    LocalVectorType& rRHS,
    ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    LocalMatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);

    CalculateGaussPointData(rData, rN);

    this->AddWaveTerms(lhs, rRHS, rData, rN, rDN_DX, Weight);
    this->AddFrictionTerms(lhs, rRHS, rData, rN, rDN_DX, Weight);
    this->AddDispersiveTerms(rRHS, rData, rN, rDN_DX, Weight);
    this->AddArtificialViscosityTerms(lhs, rData, rDN_DX, Weight);

    noalias(rRHS) -= prod(lhs, this->GetUnknownVector(rData));
}


template<>
void BoussinesqElement<3>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    auto& r_geom = this->GetGeometry();

    LocalVectorType f0 = ZeroVector(mLocalSize);
    LocalVectorType f1 = ZeroVector(mLocalSize);
    LocalVectorType f2 = ZeroVector(mLocalSize);
    LocalVectorType f3 = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, area);

    GetNodalData(data, r_geom, 0);
    AddRightHandSide(f0, data, N, DN_DX, area);

    GetNodalData(data, r_geom, 1);
    AddRightHandSide(f1, data, N, DN_DX, area);

    GetNodalData(data, r_geom, 2);
    AddRightHandSide(f2, data, N, DN_DX, area);

    GetNodalData(data, r_geom, 3);
    AddRightHandSide(f3, data, N, DN_DX, area);

    noalias(rRightHandSideVector) = (9*f0 + 19*f1 - 5*f2 + f3) / 24.0;
}


template<>
void BoussinesqElement<3>::AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo)
{
    auto& r_geom = this->GetGeometry();

    LocalVectorType f1 = ZeroVector(mLocalSize);
    LocalVectorType f2 = ZeroVector(mLocalSize);
    LocalVectorType f3 = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, area);

    GetNodalData(data, r_geom, 1);
    AddRightHandSide(f1, data, N, DN_DX, area);

    GetNodalData(data, r_geom, 2);
    AddRightHandSide(f2, data, N, DN_DX, area);

    GetNodalData(data, r_geom, 3);
    AddRightHandSide(f3, data, N, DN_DX, area);

    LocalVectorType increment = (23*f1 - 16*f2 + 5*f3) / 12.0;
    array_1d<double,3> nodal_increment;
    for (std::size_t i = 0; i < r_geom.size(); ++i)
    {
        std::size_t block = 3*i;
        nodal_increment[0] = increment[block];
        nodal_increment[1] = increment[block + 1];
        nodal_increment[2] = increment[block + 2];

        r_geom[i].SetLock();
        r_geom[i].FastGetSolutionStepValue(RESIDUAL_VECTOR) += nodal_increment;
        r_geom[i].UnSetLock();
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if(rLumpedMassVector.size() != mLocalSize)
        rLumpedMassVector.resize(mLocalSize, false);

    const double area = this->GetGeometry().Area();
    const double lump_factor = 1.0 / static_cast<double>(TNumNodes);
    rLumpedMassVector = Vector(mLocalSize, area * lump_factor);
}


template class BoussinesqElement<3>;

} // namespace Kratos
