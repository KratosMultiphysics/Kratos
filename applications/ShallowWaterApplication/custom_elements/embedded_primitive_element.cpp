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
//                   Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "embedded_primitive_element.h"
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

// Application includes
#include "custom_utilities/phase_function.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
void EmbeddedPrimitiveElement<TNumNodes>::UpdateGaussPointData(ElementData& rData, const array_1d<double,TNumNodes>& rN)
{
    const double h = inner_prod(rData.nodal_h, rN);
    const double g = rData.gravity;
    const array_1d<double,3> v = BaseType::VectorProduct(rData.nodal_v, rN);

    rData.height = h;
    rData.velocity = v;

    /**
     * A_0 = {{ u_0  0     g},
     *        { 0    u_0   0},
     *        { h    0     u_0}}
     */
    rData.A1(0,0) = v[0];
    rData.A1(0,1) = 0;
    rData.A1(0,2) = g;
    rData.A1(1,0) = 0;
    rData.A1(1,1) = v[0];
    rData.A1(1,2) = 0;
    rData.A1(2,0) = h;
    rData.A1(2,1) = 0;
    rData.A1(2,2) = v[0];

    /*
     * A_1 = {{u_1    0      0},
     *        { 0    u_1     g},
     *        { 0     h      u_1}}
     */
    rData.A2(0,0) = v[1];
    rData.A2(0,1) = 0;
    rData.A2(0,2) = 0;
    rData.A2(1,0) = 0;
    rData.A2(1,1) = v[1];
    rData.A2(1,2) = g;
    rData.A2(2,0) = 0;
    rData.A2(2,1) = h;
    rData.A2(2,2) = v[1];

    /// b_1
    rData.b1[0] = g;
    rData.b1[1] = 0;
    rData.b1[2] = 0;

    /// b_2
    rData.b2[0] = 0;
    rData.b2[1] = g;
    rData.b2[2] = 0;
}

template<std::size_t TNumNodes>
double EmbeddedPrimitiveElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    // const double threshold = rData.relative_dry_height * rData.length;
    const double threshold = 1e-6 * rData.length;
    const double height = std::max(rData.height, threshold);
    const double lambda = std::sqrt(rData.gravity * height) + norm_2(rData.velocity);
    return rData.length * rData.stab_factor / lambda;
}

template<std::size_t TNumNodes>
int EmbeddedPrimitiveElement<TNumNodes>::CalculateDistances(Vector& rDistances) const
{
    if (rDistances.size() != TNumNodes) {
        rDistances.resize(TNumNodes, false);
    }
    std::size_t n_pos = 0;
    std::size_t n_neg = 0;
    auto& r_geom = this->GetGeometry();
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        const double d = r_geom[i].FastGetSolutionStepValue(DISTANCE);
        rDistances[i] = d;
        if (d > 0.0) {
            n_pos++;
        } else {
            n_neg++;
        }
    }
    if (n_pos == TNumNodes) {
        return 0;
    } else if (n_pos > 0 && n_neg > 0) {
        return 1;
    } else {
        return 2;
    }
}

//TODO: Remove the rGeometry argument (before it was static but now it's no longer required)
template<std::size_t TNumNodes>
void EmbeddedPrimitiveElement<TNumNodes>::CalculateGeometryData(
    const GeometryType& rGeometry,
    Vector &rGaussWeights,
    Matrix &rNContainer,
    ShapeFunctionsGradientsType &rDN_DX)
{
    Vector distances(TNumNodes);
    int is_cut = CalculateDistances(distances);

    if (is_cut == 0) {
        BaseType::CalculateGeometryData(rGeometry, rGaussWeights, rNContainer, rDN_DX);
    }
    else if (is_cut == 1) {
        auto p_mod_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), distances);

        p_mod_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            rNContainer,
            rDN_DX,
            rGaussWeights,
            GeometryData::IntegrationMethod::GI_GAUSS_2);
    }
}

template<std::size_t TNumNodes>
void EmbeddedPrimitiveElement<TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Add volumetric contribution from the base implementation
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Add the weak BC imposition over the interface
    //TODO: We're temporarily doing the splitting twice
    // Set the modified shape functions pointer
    const auto& r_geom = this->GetGeometry();
    Vector distances(TNumNodes);
    int is_cut = CalculateDistances(distances);

    if (is_cut == 1) {
        // Set element data container
        ElementData data;
        BaseType::InitializeData(data, rCurrentProcessInfo);
        BaseType::GetNodalData(data, r_geom);

        auto p_mod_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), distances);

        // Fluid side interface
        Vector interface_w;
        Matrix interface_N;
        ShapeFunctionsGradientsType interface_DN_DX;
        p_mod_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
            interface_N,
            interface_DN_DX,
            interface_w,
            GeometryData::IntegrationMethod::GI_GAUSS_2);

        // Fluid side interface normals
        ModifiedShapeFunctions::AreaNormalsContainerType interface_unit_n;
        p_mod_sh_func->ComputePositiveSideInterfaceAreaNormals(
            interface_unit_n,
            GeometryData::IntegrationMethod::GI_GAUSS_2);

        // Normalize the normals
        double h = ElementSizeCalculator<2,TNumNodes>::MinimumElementSize(this->GetGeometry());
        const double tolerance = 1.0e-3 * h;
        for (unsigned int i = 0; i < interface_unit_n.size(); ++i) {
            double norm = norm_2(interface_unit_n[i]);
            interface_unit_n[i] /= std::max(norm,tolerance);
        }

        // Nitsche penalty term imposition
        const double kappa = 1.0e-4;
        const double aux_pen = kappa / data.length;
        const double h_imposed = 1.0e-12;
        array_1d<double,3> aux_v;
        array_1d<double,3> aux_normal;
        array_1d<double,TNumNodes> aux_N;
        BoundedMatrix<double, TNumNodes, 2> aux_DN_DX;
        std::size_t n_gauss_int = interface_w.size();
        for (std::size_t g = 0; g < n_gauss_int; ++g) {
            // Get Gauss pt. data
            const double w_g = interface_w[g];
            noalias(aux_N) = row(interface_N, g);
            noalias(aux_DN_DX) = interface_DN_DX[g];
            noalias(aux_normal) = interface_unit_n[g];

            // Compute velocity norm
            noalias(aux_v) = ZeroVector(3);
            for (std::size_t j = 0; j < TNumNodes; ++j) {
                aux_v += aux_N(j) * r_geom[j].FastGetSolutionStepValue(VELOCITY);
            }

            // Assemble Nitsche terms
            const double penalty = aux_pen * norm_2(aux_v);
            for (std::size_t i = 0; i < TNumNodes; ++i) {
                std::size_t i_block = 3*i;
                for (std::size_t j = 0; j < TNumNodes; ++j) {
                    std::size_t j_block = 3*j;
                    const double h_j = r_geom[j].FastGetSolutionStepValue(HEIGHT);
                    // Height equation penalty term
                    rRightHandSideVector(i_block + 2) += w_g * penalty * aux_N(i) * (aux_N(j) * h_j - h_imposed);
                    rLeftHandSideMatrix(i_block + 2, j_block + 2) -= w_g * penalty * aux_N(i) * aux_N(j);
                    // Velocity equation stabilization term
                    rRightHandSideVector(i_block) += w_g * aux_N(i) * aux_normal[0] * (aux_N(j) * h_j - h_imposed);
                    rRightHandSideVector(i_block + 1) += w_g * aux_N(i) * aux_normal[1] * (aux_N(j) * h_j - h_imposed);
                    rLeftHandSideMatrix(i_block, j_block + 2) -= w_g * aux_N(i) * aux_normal[0] * aux_N(j);
                    rLeftHandSideMatrix(i_block + 1, j_block + 2) -= w_g * aux_N(i) * aux_normal[1] * aux_N(j);
                }
            }
        }
    }
}


template class EmbeddedPrimitiveElement<3>;

} // namespace Kratos
