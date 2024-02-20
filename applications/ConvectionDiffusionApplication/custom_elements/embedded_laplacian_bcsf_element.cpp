// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"
#include "utilities/element_size_calculator.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

// Application includes
#include "custom_elements/embedded_laplacian_bcsf_element.h"
#include "convection_diffusion_application_variables.h"

namespace Kratos
{

template<std::size_t TTDim>
EmbeddedLaplacianBCSFElement<TTDim>::EmbeddedLaplacianBCSFElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : LaplacianElement(
        NewId,
        pGeometry)
{
}

template<std::size_t TTDim>
EmbeddedLaplacianBCSFElement<TTDim>::EmbeddedLaplacianBCSFElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : LaplacianElement(
        NewId,
        pGeometry,
        pProperties)
{
}

template<std::size_t TTDim>
Element::Pointer EmbeddedLaplacianBCSFElement<TTDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<EmbeddedLaplacianBCSFElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TTDim>
Element::Pointer EmbeddedLaplacianBCSFElement<TTDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<EmbeddedLaplacianBCSFElement>(NewId, pGeom, pProperties);
}

template<std::size_t TTDim>
EmbeddedLaplacianBCSFElement<TTDim>::~EmbeddedLaplacianBCSFElement()
{
}

template<std::size_t TTDim>
void EmbeddedLaplacianBCSFElement<TTDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    EmbeddedElementData data;
    data.Initialize(*this);

    // Check if the element belongs to the intersected ones
    if (data.IsSplit()) {

        // Get nodal distances and set splitting as well as shape functions
        InitializeGeometryData(data);

        // Resizing and resetting the LHS
        if(rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
            rLeftHandSideMatrix.resize(NumNodes,NumNodes,false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes,NumNodes);

        // Resizing and resetting the RHS
        if(rRightHandSideVector.size() != NumNodes)
            rRightHandSideVector.resize(NumNodes,false);
        noalias(rRightHandSideVector) = ZeroVector(NumNodes);

        // Calculate and add local system for the positive side of the element
        AddPositiveElementSide(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, data);

        // Calculate and add interface terms
        AddPositiveInterfaceTerms(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, data);

        // Calculate and add Nitsche terms for weak imposition of boundary condition
        //AddNitscheBoundaryTerms(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, data);

    } else {
        // Add base Laplacian contribution (standard Galerkin)
        BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

template<std::size_t TTDim>
void EmbeddedLaplacianBCSFElement<TTDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType aux(0);
    CalculateLocalSystem(rLeftHandSideMatrix, aux, rCurrentProcessInfo);
}

template<std::size_t TTDim>
void EmbeddedLaplacianBCSFElement<TTDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType aux(0,0);
    CalculateLocalSystem(aux, rRightHandSideVector, rCurrentProcessInfo);
}

template<std::size_t TTDim>
int EmbeddedLaplacianBCSFElement<TTDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Base Laplacian element check
    return BaseType::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template<std::size_t TTDim>
void EmbeddedLaplacianBCSFElement<TTDim>::InitializeGeometryData(
    EmbeddedElementData& rData)
{
    // Get shape function calculator
    ModifiedShapeFunctions::Pointer p_calculator =
        EmbeddedLaplacianBCSFInternals::GetContinuousShapeFunctionCalculator<TTDim, NumNodes>(
            *this,
            rData.NodalDistances);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN,
        rData.PositiveSideDNDX,
        rData.PositiveSideDNDX_unc,
        rData.PositiveSideWeights,
        this->GetIntegrationMethod());

    // Positive side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN,
        rData.PositiveInterfaceDNDX,
        rData.PositiveInterfaceWeights,
        this->GetIntegrationMethod());

    // Positive side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals,
        this->GetIntegrationMethod());

    // Normalize the interface normals
    const double h = ElementSizeCalculator<TTDim,NumNodes>::AverageElementSize(this->GetGeometry());
    const double tolerance = std::pow(1e-3 * h, TTDim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
}

template<std::size_t TTDim>
void EmbeddedLaplacianBCSFElement<TTDim>::NormalizeInterfaceNormals(
    typename EmbeddedElementData::InterfaceNormalsType& rNormals,
    double Tolerance) const
{
    for (std::size_t i = 0; i < rNormals.size(); ++i) {
        double norm = norm_2(rNormals[i]);
        rNormals[i] /= std::max(norm,Tolerance);
    }
}

template<std::size_t TTDim>
void EmbeddedLaplacianBCSFElement<TTDim>::AddPositiveElementSide(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const EmbeddedElementData& rData)
{
    const auto& r_geom = GetGeometry();
    auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const auto& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    // Get heat flux, conductivity and temp (RHS = ExternalSources - K*temp) nodal vectors
    Vector heat_flux_local(NumNodes);
    Vector nodal_conductivity(NumNodes);
    Vector temp(NumNodes+3);
    for(std::size_t n = 0; n < NumNodes; ++n) {
        heat_flux_local[n] = r_geom[n].FastGetSolutionStepValue(r_volume_source_var);
        nodal_conductivity[n] = r_geom[n].FastGetSolutionStepValue(r_diffusivity_var);
        temp[n] = r_geom[n].GetSolutionStepValue(r_unknown_var);
    }
    const double value_interface_bc = 2.0;
    for(std::size_t i_edge = NumNodes; i_edge < NumNodes+3; ++i_edge) {
        temp[i_edge] = value_interface_bc;
    }

    BoundedMatrix<double, NumNodes, NumNodes+3> LHS_upperPart = ZeroMatrix(NumNodes, NumNodes+3);
    // Iterate over the positive side volume integration points
    // = number of integration points * number of subdivisions on positive side of element
    const std::size_t number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    for (std::size_t g = 0; g < number_of_positive_gauss_points; ++g) {

        const auto& N = row(rData.PositiveSideN, g);
        const auto& DN_DX = rData.PositiveSideDNDX[g];
        const auto& DN_DX_unc = rData.PositiveSideDNDX_unc[g];
        const double weight_gauss = rData.PositiveSideWeights[g];

        //Calculate the local conductivity
        const double conductivity_gauss = inner_prod(N, nodal_conductivity);

        // TODO splice DN_DX_unc?
        BoundedMatrix<double, NumNodes, TTDim> DN_DX_unc_nodalPart;
        for(std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
            for(std::size_t d = 0; d < TTDim; ++d) {
                DN_DX_unc_nodalPart(i_node,d) = DN_DX_unc(i_node,d);
            }
        }
        noalias(rLeftHandSideMatrix) += weight_gauss * conductivity_gauss * prod(DN_DX_unc_nodalPart, trans(DN_DX_unc_nodalPart));
        LHS_upperPart += weight_gauss * conductivity_gauss * prod(DN_DX_unc_nodalPart, trans(DN_DX_unc));

        // Calculate the local RHS (external source)
        const double q_gauss = inner_prod(N, heat_flux_local);

        noalias(rRightHandSideVector) += weight_gauss * q_gauss * N;
    }

    //RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(LHS_upperPart,temp);
}

template<std::size_t TTDim>
void EmbeddedLaplacianBCSFElement<TTDim>::AddPositiveInterfaceTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const EmbeddedElementData& rData)
{
    const auto& r_geom = GetGeometry();
    auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

    // Get conductivity and temperature nodal vectors
    Vector nodal_conductivity(NumNodes);
    Vector temp(NumNodes);
    for(std::size_t n = 0; n < NumNodes; ++n) {
        nodal_conductivity[n] = r_geom[n].FastGetSolutionStepValue(r_diffusivity_var);
        temp[n] = r_geom[n].GetSolutionStepValue(r_unknown_var);
    }

    // Iterate over the positive side interface integration points
    const std::size_t number_of_positive_gauss_points = rData.PositiveInterfaceWeights.size();
    for (std::size_t g = 0; g < number_of_positive_gauss_points; ++g) {

        const auto& N = row(rData.PositiveInterfaceN, g);
        const auto& DN_DX = rData.PositiveInterfaceDNDX[g];
        const double weight_gauss = rData.PositiveInterfaceWeights[g];
        const auto& unit_normal = rData.PositiveInterfaceUnitNormals[g];

        //Calculate the local conductivity
        const double conductivity_gauss = inner_prod(N, nodal_conductivity);

        // Add interface contributions
        for (std::size_t i = 0; i < NumNodes; ++i) {
            for (std::size_t j = 0; j < NumNodes; ++j) {
                for (std::size_t d = 0; d < TTDim; ++d) {

                    const double aux_flux = weight_gauss * conductivity_gauss * N(i) * unit_normal(d) * DN_DX(j,d);

                    rLeftHandSideMatrix(i, j) -= aux_flux;
                    rRightHandSideVector(i) += aux_flux * temp(j);
                }
            }
        }
    }
}

template<std::size_t TTDim>
void EmbeddedLaplacianBCSFElement<TTDim>::AddNitscheBoundaryTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const EmbeddedElementData& rData)
{
    const auto& r_geom = GetGeometry();
    auto& r_settings = *rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    const auto& r_unknown_var = r_settings.GetUnknownVariable();
    const auto& r_diffusivity_var = r_settings.GetDiffusionVariable();

    // Get conductivity and temperature nodal vectors
    Vector nodal_conductivity(NumNodes);
    Vector temp(NumNodes);
    for(std::size_t n = 0; n < NumNodes; ++n) {
        nodal_conductivity[n] = r_geom[n].FastGetSolutionStepValue(r_diffusivity_var);
        temp[n] = r_geom[n].GetSolutionStepValue(r_unknown_var);
    }

    // Nitsche penalty constant
    const double gamma = rCurrentProcessInfo[PENALTY_DIRICHLET];
    // Dirichlet boundary value
    const double temp_bc = GetValue(EMBEDDED_SCALAR);
    // Measure of element size
    const double h = ElementSizeCalculator<TTDim,NumNodes>::AverageElementSize(r_geom);

    // Iterate over the positive side interface integration points
    const std::size_t number_of_positive_gauss_points = rData.PositiveInterfaceWeights.size();
    for (std::size_t g = 0; g < number_of_positive_gauss_points; ++g) {

        const auto& N = row(rData.PositiveInterfaceN, g);
        const auto& DN_DX = rData.PositiveInterfaceDNDX[g];
        const double weight_gauss = rData.PositiveInterfaceWeights[g];
        const auto& unit_normal = rData.PositiveInterfaceUnitNormals[g];

        //Calculate the local conductivity and temperature (at Gauss point)
        const double conductivity_gauss = inner_prod(N, nodal_conductivity);
        const double temp_gauss = inner_prod(N, temp);

        // Add Nitsche contributions
        for (std::size_t i = 0; i < NumNodes; ++i) {

            // Calculate contribution of Nitsche boundary condition - part 1
            const double aux_bc_1 = weight_gauss * conductivity_gauss * gamma / h * N(i);

            // Calculate contribution of Nitsche boundary condition - part 2
            double aux_bc_2 = 0.0;
            for (std::size_t d = 0; d < TTDim; ++d) {
                aux_bc_2 += weight_gauss * conductivity_gauss * DN_DX(i, d) * unit_normal(d);
            }

            // Add contribution of Nitsche boundary condition to LHS
            for (std::size_t j = 0; j < NumNodes; ++j) {
                rLeftHandSideMatrix(i, j) += aux_bc_1 * N(j) - aux_bc_2 * N(j);
            }

            // Add contribution of Nitsche boundary condition to RHS
            rRightHandSideVector(i) -= (aux_bc_1 - aux_bc_2) * (temp_gauss - temp_bc);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace EmbeddedLaplacianBCSFInternals {

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<2, 3>(
    const Element& rElement,
    const Vector& rNodalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Triangle2D3ModifiedShapeFunctions(rElement.pGetGeometry(), rNodalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<3, 4>(
    const Element& rElement,
    const Vector& rNodalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4ModifiedShapeFunctions(rElement.pGetGeometry(), rNodalDistances));
}

template<std::size_t TTDim>
void EmbeddedElementData<TTDim>::Initialize(
    const Element& rElement)
{
    const auto& r_geom = rElement.GetGeometry();;

    // Get nodal distances
    if (NodalDistances.size() != NumNodes) {
        NodalDistances.resize(NumNodes);
    }
    for (std::size_t i = 0; i < NumNodes; ++i) {
        NodalDistances[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
    }

    // Number and indices of positive and negative distance function values
    NumPositiveNodes = 0;
    NumNegativeNodes = 0;

    for (std::size_t i = 0; i < NumNodes; ++i){
        if (NodalDistances[i] > 0.0){
            NumPositiveNodes++;
        } else {
            NumNegativeNodes++;
        }
    }
}

template<std::size_t TTDim>
bool EmbeddedElementData<TTDim>::IsSplit()
{
    return (NumPositiveNodes > 0) && (NumNegativeNodes > 0);
}

} // Namespace EmbeddedLaplacianBCSFInternals

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class EmbeddedLaplacianBCSFElement<2>;
template class EmbeddedLaplacianBCSFElement<3>;

///////////////////////////////////////////////////////////////////////////////////////////////////

} // Namespace Kratos
