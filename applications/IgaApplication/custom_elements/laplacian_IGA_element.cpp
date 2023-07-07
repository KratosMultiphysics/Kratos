// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/laplacian_IGA_element.h"


namespace Kratos
{

template<std::size_t TDim>
LaplacianIGAElement<TDim>::LaplacianIGAElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : LaplacianElement(
        NewId,
        pGeometry)
{
}

template<std::size_t TDim>
LaplacianIGAElement<TDim>::LaplacianIGAElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : LaplacianElement(
        NewId,
        pGeometry,
        pProperties)
{
}

template<std::size_t TDim>
Element::Pointer LaplacianIGAElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianIGAElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim>
Element::Pointer LaplacianIGAElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianIGAElement>(NewId, pGeom, pProperties);
}

template<std::size_t TDim>
LaplacianIGAElement<TDim>::~LaplacianIGAElement()
{
}




template<std::size_t TDim>
void LaplacianIGAElement<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base Laplacian contribution
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template<std::size_t TDim>
void LaplacianIGAElement<TDim>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Add base Laplacian contribution
    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

}

template<std::size_t TDim>
void LaplacianIGAElement<TDim>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Add base Laplacian contribution
    BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

}

template<std::size_t TDim>
int LaplacianIGAElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Surrogate boundary checks

    // Base Laplacian element check
    return BaseType::Check(rCurrentProcessInfo);
}



// // For the Future...
// template<std::size_t TDim>
// std::vector<std::size_t> LaplacianIGAElement<TDim>::GetSurrogateFacesIds()
// {
//     const std::size_t n_faces = TDim + 1;
//     auto& r_neigh_elems = GetValue(NEIGHBOUR_ELEMENTS);
//     // Check the current element faces
//     // Note that we relly on the fact that the neighbours are sorted according to the faces
//     std::vector<std::size_t> surrogate_faces_ids;
//     for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
//         auto p_neigh_elem = r_neigh_elems(i_face).get(); 
//         if (p_neigh_elem != nullptr && p_neigh_elem->Is(VISITED)) {        // BOUNDARY (if it is cut!)
//             surrogate_faces_ids.push_back(i_face);
//         }
//     }
//     // KRATOS_WATCH(surrogate_faces_ids)

//     return surrogate_faces_ids;
// }

template class LaplacianIGAElement<2>;
template class LaplacianIGAElement<3>;

} // Namespace Kratos