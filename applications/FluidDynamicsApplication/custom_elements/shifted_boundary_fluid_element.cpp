// Project includes
#include "includes/kratos_flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "utilities/geometry_utilities.h"
#include "utilities/element_size_calculator.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

// Application includes
#include "custom_elements/shifted_boundary_fluid_element.h"
#include "custom_elements/weakly_compressible_navier_stokes.h"
#include "custom_utilities/embedded_data.h"
#include "custom_utilities/weakly_compressible_navier_stokes_data.h"
#include <cstddef>
#include <string>


namespace Kratos {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::ShiftedBoundaryFluidElement(IndexType NewId):
    TBaseElement(NewId)
{}

template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::ShiftedBoundaryFluidElement(IndexType NewId, const NodesArrayType& ThisNodes):
    TBaseElement(NewId,ThisNodes)
{}


template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::ShiftedBoundaryFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry):
    TBaseElement(NewId,pGeometry)
{}

template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::ShiftedBoundaryFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties):
    TBaseElement(NewId,pGeometry,pProperties)
{}


template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::~ShiftedBoundaryFluidElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TBaseElement >
Element::Pointer ShiftedBoundaryFluidElement<TBaseElement>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ShiftedBoundaryFluidElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TBaseElement >
Element::Pointer ShiftedBoundaryFluidElement<TBaseElement>::Create(
    IndexType NewId,
    Geometry<NodeType>::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ShiftedBoundaryFluidElement>(NewId, pGeom, pProperties);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Call the base element initialize method to set the constitutive law
    TBaseElement::Initialize(rCurrentProcessInfo);

    // Initialize the ELEMENTAL_DISTANCES variable (make it threadsafe)
    if (!this->Has(ELEMENTAL_DISTANCES)) {
        Vector zero_vector(NumNodes, 0.0);
        this->SetValue(ELEMENTAL_DISTANCES, zero_vector);
    }

    KRATOS_CATCH("");
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base fluid contribution (volume integration)
    TBaseElement::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface.
    // Note that the INTERFACE flag is assumed to be set in the 1st layer of elements attached to the surrogate boundary e.g. by the ShiftedBoundaryMeshlessInterfaceUtility.
    // At the faces of these INTERFACE elements which are attached to BOUNDARY elements (surrogate interface gamma_tilde), the boundary flux contributions is added as a surrogate.
    if (this->Is(INTERFACE)) {

        // Initialize the element data
        ShiftedBoundaryElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);

        // Get the surrogate faces local IDs.
        // Note that it might happen that an INTERFACE element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();
        if (sur_bd_ids_vect.size() != 0) {

            // Get the parent geometry data
            double size_parent;
            const auto& r_parent_geom = this->GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, Dim> DN_DX_parent;  // a node's shape function derivatives at the respective node
            GeometryUtils::CalculateGeometryData(r_parent_geom, DN_DX_parent, N_parent, size_parent);
            const auto& r_boundaries = r_parent_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_parent_geom.NodesInFaces(nodes_in_faces);

            // Loop the surrogate faces of the element
            // NOTE that there is the chance that the surrogate face is not unique
            const std::size_t n_volume_gauss_points = data.PositiveSideWeights.size();
            std::size_t surrogate_pt_index = 0;
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_bd_geom = r_boundaries[sur_bd_id];
                //const unsigned int n_bd_points = r_bd_geom.PointsNumber();  //NOTE 2 nodes for a triangle parent (EdgeType is Line2D2N)
                const DenseVector<std::size_t> n_bd_local_ids = column(nodes_in_faces, sur_bd_id);  //ERROR?? face is column not row, for triangle row works as well because of symmetry, first node is the contrary one for tri and tetra
                //auto& r_sur_bd_N = r_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);   //NOTE returns matrix of SF values F_{ij}, where i is integration pt index and j is SF index
                const auto& r_integration_points = r_bd_geom.IntegrationPoints(r_bd_geom.GetDefaultIntegrationMethod());
                //const std::size_t n_int_pts = r_bd_geom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_1); // Number of Gauss pts. per subdivision

                // Get the gradient of the node contrary to the surrogate face
                // NOTE that this is used to calculate the normal as n = - DN_DX_contrary_node / norm_2(DN_DX_contrary_node)
                // NOTE that DN_DX of a node is calculated as the normal of the opposing face (cross product for Tetrahedra3D4N)
                BoundedVector<double,Dim> DN_DX_contrary_node = row(DN_DX_parent, n_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(DN_DX_contrary_node);
                BoundedVector<double,Dim> normal_sur_bd = - DN_DX_contrary_node * h_sur_bd;

                // Calculate a scaling factor to eliminate the change of area from global to local coordinate system
                /*double scaling_factor = 0.0;
                for (auto& r_int_pt : r_integration_points) {
                    scaling_factor += r_int_pt.Weight();
                }*/

                // NOTE that the integration weight is calculated as detJ of the parent multiplied by the global size of the surrogate boundary face: Dim * Parent domain volume * norm(DN_DX_contrary_node)
                // NOTE that detJ is only constant inside the element for simplex elements and only for Triangle2D3N and Tetrahedra3D4N Dim*size_parent results in the correct multiplier for the norm
                // Triangle2D3N: 2*1/2*detJ*length_of_line | Tetrahedra3D4N: 3*1/6*detJ*area_of_parallelogram
                //double weight = Dim * size_parent / h_sur_bd / scaling_factor;
                // NOTE that the integration weight is calculated as the global size of the surrogate boundary face: Triangle2D3N: length_of_line | Tetrahedra3D4N: 1/2*area_of_parallelogram
                double weight = (Dim==2) ? 1.0/h_sur_bd : 0.5/h_sur_bd;

                // Get detJ for all integration points of the surrogate boundary face
                Vector int_pt_detJs;
                r_bd_geom.DeterminantOfJacobian(int_pt_detJs, r_bd_geom.GetDefaultIntegrationMethod());

                // Loop over the integration points of the surrogate boundary face for the numerical integration of the surrogate boundary flux
                for (std::size_t i_int_pt = 0; i_int_pt < r_integration_points.size(); i_int_pt++) {
                    // Integration weight of one integration point is calculated by multiplying the scaled weight by the respective integration point weight
                    //const double int_pt_weight = weight * r_int_pt.Weight();
                    const double int_pt_weight = weight * int_pt_detJs[i_int_pt] * r_integration_points[i_int_pt].Weight();

                    // Compute the local coordinates of the integration point in the parent element's geometry
                    Geometry<Node>::CoordinatesArrayType aux_global_coords = ZeroVector(3);
                    r_bd_geom.GlobalCoordinates(aux_global_coords, i_int_pt);
                    Geometry<Node>::CoordinatesArrayType int_pt_local_coords_parent = ZeroVector(3);
                    r_parent_geom.PointLocalCoordinates(int_pt_local_coords_parent, aux_global_coords);

                    // Get N of the element at the integration point
                    Matrix int_pt_N_parent = ZeroMatrix(1, NumNodes);
                    Vector aux_N_parent;
                    r_parent_geom.ShapeFunctionsValues(aux_N_parent, int_pt_local_coords_parent);
                    for (std::size_t i_node = 0; i_node < NumNodes; i_node++) {
                        int_pt_N_parent(0, i_node) = aux_N_parent(i_node);
                    }

                    // Get DN_DX of the element at the integration point
                    Matrix int_pt_DN_DX_parent = ZeroMatrix(NumNodes, Dim), aux_DN_DXi_parent, aux_J_parent, aux_J_inv_parent;
                    double aux_detJ_parent;
                    r_parent_geom.ShapeFunctionsLocalGradients(aux_DN_DXi_parent, int_pt_local_coords_parent);
                    r_parent_geom.Jacobian(aux_J_parent, int_pt_local_coords_parent);
                    MathUtils<double>::InvertMatrix(aux_J_parent, aux_J_inv_parent, aux_detJ_parent);
                    int_pt_DN_DX_parent = prod(aux_DN_DXi_parent, aux_J_inv_parent);

                    // Update the element's data
                    this->UpdateIntegrationPointData(data, n_volume_gauss_points+surrogate_pt_index++, int_pt_weight, row(int_pt_N_parent,0), int_pt_DN_DX_parent);

                    // Calculate the surrogate boundary traction as t_i = tau_ij - p_h * I_ij) * n_j with tau_ij = C:delta^s u_h taking N, DN_DX and Weight from the element's data
                    this->AddBoundaryTraction(data, normal_sur_bd, rLeftHandSideMatrix, rRightHandSideVector);
                }
            }
        }
    }

    KRATOS_CATCH("")
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType rRightHandSideVector;
    this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType rLeftHandSideMatrix;
    this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<double> &rVariable,
    double& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == CUTTED_AREA) {
        // Initialize the embedded element data
        ShiftedBoundaryElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        // Calculate the intersection area as the Gauss weights summation
        const unsigned int n_int_pos_gauss = data.PositiveInterfaceWeights.size();
        rOutput = 0.0;
        for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {
            rOutput += data.PositiveInterfaceWeights[g];
        }
    } else {
        TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);

    // If the element is split, integrate sigma.n over the interface
    // Note that in the Ausas formulation (discontinuous for thin-walled), both interface sides need to be integrated
    if (rVariable == DRAG_FORCE) {
        ShiftedBoundaryElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        // Calculate the drag force
        this->CalculateDragForce(data, rOutput);
    } else if (rVariable == DRAG_FORCE_CENTER) {
        ShiftedBoundaryElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        // Calculate the drag force location
        this->CalculateDragForceCenter(data, rOutput);
    } else {
        TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<Vector> &rVariable,
    Vector& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<Matrix> &rVariable,
    Matrix& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <class TBaseElement>
int ShiftedBoundaryFluidElement<TBaseElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int out = ShiftedBoundaryElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    return TBaseElement::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TBaseElement>
const Parameters ShiftedBoundaryFluidElement<TBaseElement>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["VELOCITY","PRESSURE"],
            "nodal_non_historical"   : ["EMBEDDED_VELOCITY"],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","PRESSURE","MESH_VELOCITY","MESH_DISPLACEMENT"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Tetrahedra3D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["Newtonian2DLaw","Newtonian3DLaw","NewtonianTemperatureDependent2DLaw","NewtonianTemperatureDependent3DLaw","Euler2DLaw","Euler3DLaw"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element is based on the Shifted-Boundary Method using MLS shape functions. Therfore, it adds surrogate boundary flux contributions for elements marked with the INTERFACE flag. The element is meant to be used together with the NavierStokesShiftedBoundarySolver, which sets up the ShiftedBoundaryMeshlessInterfaceUtility creating ShiftedBoundaryWallConditions for the interface."
    })");

    if (Dim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template <class TBaseElement>
std::string ShiftedBoundaryFluidElement<TBaseElement>::Info() const
{
    std::stringstream buffer;
    buffer << "ShiftedBoundaryFluidElement #" << this->Id();
    return buffer.str();
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ShiftedBoundaryFluidElement" << Dim << "D" << NumNodes << "N"
             << std::endl
             << "on top of ";
    TBaseElement::PrintInfo(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// Operations

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::InitializeGeometryData(ShiftedBoundaryElementData& rData) const
{
    rData.PositiveIndices.clear();
    rData.NegativeIndices.clear();

    // Number of positive and negative distance function values
    for (std::size_t i = 0; i < ShiftedBoundaryElementData::NumNodes; ++i){
        if (rData.Distance[i] > 0.0) {
            rData.NumPositiveNodes++;
            rData.PositiveIndices.push_back(i);
        }
        else {
            rData.NumNegativeNodes++;
            rData.NegativeIndices.push_back(i);
        }
    }

    if ( rData.IsCut() ) {
        this->DefineCutGeometryData(rData);
    }
    else {
        this->DefineStandardGeometryData(rData);
    }
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::DefineStandardGeometryData(ShiftedBoundaryElementData& rData) const
{
    rData.NumNegativeNodes = 0;
    rData.NumPositiveNodes = NumNodes;
    this->CalculateGeometryData(rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::DefineCutGeometryData(ShiftedBoundaryElementData& rData) const
{
    // Auxiliary distance vector for the element subdivision utility
    Vector distances = rData.Distance;

    ModifiedShapeFunctions::Pointer p_calculator = ShiftedBoundaryInternals::GetShapeFunctionCalculator<
        ShiftedBoundaryElementData::Dim,
        ShiftedBoundaryElementData::NumNodes>(*this, distances);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN, rData.PositiveSideDNDX, rData.PositiveSideWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Positive side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN, rData.PositiveInterfaceDNDX, rData.PositiveInterfaceWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Positive side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Normalize the normals
    // Note: we calculate h here (and we don't use the value in rData.ElementSize)
    // because rData.ElementSize might still be uninitialized: some data classes define it at the Gauss point.
    double h = ElementSizeCalculator<Dim,NumNodes>::MinimumElementSize(this->GetGeometry());
    const double tolerance = std::pow(1e-3 * h, Dim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::NormalizeInterfaceNormals(
    typename ShiftedBoundaryElementData::InterfaceNormalsType& rNormals,
    double Tolerance) const
{
    for (std::size_t i = 0; i < rNormals.size(); ++i) {
        double norm = norm_2(rNormals[i]);
        rNormals[i] /= std::max(norm,Tolerance);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template<class TBaseElement>
std::vector<std::size_t> ShiftedBoundaryFluidElement<TBaseElement>::GetSurrogateFacesIds()
{
    const std::size_t n_faces = Dim + 1;
    auto& r_neigh_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    // Check the current element faces
    // Note that we rely on the fact that the neighbors are sorted according to the faces
    std::vector<std::size_t> surrogate_faces_ids;
    for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
        auto p_neigh_elem = r_neigh_elems(i_face).get();
        if (p_neigh_elem != nullptr && p_neigh_elem->Is(BOUNDARY)) {
            surrogate_faces_ids.push_back(i_face);
        }
    }

    return surrogate_faces_ids;
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::CalculateDragForce(
    ShiftedBoundaryElementData& rData,
    array_1d<double,3>& rDragForce) const
{
    const std::size_t n_pos_gauss = rData.PositiveSideWeights.size();

    if (rData.IsCut()) {
        // Integrate positive interface side drag
        const unsigned int n_int_pos_gauss = rData.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {

            // Update the Gauss pt. rData and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + n_pos_gauss,
                rData.PositiveInterfaceWeights[g],
                row(rData.PositiveInterfaceN, g),
                rData.PositiveInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i) {
                rDragForce(i) -= shear_proj(i);
            }
            rDragForce += rData.Weight * p_gauss * aux_unit_normal;
        }
    }
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::CalculateDragForceCenter(
    ShiftedBoundaryElementData& rData,
    array_1d<double,3>& rDragForceLocation) const
{
    const auto &r_geometry = this->GetGeometry();
    array_1d<double,3> tot_drag = ZeroVector(3);
    const unsigned int n_pos_gauss = rData.PositiveSideWeights.size();

    if (rData.IsCut()) {
        // Integrate positive interface side drag
        const unsigned int n_int_pos_gauss = rData.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {
            // Calculate the Gauss pt. coordinates
            array_1d<double,3> g_coords = ZeroVector(3);
            const auto g_shape_functions = row(rData.PositiveInterfaceN, g);
            for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
                g_coords += g_shape_functions[i_node] * r_geometry[i_node].Coordinates();
            }

            // Update the Gauss pt. rData and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + n_pos_gauss,
                rData.PositiveInterfaceWeights[g],
                g_shape_functions,
                rData.PositiveInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, 3> p_proj = rData.Weight * p_gauss * aux_unit_normal;
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i) {
                tot_drag(i) -= shear_proj(i);
                rDragForceLocation(i) += g_coords(i) * p_proj(i);
                rDragForceLocation(i) -= g_coords(i) * shear_proj(i);
            }
            tot_drag += p_proj;
        }

        // Divide the obtained result by the total drag
        rDragForceLocation(0) /= tot_drag(0);
        rDragForceLocation(1) /= tot_drag(1);
        if (Dim == 3) {
            rDragForceLocation(2) /= tot_drag(2);
        }
    }
}

// serializer

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace ShiftedBoundaryInternals {

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<2, 3>(
    const Element& rElement, const Vector& rDistance) {
    return ModifiedShapeFunctions::Pointer(new Triangle2D3ModifiedShapeFunctions(rElement.pGetGeometry(),rDistance));
}

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<3, 4>(
    const Element& rElement, const Vector& rDistance) {
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4ModifiedShapeFunctions(rElement.pGetGeometry(),rDistance));
}

}  // namespace ShiftedBoundaryInternals

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class ShiftedBoundaryFluidElement< WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> > >;
template class ShiftedBoundaryFluidElement< WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace Kratos
