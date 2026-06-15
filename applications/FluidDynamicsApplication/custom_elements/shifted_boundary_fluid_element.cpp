// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "input_output/logger.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/shifted_boundary_fluid_element.h"
#include "custom_elements/weakly_compressible_navier_stokes.h"
#include "custom_elements/data_containers/weakly_compressible_navier_stokes/weakly_compressible_navier_stokes_data.h"
#include <cstddef>
#include <ostream>
#include <string>
#include <sys/types.h>


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
    // Call the base element initialize method to set the constitutive law
    TBaseElement::Initialize(rCurrentProcessInfo);
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
    // Note that the SBM_INTERFACE flag is assumed to be set in the 1st layer of elements attached to the surrogate boundary (BOUNDARY elements) e.g. by the ShiftedBoundaryUtility.
    // One or multiple faces of these SBM_INTERFACE elements are attached to SBM_BOUNDARY elements.
    // These faces are therefore a part of the surrogate interface gamma_tilde, so the boundary flux contributions is added.
    // Integration is done over the boundary using the faces' geometry and its GI_GAUSS_2 integration points.
    // NOTE that this is the safer alternative because integration points are used directly and do not need to be found inside the element.
    if (this->Is(SBM_INTERFACE)) {
        // Initialize the element data
        ShiftedBoundaryElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Get the local IDs of the SBM_INTERFACE elements of the faces attached to SBM_BOUNDARY elements.
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();

        // Set integration method for surrogate boundary face
        // NOTE that second order is required here
        const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;

        // Get the parent geometry data. It is calculated at the element midpoint.
        const auto& r_geom = this->GetGeometry();
        double size_parent;
        array_1d<double, NumNodes> N_parent;
        BoundedMatrix<double, NumNodes, Dim> DN_DX_parent;
        GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, size_parent);

        // Calculate material response in order to update C matrix (constant)
        this->CalculateMaterialResponse(data);

        // Get strain matrix
        // NOTE that here we are assuming constant strain kinematics (which is correct for the linear shape functions of simplex elements)
        BoundedMatrix<double, StrainSize, LocalSize> B_matrix;
        FluidElementUtilities<NumNodes>::GetStrainMatrix(DN_DX_parent, B_matrix);

        // Generate boundary topology
        const auto &r_boundaries = r_geom.GenerateBoundariesEntities();
        DenseMatrix<unsigned int> nodes_in_faces;
        r_geom.NodesInFaces(nodes_in_faces);

        // Auxilary LHS matrix for summation over the element's surrogate boundaries
        BoundedMatrix<double, LocalSize, LocalSize> aux_LHS = ZeroMatrix(LocalSize, LocalSize);

        // Loop the surrogate faces
        // NOTE that more than one face of the element might be part of the surrogate boundary
        for (std::size_t sur_bd_id : sur_bd_ids_vect) {

            const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
            const std::size_t n_bd_points = r_sur_bd_geom.PointsNumber();  // number of nodes of the surrogate face
            const DenseVector<std::size_t> sur_bd_local_ids = column(nodes_in_faces, sur_bd_id);

            // Get integration points and their global shape function values of the surrogate boundary face
            const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(integration_method);  // matrix of SF values F_{ij}, where i is integration pt index and j is SF index
            const std::vector<IntegrationPoint<3>> & r_integration_points = r_sur_bd_geom.IntegrationPoints(integration_method);

            // Get values of determinant of Jacobian for all integration points of the surrogate boundary face
            VectorType sur_bd_detJ;
            r_sur_bd_geom.DeterminantOfJacobian(sur_bd_detJ, integration_method);

            // Get the gradient of the node opposite of the surrogate face
            // NOTE that the gradient is used to calculate the normal as n = - DN_DX_opposite_node / norm_2(DN_DX_opposite_node)
            const BoundedVector<double,Dim> DN_DX_opposite_node = row(DN_DX_parent, sur_bd_local_ids[0]);
            const double h_sur_bd = 1.0 / norm_2(DN_DX_opposite_node);
            BoundedVector<double,3> normal_sur_bd = ZeroVector(3);
            for (std::size_t d = 0; d < Dim; ++d) {
                normal_sur_bd(d) = - DN_DX_opposite_node(d) * h_sur_bd;
            }

            // Compute the required projections using the boundary's normal,
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_projection_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(normal_sur_bd, voigt_normal_projection_matrix);

            // Compute auxilary matrices for the shear stress with the elements constitutive matrix (C) and the strain matrix (B)
            const BoundedMatrix<double, Dim, StrainSize> AC = prod(voigt_normal_projection_matrix, data.C);
            const BoundedMatrix<double, Dim, LocalSize> ACB = prod(AC, B_matrix);

            // Precompute node global DOF row bases for boundary nodes once per face
            // sur_bd_local_ids[0] is the opposite node; [1..n_bd_points] are the nodes of the face
            std::array<std::size_t, NumNodes-1> bd_node_base;  // max face nodes = NumNodes-1
            for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; ++i_bd_node) {
                bd_node_base[i_bd_node] = sur_bd_local_ids[i_bd_node + 1] * BlockSize;
            }

            // Loop over the integration points of the surrogate boundary face for the numerical integration of the surrogate boundary flux
            for (std::size_t i_int_pt = 0; i_int_pt < r_integration_points.size(); ++i_int_pt) {

                // Calculate integration point weight by multiplying its detJ with its gauss weight
                const double int_pt_weight = sur_bd_detJ[i_int_pt] * r_integration_points[i_int_pt].Weight();

                for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; ++i_bd_node) {
                    const double w_Ni = int_pt_weight * r_sur_bd_N(i_int_pt, i_bd_node);
                    const std::size_t base_i = bd_node_base[i_bd_node];

                    // --- Shear stress contribution ---
                    // aux_LHS -= w * N_u^T * ACB
                    for (std::size_t di = 0; di < Dim; ++di) {
                        const std::size_t row = base_i + di;
                        for (std::size_t col = 0; col < LocalSize; ++col) {
                            aux_LHS(row, col) -= w_Ni * ACB(di, col);
                        }
                    }

                    // --- Pressure contribution ---
                    // aux_LHS += w * N_u^T * n_d * N_p
                    for (std::size_t di = 0; di < Dim; ++di) {
                        const std::size_t row = base_i + di;
                        const double w_Ni_nd = w_Ni * normal_sur_bd(di);
                        for (std::size_t j_bd_node = 0; j_bd_node < n_bd_points; ++j_bd_node) {
                            const std::size_t pres_col = bd_node_base[j_bd_node] + Dim;
                            aux_LHS(row, pres_col) += w_Ni_nd * r_sur_bd_N(i_int_pt, j_bd_node);
                        }
                    }
                }
            }
        }

        // Add boundary traction of the element's surrogate boundaries to the system
        array_1d<double,LocalSize> values;
        this->GetCurrentValuesVector(data, values);
        rLeftHandSideMatrix  += aux_LHS;
        rRightHandSideVector -= prod(aux_LHS, values);
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
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<VectorType> &rVariable,
    VectorType& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<MatrixType> &rVariable,
    MatrixType& rOutput,
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
        "required_variables"         : ["VELOCITY","PRESSURE","MESH_VELOCITY","MESH_DISPLACEMENT","DISTANCE"],
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
    //NOTE necessary for (dis-)continuous level set with EmbeddedData
    /*rData.PositiveIndices.clear();
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
    rData.NumNegativeNodes = 0;
    rData.NumPositiveNodes = NumNodes;
    this->CalculateGeometryData(rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);*/
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template<class TBaseElement>
std::vector<std::size_t> ShiftedBoundaryFluidElement<TBaseElement>::GetSurrogateFacesIds()
{
    const std::size_t n_faces = Dim + 1;  // NOTE this is only valid for tri and tetra
    auto& r_neigh_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    // Check the current element faces
    // NOTE that we rely on the fact that the neighbors are sorted according to the order of faces in NEIGHBOUR_ELEMENTS
    std::vector<std::size_t> surrogate_faces_ids;
    for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
        auto p_neigh_elem = r_neigh_elems(i_face).get();
        if (p_neigh_elem != nullptr && p_neigh_elem->Is(SBM_BOUNDARY)) {
            surrogate_faces_ids.push_back(i_face);
        }
    }

    return surrogate_faces_ids;
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
// Class template instantiation

template class ShiftedBoundaryFluidElement< WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> > >;
template class ShiftedBoundaryFluidElement< WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace Kratos
