#include "includes/kratos_flags.h"

#include "custom_elements/embedded_bcsf_fluid_element.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/weakly_compressible_navier_stokes.h"

#include "utilities/element_size_calculator.h"
#include "custom_utilities/embedded_discontinuous_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/weakly_compressible_navier_stokes_data.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions.h"

namespace Kratos {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TBaseElement >
EmbeddedBCSFFluidElement<TBaseElement>::EmbeddedBCSFFluidElement(IndexType NewId):
    TBaseElement(NewId)
{}

template< class TBaseElement >
EmbeddedBCSFFluidElement<TBaseElement>::EmbeddedBCSFFluidElement(IndexType NewId, const NodesArrayType& ThisNodes):
    TBaseElement(NewId,ThisNodes)
{}


template< class TBaseElement >
EmbeddedBCSFFluidElement<TBaseElement>::EmbeddedBCSFFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry):
    TBaseElement(NewId,pGeometry)
{}

template< class TBaseElement >
EmbeddedBCSFFluidElement<TBaseElement>::EmbeddedBCSFFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties):
    TBaseElement(NewId,pGeometry,pProperties)
{}


template< class TBaseElement >
EmbeddedBCSFFluidElement<TBaseElement>::~EmbeddedBCSFFluidElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TBaseElement >
Element::Pointer EmbeddedBCSFFluidElement<TBaseElement>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<EmbeddedBCSFFluidElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TBaseElement >
Element::Pointer EmbeddedBCSFFluidElement<TBaseElement>::Create(
    IndexType NewId,
    Geometry<NodeType>::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<EmbeddedBCSFFluidElement>(NewId, pGeom, pProperties);
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Call the base element initialize method to set the constitutive law
    TBaseElement::Initialize(rCurrentProcessInfo);

    // Initialize the ELEMENTAL_DISTANCES variable (make it threadsafe)
    if (!this->Has(ELEMENTAL_DISTANCES)) {
        Vector zero_vector(NumNodes, 0.0);
        this->SetValue(ELEMENTAL_DISTANCES, zero_vector);
    }

    // Initialize the nodal EMBEDDED_VELOCITY variable (make it threadsafe)
    const array_1d<double,3> zero_vel = ZeroVector(3);
    for (auto &r_node : this->GetGeometry()) {
        r_node.SetLock();
        if (!r_node.Has(EMBEDDED_VELOCITY)) {
            r_node.SetValue(EMBEDDED_VELOCITY, zero_vel);
        }
        r_node.UnSetLock();
    }

    KRATOS_CATCH("");
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize){
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    if (rRightHandSideVector.size() != LocalSize){
        rRightHandSideVector.resize(LocalSize, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    EmbeddedDiscontinuousElementData data;
    data.Initialize(*this, rCurrentProcessInfo);
    this->InitializeGeometryData(data);

    // Iterate of volume integration points of element or subelements.
    // If there are subelements are, assemble them and apply Dirichlet-BC or convert values for intersection points to values referring to the nodes of the element.
    AddCondensedTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);

    // If the element is cut or incised (using Ausas FE space), add the interface contributions
    if ( data.IsCut() || data.IsIncised() )
    {
        const std::size_t volume_gauss_points = data.PositiveSideWeights.size() + data.NegativeSideWeights.size();

        // Add the base element boundary contribution on the positive interface
        const std::size_t number_of_positive_interface_gauss_points = data.PositiveInterfaceWeights.size();
        for (std::size_t g = 0; g < number_of_positive_interface_gauss_points; ++g){
            const std::size_t gauss_pt_index = g + volume_gauss_points;
            this->UpdateIntegrationPointData(data, gauss_pt_index, data.PositiveInterfaceWeights[g], row(data.PositiveInterfaceN, g), data.PositiveInterfaceDNDX[g]);
            this->AddBoundaryTraction(data, data.PositiveInterfaceUnitNormals[g], rLeftHandSideMatrix, rRightHandSideVector);
        }

        // Add the base element boundary contribution on the negative interface
        const std::size_t number_of_negative_interface_gauss_points = data.NegativeInterfaceWeights.size();
        for (std::size_t g = 0; g < number_of_negative_interface_gauss_points; ++g){
            const std::size_t gauss_pt_index = g + volume_gauss_points + number_of_positive_interface_gauss_points;
            this->UpdateIntegrationPointData(data, gauss_pt_index, data.NegativeInterfaceWeights[g], row(data.NegativeInterfaceN, g), data.NegativeInterfaceDNDX[g]);
            this->AddBoundaryTraction(data, data.NegativeInterfaceUnitNormals[g], rLeftHandSideMatrix, rRightHandSideVector);
        }
    }
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::Calculate(
    const Variable<double> &rVariable,
    double& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = 0.0;
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);

    // If the element is split, integrate sigma.n over the interface
    // Note that in the Ausas formulation, both interface sides need to be integrated
    if (rVariable == DRAG_FORCE) {
        EmbeddedDiscontinuousElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        data.InitializeBoundaryConditionData(rCurrentProcessInfo);
        // Calculate the drag force
        this->CalculateDragForce(data, rOutput);
    } else if (rVariable == DRAG_FORCE_CENTER) {
        EmbeddedDiscontinuousElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        data.InitializeBoundaryConditionData(rCurrentProcessInfo);
        // Calculate the drag force location
        this->CalculateDragForceCenter(data, rOutput);
    } else {
        TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::Calculate(
    const Variable<Vector> &rVariable,
    Vector& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::Calculate(
    const Variable<Matrix> &rVariable,
    Matrix& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <class TBaseElement>
int EmbeddedBCSFFluidElement<TBaseElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int out = EmbeddedDiscontinuousElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    return TBaseElement::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TBaseElement>
const Parameters EmbeddedBCSFFluidElement<TBaseElement>::GetSpecifications() const
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
            "This element implements a Cut-FEM type (a.k.a. embedded) for a discontinuous (element-based) levelset representation. The formulation implemented by this element is specially conceived to work with thin-walled bodies as it is capable to represent the velocity and pressure discontinuities. Note that this element is understood to act as un upper-layer implementing the Cut-FEM terms of a template TBaseElement implementing the Navier-Stokeks contribution. A Navier-Slip boundary condition is imposed in the levelset intersections using the Nitsche's method. The element is able to account for the relative velocity of moving objects by defining the EMBEDDED_VELOCITY variable (this would require switching on the FM-ALE algorithm)."
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
std::string EmbeddedBCSFFluidElement<TBaseElement>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedBCSFFluidElement #" << this->Id();
    return buffer.str();
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedBCSFFluidElement" << Dim << "D" << NumNodes << "N"
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
void EmbeddedBCSFFluidElement<TBaseElement>::AddCondensedTimeIntegratedSystem(
    EmbeddedDiscontinuousElementData& rData,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    //TODO: change N and DNDX used for integration
    // --> convert values from intersection points to nodes and assemble subelements and apply Dirichlet boundary condition

    // Iterate over the positive side volume integration points
    const std::size_t number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    for (std::size_t g = 0; g < number_of_positive_gauss_points; ++g){
        const std::size_t gauss_pt_index = g;
        this->UpdateIntegrationPointData(rData, gauss_pt_index, rData.PositiveSideWeights[g], row(rData.PositiveSideN, g), rData.PositiveSideDNDX[g]);
        this->AddTimeIntegratedSystem(rData, rLeftHandSideMatrix, rRightHandSideVector);
    }

    // Iterate over the negative side volume integration points
    const std::size_t number_of_negative_gauss_points = rData.NegativeSideWeights.size();
    for (std::size_t g = 0; g < number_of_negative_gauss_points; ++g){
        const std::size_t gauss_pt_index = g + number_of_positive_gauss_points;
        this->UpdateIntegrationPointData(rData, gauss_pt_index, rData.NegativeSideWeights[g], row(rData.NegativeSideN, g), rData.NegativeSideDNDX[g]);
        this->AddTimeIntegratedSystem(rData, rLeftHandSideMatrix, rRightHandSideVector);
    }
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::InitializeGeometryData(EmbeddedDiscontinuousElementData& rData) const
{
    rData.PositiveIndices.clear();
    rData.NegativeIndices.clear();

    // Number of positive and negative distance function values
    for (std::size_t i = 0; i < EmbeddedDiscontinuousElementData::NumNodes; ++i){
        if (rData.ElementalDistances[i] > 0.0){
            rData.NumPositiveNodes++;
            rData.PositiveIndices.push_back(i);
        } else {
            rData.NumNegativeNodes++;
            rData.NegativeIndices.push_back(i);
        }
    }

    // Number of edges cut by extrapolated geometry, if not empty
    for (std::size_t i = 0; i < rData.ElementalEdgeDistancesExtrapolated.size(); ++i) {
        if (rData.ElementalEdgeDistancesExtrapolated[i] > 0.0) {
            rData.NumIntersectedEdgesExtrapolated++;
        }
    }

    // Check whether element is intersected or incised, then use Ausas incised shape functions
    if ( rData.IsCut() ) {
        this->DefineCutGeometryData(rData);
    } else if ( rData.IsIncised() ) {
        this->DefineIncisedGeometryData(rData);
    } else {
        this->DefineStandardGeometryData(rData);
    }
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::DefineStandardGeometryData(EmbeddedDiscontinuousElementData& rData) const
{
    rData.NumNegativeNodes = 0;
    rData.NumPositiveNodes = NumNodes;
    this->CalculateGeometryData(rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::DefineCutGeometryData(EmbeddedDiscontinuousElementData& rData) const
{
    // Auxiliary distance vector for the element subdivision utility
    Vector elemental_distances = rData.ElementalDistances;

    ModifiedShapeFunctions::UniquePointer p_calculator =
        EmbeddedBCSFInternals::GetShapeFunctionCalculator<EmbeddedDiscontinuousElementData::Dim, EmbeddedDiscontinuousElementData::NumNodes>(
            *this,
            elemental_distances);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN,
        rData.PositiveSideDNDX,
        rData.PositiveSideWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Negative side volume
    p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeSideN,
        rData.NegativeSideDNDX,
        rData.NegativeSideWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Positive side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN,
        rData.PositiveInterfaceDNDX,
        rData.PositiveInterfaceWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Negative side interface
    p_calculator->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeInterfaceN,
        rData.NegativeInterfaceDNDX,
        rData.NegativeInterfaceWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Positive side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Negative side interface normals
    p_calculator->ComputeNegativeSideInterfaceAreaNormals(
        rData.NegativeInterfaceUnitNormals,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Normalize the normals
    // Note: we calculate h here (and we don't use the value in rData.ElementSize)
    // because rData.ElementSize might still be uninitialized: some data classes define it at the Gauss point.
    double h = ElementSizeCalculator<Dim,NumNodes>::MinimumElementSize(this->GetGeometry());
    const double tolerance = std::pow(1e-3 * h, Dim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
    this->NormalizeInterfaceNormals(rData.NegativeInterfaceUnitNormals, tolerance);
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::DefineIncisedGeometryData(EmbeddedDiscontinuousElementData& rData) const
{
    // Auxiliary distance vector for the element subdivision utility
    Vector elemental_distances = rData.ElementalDistances;
    // Auxiliary edge distance vector of extrapolated intersecting geometry for the element subdivision utility
    Vector edge_distances_extrapolated = rData.ElementalEdgeDistancesExtrapolated;

    ModifiedShapeFunctions::UniquePointer p_calculator =
        EmbeddedBCSFInternals::GetIncisedShapeFunctionCalculator<EmbeddedDiscontinuousElementData::Dim, EmbeddedDiscontinuousElementData::NumNodes>(
            *this,
            elemental_distances,
            edge_distances_extrapolated);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN,
        rData.PositiveSideDNDX,
        rData.PositiveSideWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Negative side volume
    p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeSideN,
        rData.NegativeSideDNDX,
        rData.NegativeSideWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Positive side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN,
        rData.PositiveInterfaceDNDX,
        rData.PositiveInterfaceWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Negative side interface
    p_calculator->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeInterfaceN,
        rData.NegativeInterfaceDNDX,
        rData.NegativeInterfaceWeights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Positive side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Negative side interface normals
    p_calculator->ComputeNegativeSideInterfaceAreaNormals(
        rData.NegativeInterfaceUnitNormals,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Normalize the normals
    // Note: we calculate h here (and we don't use the value in rData.ElementSize)
    // because rData.ElementSize might still be uninitialized: some data classes define it at the Gauss point.
    double h = ElementSizeCalculator<Dim,NumNodes>::MinimumElementSize(this->GetGeometry());
    const double tolerance = std::pow(1e-3 * h, Dim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
    this->NormalizeInterfaceNormals(rData.NegativeInterfaceUnitNormals, tolerance);
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::NormalizeInterfaceNormals(
    typename EmbeddedDiscontinuousElementData::InterfaceNormalsType& rNormals,
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

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::CalculateDragForce(
    EmbeddedDiscontinuousElementData& rData,
    array_1d<double,3>& rDragForce) const
{
    // Initialize the embedded element data
    const std::size_t number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    const std::size_t number_of_negative_gauss_points = rData.NegativeSideWeights.size();
    const std::size_t volume_gauss_points = number_of_positive_gauss_points + number_of_negative_gauss_points;

    if (rData.IsCut()){
        const auto& r_geom = this->GetGeometry();

        // Integrate positive interface side drag
        const std::size_t n_int_pos_gauss = rData.PositiveInterfaceWeights.size();
        for (std::size_t g = 0; g < n_int_pos_gauss; ++g) {
            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points,
                rData.PositiveInterfaceWeights[g],
                row(rData.PositiveInterfaceN, g),
                rData.PositiveInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

            // Compute Gauss pt. values
            const double p_gauss = inner_prod(rData.N, rData.Pressure);
            const array_1d<double, Dim> v_gauss = prod(rData.N, rData.Velocity);
            array_1d<double,Dim> v_emb_gauss = ZeroVector(Dim);
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                const auto &r_i_emb_vel = r_geom[i_node].GetValue(EMBEDDED_VELOCITY);
                for (std::size_t d = 0; d < Dim; ++d) {
                    v_emb_gauss(d) += r_i_emb_vel(d) * rData.N(i_node);
                }
            }

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);
            BoundedMatrix<double, Dim, Dim> norm_proj_matrix, tang_proj_matrix;
            FluidElementUtilities<NumNodes>::SetNormalProjectionMatrix(aux_unit_normal, norm_proj_matrix);
            FluidElementUtilities<NumNodes>::SetTangentialProjectionMatrix(aux_unit_normal, tang_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            const array_1d<double, Dim> shear_proj_n = prod(shear_proj, norm_proj_matrix);
            array_1d<double, Dim> shear_proj_t = ZeroVector(Dim);
            if (rData.SlipLength > 1.0e-12) {
                const auto v_aux = v_gauss - v_emb_gauss;
                const auto v_tan = prod(v_aux, tang_proj_matrix);
                shear_proj_t = rData.Weight * (rData.DynamicViscosity / rData.SlipLength) * v_tan;
            }
            for (std::size_t i = 0; i < Dim ; ++i){
                rDragForce(i) -= shear_proj_n(i);
                rDragForce(i) += shear_proj_t(i);
            }
            rDragForce += rData.Weight * p_gauss * aux_unit_normal;
        }

        // Integrate negative interface side drag
        const std::size_t n_int_neg_gauss = rData.NegativeInterfaceWeights.size();
        for (std::size_t g = 0; g < n_int_neg_gauss; ++g) {
            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points + n_int_pos_gauss,
                rData.NegativeInterfaceWeights[g],
                row(rData.NegativeInterfaceN, g),
                rData.NegativeInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

            // Compute Gauss pt. values
            const double p_gauss = inner_prod(rData.N, rData.Pressure);
            const array_1d<double, Dim> v_gauss = prod(rData.N, rData.Velocity);
            array_1d<double,Dim> v_emb_gauss = ZeroVector(Dim);
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                const auto &r_i_emb_vel = r_geom[i_node].GetValue(EMBEDDED_VELOCITY);
                for (std::size_t d = 0; d < Dim; ++d) {
                    v_emb_gauss(d) += r_i_emb_vel(d) * rData.N(i_node);
                }
            }

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);
            BoundedMatrix<double, Dim, Dim> norm_proj_matrix, tang_proj_matrix;
            FluidElementUtilities<NumNodes>::SetNormalProjectionMatrix(aux_unit_normal, norm_proj_matrix);
            FluidElementUtilities<NumNodes>::SetTangentialProjectionMatrix(aux_unit_normal, tang_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            const array_1d<double, Dim> shear_proj_n = prod(shear_proj, norm_proj_matrix);
            array_1d<double, Dim> shear_proj_t = ZeroVector(Dim);
            if (rData.SlipLength > 1.0e-12) {
                const auto v_aux = v_gauss - v_emb_gauss;
                const auto v_tan = prod(v_aux, tang_proj_matrix);
                shear_proj_t = rData.Weight * (rData.DynamicViscosity / rData.SlipLength) * v_tan;
            }
            for (std::size_t i = 0; i < Dim ; ++i){
                rDragForce(i) -= shear_proj_n(i);
                rDragForce(i) += shear_proj_t(i);
            }
            rDragForce += rData.Weight * p_gauss * aux_unit_normal;
        }
    }
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::CalculateDragForceCenter(
    EmbeddedDiscontinuousElementData& rData,
    array_1d<double,3>& rDragForceLocation) const
{
    const auto &r_geometry = this->GetGeometry();
    array_1d<double,3> tot_drag = ZeroVector(3);
    const std::size_t number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    const std::size_t number_of_negative_gauss_points = rData.NegativeSideWeights.size();
    const std::size_t volume_gauss_points = number_of_positive_gauss_points + number_of_negative_gauss_points;

    if (rData.IsCut()){
        // Get the positive interface continuous shape functions
        // We use these ones to interpolate the position of the intersection Gauss pt.
        // Note that we take advantage of the fact that the positive and negative interface Gauss pt. coincide
        Vector pos_int_continuous_weights;
        Matrix pos_int_continuous_N;
        typename EmbeddedDiscontinuousElementData::ShapeFunctionsGradientsType pos_int_continuous_DN_DX;
        auto p_continuous_sh_func_calculator = EmbeddedBCSFInternals::GetContinuousShapeFunctionCalculator<Dim, NumNodes>(*this, rData.ElementalDistances);
        p_continuous_sh_func_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
            pos_int_continuous_N,
            pos_int_continuous_DN_DX,
            pos_int_continuous_weights,
            GeometryData::IntegrationMethod::GI_GAUSS_2);

        // Integrate positive interface side drag
        const std::size_t n_int_pos_gauss = rData.PositiveInterfaceWeights.size();
        for (std::size_t g = 0; g < n_int_pos_gauss; ++g) {
            // Obtain the Gauss pt. coordinates using the standard shape functions
            array_1d<double,3> g_coords = ZeroVector(3);
            const auto g_shape_functions = row(pos_int_continuous_N, g);
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                g_coords += g_shape_functions[i_node] * r_geometry[i_node].Coordinates();
            }

            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points,
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
            const array_1d<double, 3> p_proj = rData.Weight * p_gauss * aux_unit_normal;
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (std::size_t i = 0; i < Dim ; ++i){
                tot_drag(i) -= shear_proj(i);
                rDragForceLocation(i) += g_coords(i) * p_proj(i);
                rDragForceLocation(i) -= g_coords(i) * shear_proj(i);
            }
            tot_drag += p_proj;
        }

        // Integrate negative interface side drag
        const std::size_t n_int_neg_gauss = rData.NegativeInterfaceWeights.size();
        for (std::size_t g = 0; g < n_int_neg_gauss; ++g) {
            // Obtain the Gauss pt. coordinates using the standard shape functions
            array_1d<double,3> g_coords = ZeroVector(3);
            const auto g_shape_functions = row(pos_int_continuous_N, g);
            for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                g_coords += g_shape_functions[i_node] * r_geometry[i_node].Coordinates();
            }

            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points + n_int_pos_gauss,
                rData.NegativeInterfaceWeights[g],
                row(rData.NegativeInterfaceN, g),
                rData.NegativeInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, 3> p_proj = rData.Weight * p_gauss * aux_unit_normal;
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (std::size_t i = 0; i < Dim ; ++i){
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

template <class TBaseElement>
double EmbeddedBCSFFluidElement<TBaseElement>::AuxiliaryDensityGetter(
    const EmbeddedDiscontinuousElementData& rData,
    const std::size_t NodeIndex) const
{
    return rData.Density;
}

template <>
double EmbeddedBCSFFluidElement<WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> >>::AuxiliaryDensityGetter(
    const EmbeddedDiscontinuousElementData& rData,
    const std::size_t NodeIndex) const
{
    return rData.Density(NodeIndex);
}

template <>
double EmbeddedBCSFFluidElement<WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<3,4> >>::AuxiliaryDensityGetter(
    const EmbeddedDiscontinuousElementData& rData,
    const std::size_t NodeIndex) const
{
    return rData.Density(NodeIndex);
}

// serializer

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement);
}

template <class TBaseElement>
void EmbeddedBCSFFluidElement<TBaseElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace EmbeddedBCSFInternals {

template <>
ModifiedShapeFunctions::UniquePointer GetShapeFunctionCalculator<2, 3>(const Element& rElement, const Vector& rElementalDistances)
{
    return Kratos::make_unique<Triangle2D3AusasModifiedShapeFunctions>(
        rElement.pGetGeometry(),
        rElementalDistances);
}

template <>
ModifiedShapeFunctions::UniquePointer GetShapeFunctionCalculator<3, 4>(const Element& rElement, const Vector& rElementalDistances)
{
    return Kratos::make_unique<Tetrahedra3D4AusasModifiedShapeFunctions>(
        rElement.pGetGeometry(),
        rElementalDistances);
}

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<2, 3>(
    const Element& rElement,
    const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Triangle2D3ModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<3, 4>(
    const Element& rElement,
    const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4ModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::UniquePointer GetIncisedShapeFunctionCalculator<2, 3>(
    const Element& rElement,
    const Vector& rElementalDistancesWithExtrapolated,
    const Vector& rElementalEdgeDistancesExtrapolated)
{
    return Kratos::make_unique<Triangle2D3AusasIncisedShapeFunctions>(
        rElement.pGetGeometry(),
        rElementalDistancesWithExtrapolated,
        rElementalEdgeDistancesExtrapolated);
}

template <>
ModifiedShapeFunctions::UniquePointer GetIncisedShapeFunctionCalculator<3, 4>(
    const Element& rElement,
    const Vector& rElementalDistancesWithExtrapolated,
    const Vector& rElementalEdgeDistancesExtrapolated)
{
    return Kratos::make_unique<Tetrahedra3D4AusasIncisedShapeFunctions>(
        rElement.pGetGeometry(),
        rElementalDistancesWithExtrapolated,
        rElementalEdgeDistancesExtrapolated);
}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class EmbeddedBCSFFluidElement< WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> > >;
template class EmbeddedBCSFFluidElement< WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}
