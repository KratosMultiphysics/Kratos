#include "custom_elements/embedded_fluid_element.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/symbolic_navier_stokes.h"

#include "custom_utilities/embedded_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos {

template class EmbeddedFluidElement< QSVMS< TimeIntegratedQSVMSData<3,4> > >;
template class EmbeddedFluidElement< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::EmbeddedFluidElement(IndexType NewId):
    TBaseElement(NewId)
{}

template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::EmbeddedFluidElement(IndexType NewId, const NodesArrayType& ThisNodes):
    TBaseElement(NewId,ThisNodes)
{}


template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::EmbeddedFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry):
    TBaseElement(NewId,pGeometry)
{}

template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::EmbeddedFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties):
    TBaseElement(NewId,pGeometry,pProperties)
{}


template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::~EmbeddedFluidElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TBaseElement >
Element::Pointer EmbeddedFluidElement<TBaseElement>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Element::Pointer(new EmbeddedFluidElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}


template< class TBaseElement >
Element::Pointer EmbeddedFluidElement<TBaseElement>::Create(IndexType NewId,Geometry<NodeType>::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Element::Pointer(new EmbeddedFluidElement(NewId, pGeom, pProperties));
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {

    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    EmbeddedElementData data;
    data.Initialize(*this, rCurrentProcessInfo);
    this->InitializeGeometryData(data);

    // Iterate over integration points on the volume
    const unsigned int number_of_positive_gauss_points =
        data.PositiveSideWeights.size();
    for (unsigned int g = 0; g < number_of_positive_gauss_points; g++) {
        data.UpdateGeometryValues(data.PositiveSideWeights[g],
            row(data.PositiveSideN, g), data.PositiveSideDNDX[g]);

        this->AddTimeIntegratedSystem(
            data, rLeftHandSideMatrix, rRightHandSideVector);
    }

    // Iterate over integration points on the boundary
    const unsigned int number_of_interface_gauss_points = data.PositiveInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_interface_gauss_points; g++) {
        data.UpdateGeometryValues(data.PositiveInterfaceWeights[g], row(data.PositiveInterfaceN,g), data.PositiveSideDNDX[g]);

        // Set the normal too?
        // Actually add the terms
    }

    // Add terms specific of this embedded formulation

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <class TBaseElement>
int EmbeddedFluidElement<TBaseElement>::Check(
    const ProcessInfo& rCurrentProcessInfo) {

    EmbeddedElementData::Check(*this,rCurrentProcessInfo);
    return TBaseElement::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TBaseElement>
std::string EmbeddedFluidElement<TBaseElement>::Info() const {
    std::stringstream buffer;
    buffer << "EmbeddedFluidElement #" << this->Id();
    return buffer.str();
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::PrintInfo(
    std::ostream& rOStream) const {
    rOStream << "EmbeddedFluidElement" << Dim << "D" << NumNodes << "N"
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
void EmbeddedFluidElement<TBaseElement>::InitializeGeometryData(
    EmbeddedElementData& rData) const {

    rData.PositiveIndices.clear();
    rData.NegativeIndices.clear();

    // Number of positive and negative distance function values
    for (size_t i = 0; i < EmbeddedElementData::NumNodes; ++i) {

        if (rData.Distance[i] > 0.0) {
            rData.NumPositiveNodes++;
            rData.PositiveIndices.push_back(i);
        }
        else {
            rData.NumNegativeNodes++;
            rData.NegativeIndices.push_back(i);
        }
    }

    if ( (rData.NumPositiveNodes > 0) && (rData.NumNegativeNodes > 0) ) {
        this->DefineCutGeometryData(rData);
    }
    else {
        this->DefineStandardGeometryData(rData);
    }
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::DefineStandardGeometryData(
    EmbeddedElementData& rData) const {

    this->CalculateGeometryData(
        rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);
    rData.NumPositiveNodes = NumNodes;
    rData.NumNegativeNodes = 0;
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::DefineCutGeometryData(
    EmbeddedElementData& rData) const {

    // Auxiliary distance vector for the element subdivision utility
    Vector distances = rData.Distance;

    ModifiedShapeFunctions::Pointer p_calculator =
        Internals::GetShapeFunctionCalculator<EmbeddedElementData::Dim,
            EmbeddedElementData::NumNodes>(*this, distances);

    // Fluid side
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN, rData.PositiveSideDNDX, rData.PositiveSideWeights,
        GeometryData::GI_GAUSS_2);

    // Fluid side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN, rData.PositiveInterfaceDNDX,
        rData.PositiveInterfaceWeights, GeometryData::GI_GAUSS_2);

    // Fluid side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals, GeometryData::GI_GAUSS_2);

    // Normalize the normals
    const double tolerance = std::pow(1e-3 * this->ElementSize(),Dim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::NormalizeInterfaceNormals(
    typename EmbeddedElementData::InterfaceNormalsType& rNormals,
    double Tolerance) const {
    for (unsigned int i = 0; i < rNormals.size(); ++i) {
        double norm = norm_2(rNormals[i]);
        rNormals[i] /= std::max(norm,Tolerance);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement);
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace Internals {

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

}

}