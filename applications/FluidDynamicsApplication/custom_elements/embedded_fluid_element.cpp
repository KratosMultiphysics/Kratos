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
    this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);

/*
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            data.UpdateGeometryValues(gauss_weights[g], row(shape_functions, g),
                shape_derivatives[g]);

            this->AddTimeIntegratedSystem(
                data, rLeftHandSideMatrix, rRightHandSideVector);
        }*/
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

    /*this->CalculateGeometryData(
        rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);*/
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
        rData.PositiveInterfaceUnitNormals, GeometryData::GI_GAUSS_2); // This name needs to change? (JC)
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