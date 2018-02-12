#include "fluid_element_data.h"

namespace Kratos
{

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
FluidElementData<TDim,TNumNodes, TElementIntegratesInTime>::FluidElementData():
    Weight(0.0)
{
}

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
FluidElementData<TDim,TNumNodes, TElementIntegratesInTime>::~FluidElementData()
{
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) {

    ConstitutiveLawValues = ConstitutiveLaw::Parameters(rElement.GetGeometry(),rElement.GetProperties(),rProcessInfo);

    if (StrainRate.size() != StrainSize) {
        StrainRate.resize(StrainSize);
    }
    if (Stress.size() != StrainSize) {
        Stress.resize(StrainSize);
    }
    if (C.size1() != StrainSize || C.size2() != StrainSize) {
        C.resize(StrainSize,StrainSize);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
int FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::Check(
    const Element& rElement, const ProcessInfo& rProcessInfo)
{
    return 0;
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::UpdateGeometryValues(double NewWeight,
    const boost::numeric::ublas::matrix_row<Kratos::Matrix> rN,
    const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rDN_DX)
{
    this->Weight = NewWeight;
    noalias(this->N) = rN;
    noalias(this->DN_DX) = rDN_DX;
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromNodalData(
    NodalScalarData& rData, const Variable<double>& rVariable,
    const Geometry<Node<3>>& rGeometry)
{
    noalias(rData) = array_1d<double,TNumNodes>(TNumNodes,0.0);
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromNodalData(NodalVectorData& rData,
    const Variable<array_1d<double, 3>>& rVariable,
    const Geometry<Node<3>>& rGeometry) 
{
    for (size_t i = 0; i < TNumNodes; i++) {
        const array_1d<double, 3>& r_nodal_values =
            rGeometry[i].FastGetSolutionStepValue(rVariable);
        for (size_t j = 0; j < rData.size2(); j++) {
            rData(i, j) = r_nodal_values[j];
        }
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalScalarData& rData, const Variable<double>& rVariable,
    const Geometry<Node<3>>& rGeometry, const unsigned int Step)
{
    noalias(rData) = array_1d<double,TNumNodes>(TNumNodes,0.0);
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable,Step);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalVectorData& rData, const Variable<array_1d<double, 3>>& rVariable,
    const Geometry<Node<3>>& rGeometry, const unsigned int Step)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        const array_1d<double, 3>& r_nodal_values =
            rGeometry[i].FastGetSolutionStepValue(rVariable,Step);
        for (size_t j = 0; j < rData.size2(); j++) {
            rData(i, j) = r_nodal_values[j];
        }
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(double& rData,
    const Variable<double>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(int& rData,
    const Variable<int>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(double& rData,
    const Variable<double>& rVariable, const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProperties(double& rData,
    const Variable<double>& rVariable, const Properties& rProperties)
{
    rData = rProperties.GetValue(rVariable);
}

// Triangles
template class FluidElementData<2,3,false>;
template class FluidElementData<2,3,true>;

// Quadrilaterals
template class FluidElementData<2,4,false>;

// Tetrahedra
template class FluidElementData<3,4,false>;
template class FluidElementData<3,4,true>;

// Hexahedra
template class FluidElementData<3,8,false>;

}