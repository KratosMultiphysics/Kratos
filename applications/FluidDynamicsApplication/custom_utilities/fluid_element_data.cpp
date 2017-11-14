#include "fluid_element_data.h"

namespace Kratos
{

template< size_t TDim, size_t TNumNodes >
FluidElementData<TDim,TNumNodes>::FluidElementData()
{
}

template< size_t TDim, size_t TNumNodes >
FluidElementData<TDim,TNumNodes>::~FluidElementData()
{
}

template <size_t TDim, size_t TNumNodes>
void FluidElementData<TDim, TNumNodes>::UpdateGeometryValues(double NewWeight,
    boost::numeric::ublas::matrix_row<Kratos::Matrix> rN,
    const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rDN_DX)
{
    this->Weight = NewWeight;
    noalias(this->N) = rN;
    noalias(this->DN_DX) = rDN_DX;
}

template <size_t TDim, size_t TNumNodes>
void FluidElementData<TDim, TNumNodes>::FillFromNodalData(
    NodalScalarData& rData, const Variable<double>& rVariable,
    const Geometry<Node<3>>& rGeometry)
{
    noalias(rData) = array_1d<double,TNumNodes>(TNumNodes,0.0);
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable);
    }
}

template <size_t TDim, size_t TNumNodes>
void FluidElementData<TDim, TNumNodes>::FillFromNodalData(NodalVectorData& rData,
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

template <size_t TDim, size_t TNumNodes>
void FluidElementData<TDim, TNumNodes>::FillFromHistoricalNodalData(
    NodalScalarData& rData, const Variable<double>& rVariable,
    const Geometry<Node<3>>& rGeometry, const unsigned int Step)
{
    noalias(rData) = array_1d<double,TNumNodes>(TNumNodes,0.0);
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable,Step);
    }
}

template <size_t TDim, size_t TNumNodes>
void FluidElementData<TDim, TNumNodes>::FillFromHistoricalNodalData(
    NodalVectorData& rData, const Variable<array_1d<double, 3>>& rVariable,
    const Geometry<Node<3>>& rGeometry, const unsigned int Step)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        const array_1d<double, 3>& r_nodal_values =
            rGeometry[i].FastGetSolutionStepValue(rVariable);
        for (size_t j = 0; j < rData.size2(); j++) {
            rData(i, j) = r_nodal_values[j];
        }
    }
}

template <size_t TDim, size_t TNumNodes>
void FluidElementData<TDim, TNumNodes>::FillFromProcessInfo(double& rData,
    const Variable<double>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes>
void FluidElementData<TDim, TNumNodes>::FillFromElementData(double& rData,
    const Variable<double>& rVariable, const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template class FluidElementData<2,3>;
template class FluidElementData<3,4>;

}