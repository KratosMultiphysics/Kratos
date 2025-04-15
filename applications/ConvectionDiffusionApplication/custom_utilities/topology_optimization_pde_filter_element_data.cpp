#include "topology_optimization_pde_filter_element_data.h"

namespace Kratos
{

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
TopologyOptimizationPdeFilterElementData<TDim,TNumNodes, TElementIntegratesInTime>::TopologyOptimizationPdeFilterElementData():
    Weight(0.0)
{
}

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
TopologyOptimizationPdeFilterElementData<TDim,TNumNodes, TElementIntegratesInTime>::~TopologyOptimizationPdeFilterElementData()
{
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) {

    // IMPORTANT PARAMETERS AND VARAIBLES INITIALIZATION
    const Geometry< Node >& r_geometry = rElement.GetGeometry();
    const Properties& r_properties = rElement.GetProperties();

    // COMMON PHYSICAL QUANTITIES
    // this->FillFromProperties(Conductivity,CONDUCTIVITY,r_properties);
    // this->FillFromProperties(Decay,DECAY,r_properties);
    // this->FillFromProperties(ConvectionCoefficient,CONVECTION_COEFFICIENT,r_properties);
    this->FillFromNonHistoricalNodalData(PdeFilterDiffusion,PDE_FILTER_DIFFUSION,r_geometry); 
    this->FillFromNonHistoricalNodalData(PdeFilterReaction,PDE_FILTER_REACTION,r_geometry); 

    // Calculate element characteristic size
    ElementSize = ElementSizeCalculator<TDim,TNumNodes>::MinimumElementSize(r_geometry);

    this->FillFromHistoricalNodalData(PdeFilterResult,PDE_FILTER_RESULT,r_geometry);
    this->FillFromHistoricalNodalData(PdeFilterForcing, PDE_FILTER_FORCING, r_geometry);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
int TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::Check(
    const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node >& r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PDE_FILTER_RESULT,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PDE_FILTER_FORCING,r_geometry[i]);
    }

    return 0;
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::UpdateGeometryValues(
    unsigned int IntegrationPointIndex,
    double NewWeight,
    const MatrixRowType& rN,
    const ShapeDerivativesType& rDN_DX)
{
    this->IntegrationPointIndex = IntegrationPointIndex;
    this->Weight = NewWeight;
    noalias(this->N) = rN;
    noalias(this->DN_DX) = rDN_DX;
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalScalarData &rData,
    const Variable<double> &rVariable,
    const Geometry<Node> &rGeometry)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalVectorData& rData,
    const Variable<array_1d<double, 3>>& rVariable,
    const Geometry<Node>& rGeometry)
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
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalTensorData& rData,
    const Variable<Matrix>& rVariable,
    const Geometry<Node>& rGeometry)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        const Matrix& r_nodal_values =
            rGeometry[i].FastGetSolutionStepValue(rVariable);
        rData[i] = r_nodal_values;
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalScalarData& rData, const Variable<double>& rVariable,
    const Geometry<Node>& rGeometry, const unsigned int Step)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable,Step);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalVectorData& rData, const Variable<array_1d<double, 3>>& rVariable,
    const Geometry<Node>& rGeometry, const unsigned int Step)
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
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromNonHistoricalNodalData(
    NodalScalarData& rData,
    const Variable<double>& rVariable,
    const Geometry<Node>& rGeometry)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].GetValue(rVariable);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromNonHistoricalNodalData(
    NodalVectorData& rData,
    const Variable<array_1d<double, 3>>& rVariable,
    const Geometry<Node>& rGeometry)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        const array_1d<double, 3>& r_nodal_values = rGeometry[i].GetValue(rVariable);
        for (size_t j = 0; j < rData.size2(); j++) {
            rData(i, j) = r_nodal_values[j];
        }
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(double& rData,
    const Variable<double>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(Vector& rData,const Variable<Vector>& rVariable,const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(int& rData,
    const Variable<int>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(double& rData,
    const Variable<double>& rVariable, const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(
    Vector& rData,
    const Variable<Vector>& rVariable,
    const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(
    NodalScalarData& rData,
    const Variable<Vector>& rVariable,
    const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TopologyOptimizationPdeFilterElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProperties(double& rData,
    const Variable<double>& rVariable, const Properties& rProperties)
{
    rData = rProperties.GetValue(rVariable);
}

// Triangles
// template class TopologyOptimizationPdeFilterElementData<2,3,false>;
template class TopologyOptimizationPdeFilterElementData<2,3,true>;

// Quadrilaterals
// template class TopologyOptimizationPdeFilterElementData<2,4,false>;
// template class TopologyOptimizationPdeFilterElementData<2,6,false>;
// template class TopologyOptimizationPdeFilterElementData<2,9,false>;
template class TopologyOptimizationPdeFilterElementData<2,4,true>;

// Tetrahedra
// template class TopologyOptimizationPdeFilterElementData<3,4,false>;
template class TopologyOptimizationPdeFilterElementData<3,4,true>;

// Prism
// template class TopologyOptimizationPdeFilterElementData<3,6,false>;
// template class TopologyOptimizationPdeFilterElementData<3,6,true>;

// Hexahedra
// template class TopologyOptimizationPdeFilterElementData<3,8,false>;
// template class TopologyOptimizationPdeFilterElementData<3,8,true>;
// template class TopologyOptimizationPdeFilterElementData<3,27,false>;

} // namespace Kratos
