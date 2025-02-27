#include "transport_topology_optimization_element_data.h"

namespace Kratos
{

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
TransportTopologyOptimizationElementData<TDim,TNumNodes, TElementIntegratesInTime>::TransportTopologyOptimizationElementData():
    Weight(0.0)
{
}

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
TransportTopologyOptimizationElementData<TDim,TNumNodes, TElementIntegratesInTime>::~TransportTopologyOptimizationElementData()
{
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) {

    // IMPORTANT PARAMETERS AND VARAIBLES INITIALIZATION
    const Geometry< Node >& r_geometry = rElement.GetGeometry();
    const Properties& r_properties = rElement.GetProperties();

    // COMMON PHYSICAL QUANTITIES
    // this->FillFromProperties(Conductivity,CONDUCTIVITY,r_properties);
    // this->FillFromProperties(Decay,DECAY,r_properties);
    // this->FillFromProperties(ConvectionCoefficient,CONVECTION_COEFFICIENT,r_properties);
    this->FillFromNonHistoricalNodalData(Conductivity,CONDUCTIVITY,r_geometry); 
    this->FillFromNonHistoricalNodalData(Decay,DECAY,r_geometry); 
    this->FillFromNonHistoricalNodalData(ConvectionCoefficient,CONVECTION_COEFFICIENT,r_geometry); 

    this->FillFromProcessInfo(DeltaTime,DELTA_TIME,rProcessInfo);
    // Calculate element characteristic size
    ElementSize = ElementSizeCalculator<TDim,TNumNodes>::MinimumElementSize(r_geometry);

    // COMMON STABILIZATION QUANTITIES
    DynamicTau = 1.0; //->FillFromProcessInfo(DynamicTau,1,rProcessInfo);
    this->FillFromProcessInfo(TopOptProblemStage,TRANSPORT_TOP_OPT_PROBLEM_STAGE,rProcessInfo);
    // TIME DEPENDENT PROBLEMS NOT YET IMPLEMENTED
    // const Vector& BDFVector = rProcessInfo[BDF_COEFFICIENTS];
    // bdf0 = BDFVector[0];
    // bdf1 = BDFVector[1];
    // bdf2 = BDFVector[2];
    bdf0 = 1.0;
    bdf1 = 1.0;
    bdf2 = 1.0;

    // Functionals Database
    // 0: resistance  : int_{\Omega}{alpha*||u||^2}
    // 1: strain-rate : int_{\Omega}{2*mu*||S||^2} , with S = 1/2*(grad(u)+grad(u)^T) strain-rate tensor
    // 2: vorticity   : int_{\Omega}{2*mu*||R||^2} = int_{\Omega}{mu*||curl(u)||^2} , curl(u) = vorticity vector, R = 1/2*(grad(u)-grad(u)^T) rotation-rate tensor
    // 3: ...
    // 4: ...
    // 5: ...
    // 5 is just a number big enough to contain the acutal database of functionals
    Functional_Weights.resize(20, false);       // Resize to length 5
    this->FillFromProcessInfo(Functional_Weights,FUNCTIONAL_WEIGHTS,rProcessInfo);
    this->FillFromHistoricalNodalData(Optimization_Temperature,OPTIMIZATION_TEMPERATURE,r_geometry);

    // NAVIER-STOKES VARIABLES
    this->FillFromHistoricalNodalData(Temperature,TEMPERATURE,r_geometry);
    this->FillFromHistoricalNodalData(ConvectiveVelocity,VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(MeshVelocity,MESH_VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(SourceTerm, HEAT_FLUX,r_geometry);
    // TIME DEPENDENT PROBLEMS NOT YET IMPLEMENTED
    // this->FillFromHistoricalNodalData(Temperature_OldStep1,TEMPERATURE,r_geometry,1);
    // this->FillFromHistoricalNodalData(Temperature_OldStep2,TEMPERATURE,r_geometry,2);

    // ADJOINT NAVIER-STOKES VARIABLES
    this->FillFromHistoricalNodalData(Temperature_adj,TEMPERATURE_ADJ,r_geometry);
    this->FillFromHistoricalNodalData(SourceTerm_adj, HEAT_FLUX_ADJ,r_geometry);
    // TIME DEPENDENT PROBLEMS NOT YET IMPLEMENTED
    // this->FillFromHistoricalNodalData(Temperature_adj_OldStep1,TEMPERATURE_ADJ,r_geometry,1);
    // this->FillFromHistoricalNodalData(Temperature_adj_OldStep2,TEMPERATURE_ADJ,r_geometry,2);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
int TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::Check(
    const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node >& r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        // CHECK TRANSPORT
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEAT_FLUX,r_geometry[i]);

        // CHECK CONVECTIVE VELOCITY
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,r_geometry[i]);

        // CHECK ADJ TRANSPORT
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE_ADJ,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEAT_FLUX_ADJ,r_geometry[i]);
    }

    return 0;
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::UpdateGeometryValues(
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
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalScalarData &rData,
    const Variable<double> &rVariable,
    const Geometry<Node> &rGeometry)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
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
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
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
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalScalarData& rData, const Variable<double>& rVariable,
    const Geometry<Node>& rGeometry, const unsigned int Step)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable,Step);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
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
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromNonHistoricalNodalData(
    NodalScalarData& rData,
    const Variable<double>& rVariable,
    const Geometry<Node>& rGeometry)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].GetValue(rVariable);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromNonHistoricalNodalData(
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
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(double& rData,
    const Variable<double>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(Vector& rData,const Variable<Vector>& rVariable,const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(int& rData,
    const Variable<int>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(double& rData,
    const Variable<double>& rVariable, const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(
    Vector& rData,
    const Variable<Vector>& rVariable,
    const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(
    NodalScalarData& rData,
    const Variable<Vector>& rVariable,
    const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void TransportTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProperties(double& rData,
    const Variable<double>& rVariable, const Properties& rProperties)
{
    rData = rProperties.GetValue(rVariable);
}

// Triangles
// template class TransportTopologyOptimizationElementData<2,3,false>;
template class TransportTopologyOptimizationElementData<2,3,true>;

// Quadrilaterals
// template class TransportTopologyOptimizationElementData<2,4,false>;
// template class TransportTopologyOptimizationElementData<2,6,false>;
// template class TransportTopologyOptimizationElementData<2,9,false>;
// template class TransportTopologyOptimizationElementData<2,4,true>;

// Tetrahedra
// template class TransportTopologyOptimizationElementData<3,4,false>;
template class TransportTopologyOptimizationElementData<3,4,true>;

// Prism
// template class TransportTopologyOptimizationElementData<3,6,false>;
// template class TransportTopologyOptimizationElementData<3,6,true>;

// Hexahedra
// template class TransportTopologyOptimizationElementData<3,8,false>;
// template class TransportTopologyOptimizationElementData<3,8,true>;
// template class TransportTopologyOptimizationElementData<3,27,false>;

} // namespace Kratos
