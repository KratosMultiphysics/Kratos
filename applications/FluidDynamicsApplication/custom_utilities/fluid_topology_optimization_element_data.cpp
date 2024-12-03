#include "fluid_topology_optimization_element_data.h"

namespace Kratos
{

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
FluidTopologyOptimizationElementData<TDim,TNumNodes, TElementIntegratesInTime>::FluidTopologyOptimizationElementData():
    Weight(0.0)
{
}

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
FluidTopologyOptimizationElementData<TDim,TNumNodes, TElementIntegratesInTime>::~FluidTopologyOptimizationElementData()
{
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) {
    ConstitutiveLawValues = ConstitutiveLaw::Parameters(rElement.GetGeometry(),rElement.GetProperties(),rProcessInfo);
    if (StrainRate.size() != StrainSize) {
        StrainRate.resize(StrainSize);
    }
    if (ShearStress.size() != StrainSize) {
        ShearStress.resize(StrainSize);
    }
    //  Resize for ADJOINT NS
    if (StrainRate_adj.size() != StrainSize) {
        StrainRate_adj.resize(StrainSize);
    }
    if (ShearStress_adj.size() != StrainSize) {
        ShearStress_adj.resize(StrainSize);
    }
    if (C.size1() != StrainSize || C.size2() != StrainSize) {
        C.resize(StrainSize,StrainSize,false);
    }

    // Set constitutive law flags:
    Flags& cl_options = ConstitutiveLawValues.GetOptions();
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    ConstitutiveLawValues.SetStrainVector(StrainRate);   //this is the input parameter
    ConstitutiveLawValues.SetStressVector(ShearStress);  //this is an ouput parameter

    // // Set for ADJOINT
    // ConstitutiveLawValues.SetStrainVector(StrainRate_adj);   //this is the input parameter
    // ConstitutiveLawValues.SetStressVector(ShearStress_adj);  //this is an ouput parameter

    ConstitutiveLawValues.SetConstitutiveMatrix(C);      //this is an ouput parameter

    // IMPORTANT PARAMETERS AND VARAIBLES INITIALIZATION
    const Geometry< Node >& r_geometry = rElement.GetGeometry();
    const Properties& r_properties = rElement.GetProperties();

    // COMMON PHYSICAL QUANTITIES
    this->FillFromProperties(Density,DENSITY,r_properties);
    // this->FillFromNonHistoricalNodalData(SoundVelocity, SOUND_VELOCITY, r_geometry);
    this->FillFromProperties(DynamicViscosity,DYNAMIC_VISCOSITY,r_properties); //TODO: This needs to be retrieved from the EffectiveViscosity of the constitutive law
    // Resistance has to be used as a nodal value, since in this way evaluating the functional derivatives
    // w.r.t the design parameters is easier.
    this->FillFromNonHistoricalNodalData(Resistance, RESISTANCE, r_geometry); 
    this->FillFromProcessInfo(DeltaTime,DELTA_TIME,rProcessInfo);
    // Calculate element characteristic size
    ElementSize = ElementSizeCalculator<TDim,TNumNodes>::MinimumElementSize(r_geometry);

    // COMMON STABILIZATION QUANTITIES
    this->FillFromProcessInfo(DynamicTau,DYNAMIC_TAU,rProcessInfo);
    this->FillFromProcessInfo(TopOptProblemStage,FLUID_TOP_OPT_PROBLEM_STAGE,rProcessInfo);
    // TIME DEPENDENT PROBLEMS NOT YET IMPLEMENTED
    // const Vector& BDFVector = rProcessInfo[BDF_COEFFICIENTS];
    // bdf0 = BDFVector[0];
    // bdf1 = BDFVector[1];
    // bdf2 = BDFVector[2];
    bdf0 = 1.0;
    bdf1 = 1.0;
    bdf2 = 1.0;

    // NAVIER-STOKES VARIABLES
    this->FillFromHistoricalNodalData(Velocity,VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(MeshVelocity,MESH_VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(BodyForce,BODY_FORCE,r_geometry);
    this->FillFromHistoricalNodalData(Pressure,PRESSURE,r_geometry);
    // TIME DEPENDENT PROBLEMS NOT YET IMPLEMENTED
    // this->FillFromHistoricalNodalData(Velocity_OldStep1,VELOCITY,r_geometry,1);
    // this->FillFromHistoricalNodalData(Velocity_OldStep2,VELOCITY,r_geometry,2);
    // this->FillFromHistoricalNodalData(Pressure_OldStep1,PRESSURE,r_geometry,1);
    // this->FillFromHistoricalNodalData(Pressure_OldStep2,PRESSURE,r_geometry,2);

    // ADJOINT NAVIER-STOKES VARIABLES
    this->FillFromHistoricalNodalData(Velocity_adj,VELOCITY_ADJ,r_geometry);
    this->FillFromHistoricalNodalData(MeshVelocity_adj,MESH_VELOCITY_ADJ,r_geometry);
    this->FillFromHistoricalNodalData(BodyForce_adj,BODY_FORCE_ADJ,r_geometry);
    this->FillFromHistoricalNodalData(Pressure_adj,PRESSURE_ADJ,r_geometry);
    // TIME DEPENDENT PROBLEMS NOT YET IMPLEMENTED
    // this->FillFromHistoricalNodalData(Velocity_adj_OldStep1,VELOCITY_ADJ,r_geometry,1);
    // this->FillFromHistoricalNodalData(Velocity_adj_OldStep2,VELOCITY_ADJ,r_geometry,2);
    // this->FillFromHistoricalNodalData(Pressure_adj_OldStep1,PRESSURE_ADJ,r_geometry,1);
    // this->FillFromHistoricalNodalData(Pressure_adj_OldStep2,PRESSURE_ADJ,r_geometry,2);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
int FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::Check(
    const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node >& r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        // CHECK NS
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_geometry[i]);

        // CHECK ADJ NS
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_ADJ,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY_ADJ,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE_ADJ,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_ADJ,r_geometry[i]);
    }

    return 0;
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::UpdateGeometryValues(
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
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalScalarData &rData,
    const Variable<double> &rVariable,
    const Geometry<Node> &rGeometry)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
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
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
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
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
    NodalScalarData& rData, const Variable<double>& rVariable,
    const Geometry<Node>& rGeometry, const unsigned int Step)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].FastGetSolutionStepValue(rVariable,Step);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromHistoricalNodalData(
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
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromNonHistoricalNodalData(
    NodalScalarData& rData,
    const Variable<double>& rVariable,
    const Geometry<Node>& rGeometry)
{
    for (size_t i = 0; i < TNumNodes; i++) {
        rData[i] = rGeometry[i].GetValue(rVariable);
    }
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromNonHistoricalNodalData(
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
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(double& rData,
    const Variable<double>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProcessInfo(int& rData,
    const Variable<int>& rVariable, const ProcessInfo& rProcessInfo)
{
    rData = rProcessInfo.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(double& rData,
    const Variable<double>& rVariable, const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(
    Vector& rData,
    const Variable<Vector>& rVariable,
    const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromElementData(
    NodalScalarData& rData,
    const Variable<Vector>& rVariable,
    const Element& rElement)
{
    rData = rElement.GetValue(rVariable);
}

template <size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime>
void FluidTopologyOptimizationElementData<TDim, TNumNodes, TElementIntegratesInTime>::FillFromProperties(double& rData,
    const Variable<double>& rVariable, const Properties& rProperties)
{
    rData = rProperties.GetValue(rVariable);
}

// Triangles
template class FluidTopologyOptimizationElementData<2,3,false>;
template class FluidTopologyOptimizationElementData<2,3,true>;

// Quadrilaterals
template class FluidTopologyOptimizationElementData<2,4,false>;
template class FluidTopologyOptimizationElementData<2,6,false>;
template class FluidTopologyOptimizationElementData<2,9,false>;
template class FluidTopologyOptimizationElementData<2,4,true>;

// Tetrahedra
template class FluidTopologyOptimizationElementData<3,4,false>;
template class FluidTopologyOptimizationElementData<3,4,true>;

// Prism
template class FluidTopologyOptimizationElementData<3,6,false>;
template class FluidTopologyOptimizationElementData<3,6,true>;

// Hexahedra
template class FluidTopologyOptimizationElementData<3,8,false>;
template class FluidTopologyOptimizationElementData<3,8,true>;
template class FluidTopologyOptimizationElementData<3,27,false>;

} // namespace Kratos
