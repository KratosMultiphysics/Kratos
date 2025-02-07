//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Gianmarco Boscolo
//

#if !defined(KRATOS_TRANSPORT_TOP_OPT_ELEMENT_DATA_H)
#define KRATOS_TRANSPORT_TOP_OPT_ELEMENT_DATA_H

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/checks.h"
#include "includes/process_info.h"

// Application includes
#include "convection_diffusion_application_variables.h"
#include "utilities/element_size_calculator.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

///@brief Base class for data containers used within FluidTransportTopologyOptimizationElement and derived types.
template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) TransportTopologyOptimizationElementData
{
public:
    ///@name Type Definitions
    ///@{
    using NodalScalarData = array_1d<double,TNumNodes>;
    using NodalVectorData = BoundedMatrix<double,TNumNodes,TDim>;
    using NodalTensorData = std::array<BoundedMatrix<double,TDim,TDim>,TNumNodes>;
    using ShapeFunctionsType = array_1d<double,TNumNodes>;
    using ShapeDerivativesType = BoundedMatrix<double,TNumNodes,TDim>;
    using MatrixRowType = MatrixRow< Matrix >;
    /// Physical space dimension for the problem.
    constexpr static unsigned int Dim = TDim;
    /// Number of nodes of the element.
    constexpr static unsigned int NumNodes = TNumNodes;
    /// This lets FluidElement know wether this element requires an external time scheme or not.
    constexpr static bool ElementManagesTimeIntegration = TElementIntegratesInTime;

    ///@}
    ///@name Public Members
    ///@{

    // COMMON PHYSICAL QUANTITIES
    double Conductivity;
    double Decay;
    double DeltaTime;      // Time increment
    double ElementSize;    // Element Characteristic Length (size)

    // COMMON STABILIZATION QUANTITIES
    double DynamicTau;     // Dynamic tau considered in ASGS stabilization coefficients
    double bdf0;
    double bdf1;
    double bdf2;

    Vector Functional_Weights; // weigths of the functional terms

    int TopOptProblemStage;

    // TRANSPORT VARIABLES
    NodalScalarData Temperature;
    NodalScalarData Temperature_OldStep1;;
    NodalScalarData Temperature_OldStep2;;
    NodalVectorData ConvectiveVelocity;
    NodalVectorData MeshVelocity;
    NodalScalarData SourceTerm;
    // NS Auxiliary containers for the symbolically-generated matrices
    BoundedMatrix<double,TNumNodes,TNumNodes> lhs;
    array_1d<double,TNumNodes> rhs;

    // ADJOINT TRANSPORT VARIABLES
    NodalScalarData Temperature_adj;
    NodalScalarData Temperature_adj_OldStep1;;
    NodalScalarData Temperature_adj_OldStep2;;
    NodalScalarData SourceTerm_adj;
    // ADJ_NS Auxiliary containers for the symbolically-generated matrices
    BoundedMatrix<double,TNumNodes,TNumNodes> lhs_adj;
    array_1d<double,TNumNodes> rhs_adj;


    ///@name Life Cycle
    ///@{
    /// Default constructor
    TransportTopologyOptimizationElementData();

    /// Destructor
    virtual ~TransportTopologyOptimizationElementData();

    /// (deleted) assignment operator.
    TransportTopologyOptimizationElementData& operator=(TransportTopologyOptimizationElementData const& rOther) = delete;

    /// (deleted) copy constructor.
    TransportTopologyOptimizationElementData(TransportTopologyOptimizationElementData const& rOther) = delete;

    ///@}
    ///@name Public Operations
    ///@{

    virtual void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo);

    static int Check(const Element& rElement, const ProcessInfo& rProcessInfo);

    virtual void UpdateGeometryValues(
        unsigned int IntegrationPointIndex,
        double NewWeight,
        const MatrixRowType& rN,
        const ShapeDerivativesType& rDN_DX);

    ///@}
    ///@name Public Members
    ///@{

    unsigned int IntegrationPointIndex;

    double Weight;

    ShapeFunctionsType N;

    ShapeDerivativesType DN_DX;

    ///@}
protected:

    ///@name Protected Operations
    ///@{

    //TODO: This needs to be removed
    void FillFromNodalData(
        NodalScalarData &rData,
        const Variable<double> &rVariable,
        const Geometry<Node> &rGeometry)
    {
        KRATOS_WARNING("TransportTopologyOptimizationElementData") << "\'FillFromNodalData\' is deprecated. Use \'FillFromHistoricalNodalData\' instead." << std::endl;
        FillFromHistoricalNodalData(rData, rVariable, rGeometry);
    }

    //TODO: This needs to be removed
    void FillFromNodalData(
        NodalVectorData &rData,
        const Variable<array_1d<double, 3>> &rVariable,
        const Geometry<Node> &rGeometry)
    {
        KRATOS_WARNING("TransportTopologyOptimizationElementData") << "\'FillFromNodalData\' is deprecated. Use \'FillFromHistoricalNodalData\' instead." << std::endl;
        FillFromHistoricalNodalData(rData, rVariable, rGeometry);
    }

    void FillFromHistoricalNodalData(
        NodalScalarData &rData,
        const Variable<double> &rVariable,
        const Geometry<Node> &rGeometry);

    void FillFromHistoricalNodalData(
        NodalVectorData &rData,
        const Variable<array_1d<double, 3>> &rVariable,
        const Geometry<Node> &rGeometry);

    void FillFromHistoricalNodalData(
        NodalTensorData& rData,
        const Variable<Matrix>& rVariable,
        const Geometry<Node>& rGeometry);

    void FillFromHistoricalNodalData(NodalScalarData& rData, const Variable<double>& rVariable, const Geometry<Node>& rGeometry, const unsigned int Step);

    void FillFromHistoricalNodalData(NodalVectorData& rData, const Variable<array_1d<double,3>>& rVariable, const Geometry<Node>& rGeometry, const unsigned int Step);

    void FillFromNonHistoricalNodalData(
        NodalScalarData& rData,
        const Variable<double>& rVariable,
        const Geometry<Node>& rGeometry);

    void FillFromNonHistoricalNodalData(
        NodalVectorData& rData,
        const Variable<array_1d<double,3>>& rVariable,
        const Geometry<Node>& rGeometry);

    void FillFromProcessInfo(double& rData, const Variable<double>& rVariable, const ProcessInfo& rProcessInfo);

    void FillFromProcessInfo(int& rData, const Variable<int>& rVariable, const ProcessInfo& rProcessInfo);

    void FillFromProcessInfo(Vector& rData, const Variable<Vector>& rVariable, const ProcessInfo& rProcessInfo);

    void FillFromElementData(double& rData, const Variable<double>& rVariable, const Element& rElement);

    void FillFromElementData(Vector& rData, const Variable<Vector>& rVariable, const Element& rElement);

    void FillFromElementData(NodalScalarData& rData, const Variable<Vector>& rVariable, const Element& rElement);

    void FillFromProperties(double& rData, const Variable<double>& rVariable, const Properties& rProperties);

    ///@}
};

///@}

///@}
}

#endif // KRATOS_TRANSPORT_TOPOLOGY_OPTIMIZATION_ELEMENT_DATA_H
