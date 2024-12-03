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

#if !defined(KRATOS_FLUID_TOP_OPT_ELEMENT_DATA_H)
#define KRATOS_FLUID_TOP_OPT_ELEMENT_DATA_H

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/constitutive_law.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"
#include "utilities/element_size_calculator.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

///@brief Base class for data containers used within FluidTopologyOptimizationElement and derived types.
template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidTopologyOptimizationElementData
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
    /// Size of the strain and stress vectors (in Voigt notation) for the formulation
    constexpr static unsigned int StrainSize = (TDim-1)*3; // 3 in 2D, 6 in 3D
    /// This lets FluidElement know wether this element requires an external time scheme or not.
    constexpr static bool ElementManagesTimeIntegration = TElementIntegratesInTime;

    ///@}
    ///@name Public Members
    ///@{

    // COMMON PHYSICAL QUANTITIES
    double Density;
    NodalScalarData SoundVelocity;
    double DynamicViscosity;
    NodalScalarData Resistance;     // Darcy's law resistance 
    double DeltaTime;      // Time increment
    double ElementSize;    // Element Characteristic Length (size)

    // COMMON STABILIZATION QUANTITIES
    double DynamicTau;     // Dynamic tau considered in ASGS stabilization coefficients
    double bdf0;
    double bdf1;
    double bdf2;

    int TopOptProblemStage;

    // NAVIER-STOKES VARIABLES
    NodalVectorData Velocity;
    NodalVectorData Velocity_OldStep1;
    NodalVectorData Velocity_OldStep2;
    NodalVectorData MeshVelocity;
    NodalVectorData BodyForce;
    NodalScalarData Pressure;
    NodalScalarData Pressure_OldStep1;
    NodalScalarData Pressure_OldStep2;
    // NS Auxiliary containers for the symbolically-generated matrices
    BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> lhs;
    array_1d<double,TNumNodes*(TDim+1)> rhs;

    // ADJOINT NAVIER-STOKES VARIABLES
    NodalVectorData Velocity_adj;
    NodalVectorData Velocity_adj_OldStep1;
    NodalVectorData Velocity_adj_OldStep2;
    NodalVectorData MeshVelocity_adj;
    NodalVectorData BodyForce_adj;
    NodalScalarData Pressure_adj;
    NodalScalarData Pressure_adj_OldStep1;
    NodalScalarData Pressure_adj_OldStep2;
    // ADJ_NS Auxiliary containers for the symbolically-generated matrices
    BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> lhs_adj;
    array_1d<double,TNumNodes*(TDim+1)> rhs_adj;


    ///@name Life Cycle
    ///@{
    /// Default constructor
    FluidTopologyOptimizationElementData();

    /// Destructor
    virtual ~FluidTopologyOptimizationElementData();

    /// (deleted) assignment operator.
    FluidTopologyOptimizationElementData& operator=(FluidTopologyOptimizationElementData const& rOther) = delete;

    /// (deleted) copy constructor.
    FluidTopologyOptimizationElementData(FluidTopologyOptimizationElementData const& rOther) = delete;

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

    /// Strain rate (symmetric gradient of velocity) vector in Voigt notation.
    /** It is calculated by the constitutive law in FluidElement::ComputeMaterialResponse.*/
    Vector StrainRate;
    // For Adjoint NS problem
    Vector StrainRate_adj;

    /// Shear stress vector in Voigt notation.
    /** It is calculated by the constitutive law in FluidElement::ComputeMaterialResponse.*/
    Vector ShearStress;
    // For Adjoint NS problem
    Vector ShearStress_adj;

    /// Constitutive tensor C (expressed as a Matrix).
    /** It is calculated by the constitutive law in FluidElement::ComputeMaterialResponse.*/
    Matrix C;

    /// Constitutive law configuration (stored here to avoid re-initialization within the element).
    ConstitutiveLaw::Parameters ConstitutiveLawValues;

    /// Effective viscosity (in dynamic units) produced by the constitutive law
    double EffectiveViscosity;

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
        KRATOS_WARNING("FluidTopologyOptimizationElementData") << "\'FillFromNodalData\' is deprecated. Use \'FillFromHistoricalNodalData\' instead." << std::endl;
        FillFromHistoricalNodalData(rData, rVariable, rGeometry);
    }

    //TODO: This needs to be removed
    void FillFromNodalData(
        NodalVectorData &rData,
        const Variable<array_1d<double, 3>> &rVariable,
        const Geometry<Node> &rGeometry)
    {
        KRATOS_WARNING("FluidTopologyOptimizationElementData") << "\'FillFromNodalData\' is deprecated. Use \'FillFromHistoricalNodalData\' instead." << std::endl;
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

    void FillFromElementData(double& rData, const Variable<double>& rVariable, const Element& rElement);

    void FillFromElementData(Vector& rData, const Variable<Vector>& rVariable, const Element& rElement);

    void FillFromElementData(NodalScalarData& rData, const Variable<Vector>& rVariable, const Element& rElement);

    void FillFromProperties(double& rData, const Variable<double>& rVariable, const Properties& rProperties);

    ///@}
};

///@}

///@}
}

#endif // KRATOS_FLUID_ELEMENT_DATA_H
