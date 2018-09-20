//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_FLUID_ELEMENT_DATA_H)
#define KRATOS_FLUID_ELEMENT_DATA_H

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

///@brief Base class for data containers used within FluidElement and derived types.
template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime >
class FluidElementData
{
public:
    ///@name Type Definitions
    ///@{

    using NodalScalarData = array_1d<double,TNumNodes>;

    using NodalVectorData = BoundedMatrix<double,TNumNodes,TDim>;

    using ShapeFunctionsType = array_1d<double,TNumNodes>;

    using ShapeDerivativesType = BoundedMatrix<double,TNumNodes,TDim>;

    #ifdef KRATOS_USE_AMATRIX
    typedef AMatrix::MatrixRow< Matrix > MatrixRowType;
    #else
    typedef boost::numeric::ublas::matrix_row< Matrix > MatrixRowType;
    #endif

    /// Physical space dimension for the problem.
    constexpr static unsigned int Dim = TDim;

    /// Number of nodes of the element.
    constexpr static unsigned int NumNodes = TNumNodes;

    /// Size of the strain and stress vectors (in Voigt notation) for the formulation
    constexpr static unsigned int StrainSize = (TDim-1)*3; // 3 in 2D, 6 in 3D

    /// This lets FluidElement know wether this element requires an external time scheme or not.
    constexpr static bool ElementManagesTimeIntegration = TElementIntegratesInTime;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    FluidElementData();

    /// Destructor
    virtual ~FluidElementData();

    /// (deleted) assignment operator.
    FluidElementData& operator=(FluidElementData const& rOther) = delete;

    /// (deleted) copy constructor.
    FluidElementData(FluidElementData const& rOther) = delete;

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

    /// Shear stress vector in Voigt notation.
    /** It is calculated by the constitutive law in FluidElement::ComputeMaterialResponse.*/
    Vector ShearStress;

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

    void FillFromNodalData(NodalScalarData& rData, const Variable<double>& rVariable, const Geometry<Node<3>>& rGeometry);

    void FillFromNodalData(NodalVectorData& rData, const Variable<array_1d<double,3>>& rVariable, const Geometry<Node<3>>& rGeometry);

    void FillFromHistoricalNodalData(NodalScalarData& rData, const Variable<double>& rVariable, const Geometry<Node<3>>& rGeometry, const unsigned int Step);

    void FillFromHistoricalNodalData(NodalVectorData& rData, const Variable<array_1d<double,3>>& rVariable, const Geometry<Node<3>>& rGeometry, const unsigned int Step);

    void FillFromProcessInfo(double& rData, const Variable<double>& rVariable, const ProcessInfo& rProcessInfo);

    void FillFromProcessInfo(int& rData, const Variable<int>& rVariable, const ProcessInfo& rProcessInfo);

    void FillFromElementData(double& rData, const Variable<double>& rVariable, const Element& rElement);

    void FillFromElementData(NodalScalarData& rData, const Variable<Vector>& rVariable, const Element& rElement);

    void FillFromProperties(double& rData, const Variable<double>& rVariable, const Properties& rProperties);

    ///@}
};

///@}

///@}
}

#endif