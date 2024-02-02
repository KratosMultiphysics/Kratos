//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/two_fluid_navier_stokes_alpha_method.h"
#include "custom_utilities/fluid_element_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template< class TElementData >
class TwoFluidNavierStokesAlphaMethodDiscontinuous : public TwoFluidNavierStokesAlphaMethod<TElementData>
{
public:

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TwoFluidNavierStokesAlphaMethodDiscontinuous);

    ///@name Type Definitions
    ///@{

    using BaseType = TwoFluidNavierStokesAlphaMethod<TElementData>;

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef std::vector< Dof<double>::Pointer > DofsVectorType;
    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;
    typedef typename BaseType::ShapeFunctionsType ShapeFunctionsType;
    typedef typename BaseType::ShapeFunctionDerivativesType ShapeFunctionDerivativesType;
    typedef typename BaseType::ShapeFunctionDerivativesArrayType ShapeFunctionDerivativesArrayType;

    constexpr static SizeType Dim = BaseType::Dim;
    constexpr static SizeType NumNodes = BaseType::NumNodes;
    constexpr static SizeType BlockSize = BaseType::BlockSize;
    constexpr static SizeType LocalSize = BaseType::LocalSize;
    constexpr static SizeType StrainSize = (Dim-1)*3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constuctor.
    TwoFluidNavierStokesAlphaMethodDiscontinuous(IndexType NewId = 0)
        : TwoFluidNavierStokesAlphaMethod<TElementData>(NewId)
    {};

    /// Constructor using an array of nodes.
    TwoFluidNavierStokesAlphaMethodDiscontinuous(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
        : TwoFluidNavierStokesAlphaMethod<TElementData>(NewId, ThisNodes)
    {};

    /// Constructor using a geometry object.
    TwoFluidNavierStokesAlphaMethodDiscontinuous(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : TwoFluidNavierStokesAlphaMethod<TElementData>(NewId, pGeometry)
    {};

    /// Constuctor using geometry and properties.
    TwoFluidNavierStokesAlphaMethodDiscontinuous(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties)
        : TwoFluidNavierStokesAlphaMethod<TElementData>(NewId, pGeometry, pProperties)
    {};

    /// Copy constructor.
    TwoFluidNavierStokesAlphaMethodDiscontinuous(TwoFluidNavierStokesAlphaMethodDiscontinuous const& rOther) = delete;

    /// Destructor.
    virtual ~TwoFluidNavierStokesAlphaMethodDiscontinuous() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ComputeGaussPointLHSContribution(
        TElementData& rData,
        MatrixType& rLHS) override;

    void ComputeGaussPointRHSContribution(
        TElementData& rData,
        VectorType& rRHS) override;

	void ComputeGaussPointEnrichmentContributions(
		TElementData& rData,
		MatrixType& rV,
		MatrixType& rH,
		MatrixType& rKee,
		VectorType& rRHS_ee) override;

    void CalculateStrainRate(TElementData& rData) const override;

    double CalculateArtificialDynamicViscositySpecialization(TElementData &rData) const;

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
private:
    ///@name Member Variables
    ///@{

    array_1d<double,Dim> mVelEnrPos;
    array_1d<double,Dim> mVelEnrNeg;
    array_1d<double,NumNodes> mPresEnr;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Split shape functions computation auxiliar method
     * This method computes the standard and enrichment shape functions for a split element
     * @param rData Element data container
     * @param rStandardShapeFunctionsPos Positive side standard shape functions values
     * @param rStandardShapeFunctionsNeg Negative side standard shape functions values
     * @param rEnrichedShapeFunctionsPos Positive side enrichment shape functions values
     * @param rEnrichedShapeFunctionsNeg Negative side enrichment shape functions values
     * @param rBubbleShapeFunctionsPos Positive side enrichment bubble shape functions values
     * @param rBubbleShapeFunctionsNeg Negative side enrichment bubble shape functions values
     * @param rStandardShapeDerivativesPos  Positive side standard shape functions derivatives values
     * @param rStandardShapeDerivativesNeg  Negative side standard shape functions derivatives values
     * @param rEnrichedShapeDerivativesPos Positive side enrichment shape functions derivatives values
     * @param rEnrichedShapeDerivativesNeg Negative side enrichment shape functions derivatives values
     * @param rBubbleShapeDerivativesPos  Positive side enrichment bubble shape functions derivatives values
     * @param rBubbleShapeDerivativesNeg  Negative side enrichment bubble shape functions derivatives values
     */
    void ComputeSplitting(
        TElementData &rData,
        MatrixType &rStandardShapeFunctionsPos,
        MatrixType &rStandardShapeFunctionsNeg,
        MatrixType &rEnrichedShapeFunctionsPos,
        MatrixType &rEnrichedShapeFunctionsNeg,
        VectorType &rBubbleShapeFunctionsPos,
        VectorType &rBubbleShapeFunctionsNeg,
        GeometryType::ShapeFunctionsGradientsType &rStandardShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType &rStandardShapeDerivativesNeg,
        GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType &rEnrichedShapeDerivativesNeg,
        DenseVector<array_1d<double,Dim>> & rBubbleShapeDerivativesPos,
        DenseVector<array_1d<double,Dim>> & rBubbleShapeDerivativesNeg);

    /**
     * @brief Set up the element's data for a cut element and constitutive law for the current integration point
     * This function sets the element kinematics (standard, Ausas and enrichment FE spaces) for the current integration point
     * @param rData Element data container
     * @param IntegrationPointIndex Index of current integration point
     * @param IntegrationPointWeight Weight of current integration point
     * @param rStandardShapeFunctions Values of standard shape functions at integration point
     * @param rStandardShapeFunctionsGradients Values of standard shape functions gradients at integration point
     * @param rEnrichedShapeFunctions Values of enrichment shape functions at integration point
     * @param rEnrichedShapeFunctionsGradients Values of enrichment shape functions gradients at integration point
     * @param rBubbleShapeFunction Value of the enrichment bubble function at integration point
     * @param rBubbleShapeFunctionGradients Value of the enrichment bubble function gradients at integration point
     */
    void UpdateIntegrationPointDataDiscontinuous(
        TElementData& rData,
        IndexType IntegrationPointIndex,
        double IntegrationPointWeight,
        const typename TElementData::MatrixRowType& rStandardShapeFunctions,
        const typename TElementData::ShapeDerivativesType& rStandardShapeFunctionsGradients,
        const typename TElementData::MatrixRowType& rEnrichedShapeFunctions,
        const typename TElementData::ShapeDerivativesType& rEnrichedShapeFunctionsGradients,
        const double rBubbleShapeFunction,
        const array_1d<double, Dim>& rBubbleShapeFunctionGradients) const;

	void CondenseAndSaveEnrichmentWithContinuity(
		const TElementData& rData,
		Matrix& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		const MatrixType& rVTot,
		const MatrixType& rHTot,
		MatrixType& rKeeTot,
		const VectorType& rRHSeeTot);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class TwoFluidNavierStokesAlphaMethodDiscontinuous
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< class TElementData >
inline std::istream& operator >> (
    std::istream& rIStream,
    TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const TwoFluidNavierStokesAlphaMethodDiscontinuous<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.
