//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Uxue Chasco
//
//

#if !defined(KRATOS_TWO_FLUID_NAVIER_STOKES_ALPHA_METHOD)
#define  KRATOS_TWO_FLUID_NAVIER_STOKES_ALPHA_METHOD

// System includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
#include "custom_elements/two_fluid_navier_stokes.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "utilities/geometry_utilities.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{

/*The "TwoFluidNavierStokesAlphaMethod" element is an element based on the Variational Multiscale Stabilization technique (VMS)
* which is designed for the solution of a two fluids (air and a non-air one) problems.
*
* A distinctive feature of the element is the use of 4 LOCAL enrichment functions, which allows to model
* a discontinuity in both the pressure field and in its gradient.
* The enrichment functions are obtained by duplicating all of the degrees of freedom of the element.
* Since the enrichment is performed elementwise, a purely local static condensation step is performed,
* meaning that no extra degrees of freedom are added..
*
* Since a jump in the pressure can be considered, the element shall be able to habdle moderate changes of the viscosity
* between the two fluids to be considered*/

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
class TwoFluidNavierStokesAlphaMethod : public TwoFluidNavierStokes<TElementData>
{
public:

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TwoFluidNavierStokesAlphaMethod);

    ///@name Type Definitions
    ///@{

    using BaseType = TwoFluidNavierStokes<TElementData>;

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
    typedef typename TwoFluidNavierStokes<TElementData>::ShapeFunctionsType ShapeFunctionsType;
    typedef typename TwoFluidNavierStokes<TElementData>::ShapeFunctionDerivativesType ShapeFunctionDerivativesType;
    typedef typename TwoFluidNavierStokes<TElementData>::ShapeFunctionDerivativesArrayType ShapeFunctionDerivativesArrayType;
    constexpr static unsigned int Dim = TwoFluidNavierStokes<TElementData>::Dim;
    constexpr static unsigned int NumNodes = TwoFluidNavierStokes<TElementData>::NumNodes;
    constexpr static unsigned int BlockSize = TwoFluidNavierStokes<TElementData>::BlockSize;
    constexpr static unsigned int LocalSize = TwoFluidNavierStokes<TElementData>::LocalSize;
    constexpr static unsigned int StrainSize = (Dim - 1) * 3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constuctor.
    /**
    * @param NewId Index number of the new element (optional)
    */
    TwoFluidNavierStokesAlphaMethod(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
    * @param NewId Index of the new element
    * @param ThisNodes An array containing the nodes of the new element
    */
    TwoFluidNavierStokesAlphaMethod(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    */
    TwoFluidNavierStokesAlphaMethod(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    * @param pProperties Pointer to the element's properties
    */
    TwoFluidNavierStokesAlphaMethod(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~TwoFluidNavierStokesAlphaMethod();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
    * Returns a pointer to a new TwoFluidNavierStokesAlphaMethod element, created using given input.
    * @param NewId the ID of the new element
    * @param ThisNodes the nodes of the new element
    * @param pProperties the properties assigned to the new element
    * @return a Pointer to the new element
    */
    Element::Pointer Create(IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    /// Create a new element of this type using given geometry
    /**
    * Returns a pointer to a new FluidElement element, created using given input.
    * @param NewId the ID of the new element
    * @param pGeom a pointer to the geomerty to be used to create the element
    * @param pProperties the properties assigned to the new element
    * @return a Pointer to the new element
    */
    Element::Pointer Create(IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void CalculateOnIntegrationPoints(
        const Variable<double> &rVariable,
        std::vector<double> &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override;

    void Calculate(
        const Variable<double> &rVariable,
        double &rOutput,
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

    /**
     * @brief Computes the enriched LHS/RHS terms associated with the pressure stabilizations at the interface
     * @param rInterfaceWeightsNeg Negative side weights for the interface-gauss-points
     * @param rEnrInterfaceShapeFunctionPos Enriched shape functions at the interface-gauss-points Positive side
     * @param rEnrInterfaceShapeFunctionNeg Enriched shape functions at the interface-gauss-points Negative side
     * @param rInterfaceShapeDerivativesNeg Shape functions derivatives at the interface-gauss-points
     * @param rKeeTot Pressure enrichment contribution related to pressure enrichment DOFs
     * @param rRHSeeTot Right Hand Side vector associated to the pressure enrichment DOFs
     */
    void PressureGradientStabilization(
        const TElementData &rData,
        const Vector &rInterfaceWeights,
        const Matrix &rEnrInterfaceShapeFunctionPos,
        const Matrix &rEnrInterfaceShapeFunctionNeg,
        const GeometryType::ShapeFunctionsGradientsType &rInterfaceShapeDerivatives,
        MatrixType &rKeeTot,
        VectorType &rRHSeeTot) override;

    /**
     * @brief Computes the LHS Gauss pt. contribution
     * This method computes the contribution to the LHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     */
    void ComputeGaussPointLHSContribution(
        TElementData& rData,
        MatrixType& rLHS) override;

    /**
     * @brief Computes the RHS Gaus  pt. contribution
     * This method computes the contribution to the RHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    void ComputeGaussPointRHSContribution(
        TElementData& rData,
        VectorType& rRHS) override;

    /**
     * @brief Computes the pressure enrichment contributions
     * This method computes the pressure enrichment contributions for
     * a Gauss pt. in both the left hand side and righ hand side of the equations.
     * @param rData Reference to the element data container
     * @param rV Contribution related to the pressure enrichment DOFs in the N-S standard equations
     * @param rH Contribution related to the standard velocity and pressure DOFs in the enrichment equations
     * @param rKee Contribution related to the pressure enrichment DOFs in the enrichment equations
     * @param rRHS_ee Right Hand Side of the enrichment equations
     */
	void ComputeGaussPointEnrichmentContributions(
		TElementData& rData,
		MatrixType& rV,
		MatrixType& rH,
		MatrixType& rKee,
		VectorType& rRHS_ee) override;

    /**
     * @brief Calculate the strain rate
     * In this function we calculate the strain rate at the mid step
     * @param rData Data container with the input velocity and gradients and output strain rate vector
     */
    void CalculateStrainRate(TElementData& rData) const override;

    /**
     * @brief Calculate the artificial dynamic viscosity.
     * In this function we calculate the artificial dynamic viscosity in each gauss point.
     * @param rData Data container
     */

    double CalculateArtificialDynamicViscositySpecialization(TElementData &rData) const;
    double CalculateEnergyDissipationTermSpecialization(TElementData &rData) const;
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
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TwoFluidNavierStokesAlphaMethod& operator=(TwoFluidNavierStokesAlphaMethod const& rOther);

    /// Copy constructor.
    TwoFluidNavierStokesAlphaMethod(TwoFluidNavierStokesAlphaMethod const& rOther);

    ///@}

}; // Class TwoFluidNavierStokesAlphaMethod
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< class TElementData >
inline std::istream& operator >> (std::istream& rIStream,
    TwoFluidNavierStokesAlphaMethod<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
    const TwoFluidNavierStokesAlphaMethod<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_TWO_FLUID_NAVIER_STOKES_ALPHA_METHOD