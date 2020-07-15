//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Daniel Diez
//  Co-authors:      Ruben Zorrilla
//

#if !defined(KRATOS_TWO_FLUID_NAVIER_STOKES)
#define  KRATOS_TWO_FLUID_NAVIER_STOKES

// System includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "utilities/geometry_utilities.h"
#include "includes/deprecated_variables.h"

namespace Kratos
{

/*The "TwoFluidNavierStokes" element is an element based on the Variational Multiscale Stabilization technique (VMS)
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
class TwoFluidNavierStokes : public FluidElement<TElementData>
{
public:

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TwoFluidNavierStokes);

    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef std::vector< Dof<double>::Pointer > DofsVectorType;
    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;
    typedef typename FluidElement<TElementData>::ShapeFunctionsType ShapeFunctionsType;
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesType ShapeFunctionDerivativesType;
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesArrayType ShapeFunctionDerivativesArrayType;
    constexpr static unsigned int Dim = FluidElement<TElementData>::Dim;
    constexpr static unsigned int NumNodes = FluidElement<TElementData>::NumNodes;

    /* static constexpr int TNumNodes = NumNodes;
    static constexpr int TDim = Dim; */
    
    constexpr static unsigned int BlockSize = FluidElement<TElementData>::BlockSize;
    constexpr static unsigned int LocalSize = FluidElement<TElementData>::LocalSize;
    constexpr static unsigned int StrainSize = (Dim - 1) * 3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constuctor.
    /**
    * @param NewId Index number of the new element (optional)
    */
    TwoFluidNavierStokes(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
    * @param NewId Index of the new element
    * @param ThisNodes An array containing the nodes of the new element
    */
    TwoFluidNavierStokes(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    */
    TwoFluidNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    * @param pProperties Pointer to the element's properties
    */
    TwoFluidNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~TwoFluidNavierStokes();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
    * Returns a pointer to a new TwoFluidNavierStokes element, created using given input.
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

    /// Computes the elemental LHS and RHS elemental contributions
    /**
     * Given a distance function, computes the time integrated Left Hand Side (LHS)
     * and Right Hand Side elemental contributions for the two-fluid element.
     * @param rLeftHandSideMatrix elemental stiffness matrix
     * @param rRightHandSideVector elemental residual vector
     * @param rCurrentProcessInfo reference to the current process info
     */
    void CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo) override;

    /// Computes the elemental RHS elemental contribution
    /**
     * Given a distance function, computes the time integrated Right Hand Side (RHS)
     * elemental contribution for the two-fluid element.
     * @param rRightHandSideVector elemental residual vector
     * @param rCurrentProcessInfo reference to the current process info
     */
    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo) override;


    /**
     * @brief MassMatrix Calculate the local mass matrix.
     * @param rFluidStress Viscous stress in the fluid given in Voigt notation
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void Calculate( const Variable<Vector>& rVariable,
                    Vector& rOutput,
                    const ProcessInfo& rCurrentProcessInfo) override;

    /// Auxiliar element check function
    /**
     * This function calls the base element check method and adds the
     * current element check implementations
     * @param rCurrentProcessInfo reference to the current process info
     */
    int Check(const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Function to visualize the divergence field and the contact angle

    /**
     * @brief Get the Value On Integration Points object (used to visualize the divergence field and the contact angle)
     *
     * @param rVariable Variable to be retrieved (implementation supports DIVERGENCE and CONTACT_ANGLE)
     * @param rValues Vector for the values at the Gauss integration points
     * @param rCurrentProcessInfo ProcessInfo object
     */
    void GetValueOnIntegrationPoints(   const Variable<double> &rVariable,
                                        std::vector<double> &rValues,
                                        const ProcessInfo &rCurrentProcessInfo ) override;

    /**
     * @brief Calculate the Value On Integration Points object (used to visualize the contact angle)
     * CalculateOnIntegrationPoints seems to be deprecated!
     *
     * @param rVariable Variable to be retrieved (implementation supports CONTACT_ANGLE)
     * @param rOutput Vector for the values at the Gauss integration points
     * @param rCurrentProcessInfo ProcessInfo object
     */

    void CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo ) override;

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
     * @brief Computes time integrated LHS and RHS arrays
     * This method computes both the Left Hand Side and
     * Right Hand Side time integrated contributions.
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    void AddTimeIntegratedSystem(
        TElementData& rData,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    /**
     * @brief Computes the time integrated LHS matrix
     * This method computes the Left Hand Side time integrated contribution
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     */
    void AddTimeIntegratedLHS(
        TElementData& rData,
        MatrixType& rLHS) override;

    /**
     * @brief Computes the time integrated RHS vector
     * This method computes the Right Hand Side time integrated contribution
     * @param rData Reference to the element data container
     * @param rRHS Reference to the Right Hand Side matrix to be filled
     */
    void AddTimeIntegratedRHS(
        TElementData& rData,
        VectorType& rRHS) override;

    /**
     * @brief Computes the LHS Gauss pt. contribution
     * This method computes the contribution to the LHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     */
    void ComputeGaussPointLHSContribution(
        TElementData& rData,
        MatrixType& rLHS);

    /**
     * @brief Computes the RHS Gaus  pt. contribution
     * This method computes the contribution to the RHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    void ComputeGaussPointRHSContribution(
        TElementData& rData,
        VectorType& rRHS);

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
		VectorType& rRHS_ee);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Computes the LHS Gauss pt. contribution
     * This method computes the contribution to the LHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     */
    void ComputeGaussPointLHSContributionCut(
        TElementData& rData,
        MatrixType& rLHS);

    /**
     * @brief Computes the RHS Gaus  pt. contribution
     * This method computes the contribution to the RHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    void ComputeGaussPointRHSContributionCut(
        TElementData& rData,
        VectorType& rRHS);

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
	void ComputeGaussPointEnrichmentContributionsCut(
		TElementData& rData,
		MatrixType& rV,
		MatrixType& rH,
		MatrixType& rKee,
		VectorType& rRHS_ee);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /// Set up the element's data and constitutive law for the current integration point.
    /** @param[in/out] rData Container for the current element's data.
     *  @param[in] Weight Integration point weight.
     *  @param[in] rN Values of nodal shape functions at the integration point.
     *  @param[in] rDN_DX Values of nodal shape function gradients at the integration point.
     */
    void UpdateIntegrationPointData(
        TElementData& rData,
        unsigned int IntegrationPointIndex,
        double Weight,
        const typename TElementData::MatrixRowType& rN,
        const typename TElementData::ShapeDerivativesType& rDN_DX) const override;

    /// Set up the element's data for a cut element and constitutive law for the current integration point.
    /** @param[in/out] rData Container for the current element's data.
     *  @param[in] Weight Integration point weight.
     *  @param[in] rN Values of nodal shape functions at the integration point.
     *  @param[in] rDN_DX Values of nodal shape function gradients at the integration point.
     *  @param[in] rNenr Values of nodal enriched shape functions at the integration point.
     *  @param[in] rDN_DXenr Values of nodal enriched shape functions gradients at the integration point.
     */
    void UpdateIntegrationPointData(
        TElementData& rData,
        unsigned int IntegrationPointIndex,
        double Weight,
        const typename TElementData::MatrixRowType& rN,
        const typename TElementData::ShapeDerivativesType& rDN_DX,
        const typename TElementData::MatrixRowType& rNenr,
        const typename TElementData::ShapeDerivativesType& rDN_DXenr) const;

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

    /**
     * @brief Split shape functions computation auxiliar method
     * This method computes the standard and enrichment shape functions for a split element
     * @param rData Element data container
     * @param rShapeFunctionsPos Positive side shape functions values
     * @param rShapeFunctionsNeg Negative side shape functions values
     * @param rEnrichedShapeFunctionsPos Positive side enrichment shape functions values
     * @param rEnrichedShapeFunctionsNeg Negative side enrichment shape functions values
     * @param rShapeDerivativesPos  Positive side shape functions derivatives values
     * @param rShapeDerivativesNeg  Negative side shape functions derivatives values
     * @param rEnrichedShapeDerivativesPos Positive side enrichment shape functions derivatives values
     * @param rEnrichedShapeDerivativesNeg Negative side enrichment shape functions derivatives values
     */
    void ComputeSplitting(
		TElementData& rData,
		MatrixType& rShapeFunctionsPos,
        MatrixType& rShapeFunctionsNeg,
        MatrixType& rEnrichedShapeFunctionsPos,
        MatrixType& rEnrichedShapeFunctionsNeg,
        GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesNeg,
        GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesNeg);

    /**
     * @brief Split shape functions computation auxiliar method
     * This method computes the standard and enrichment shape functions for a split element and interfaces
     * @param rData Element data container
     * @param rShapeFunctionsPos Positive side shape functions values
     * @param rShapeFunctionsNeg Negative side shape functions values
     * @param rEnrichedShapeFunctionsPos Positive side enrichment shape functions values
     * @param rEnrichedShapeFunctionsNeg Negative side enrichment shape functions values
     * @param rShapeDerivativesPos  Positive side shape functions derivatives values
     * @param rShapeDerivativesNeg  Negative side shape functions derivatives values
     * @param rEnrichedShapeDerivativesPos Positive side enrichment shape functions derivatives values
     * @param rEnrichedShapeDerivativesNeg Negative side enrichment shape functions derivatives values
     * @param rInterfaceShapeDerivativesNeg Negative side shape functions derivatives at the interface-gauss-points
     * @param rInterfaceWeightsNeg Negative side weights for the interface-gauss-points
     * @param rInterfaceNormalsNeg Negative side normal vectors for the interface-gauss-points
     */
    void ComputeSplitting(
		TElementData& rData,
		MatrixType& rShapeFunctionsPos,
        MatrixType& rShapeFunctionsNeg,
        MatrixType& rEnrichedShapeFunctionsPos,
        MatrixType& rEnrichedShapeFunctionsNeg,
        GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesNeg,
        GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesNeg,
        MatrixType& rInterfaceShapeFunctionNeg,
        MatrixType& rEnrInterfaceShapeFunctionPos,
        MatrixType& rEnrInterfaceShapeFunctionNeg,
        GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        Kratos::Vector& rInterfaceWeightsNeg,
        std::vector<Vector>& rInterfaceNormalsNeg);

    /**
     * @brief Split shape functions computation auxiliar method
     * This method computes the standard and enrichment shape functions for a split element and interfaces
     * @param rData Element data container
     * @param rShapeFunctionsPos Positive side shape functions values
     * @param rShapeFunctionsNeg Negative side shape functions values
     * @param rEnrichedShapeFunctionsPos Positive side enrichment shape functions values
     * @param rEnrichedShapeFunctionsNeg Negative side enrichment shape functions values
     * @param rShapeDerivativesPos  Positive side shape functions derivatives values
     * @param rShapeDerivativesNeg  Negative side shape functions derivatives values
     * @param rEnrichedShapeDerivativesPos Positive side enrichment shape functions derivatives values
     * @param rEnrichedShapeDerivativesNeg Negative side enrichment shape functions derivatives values
     * @param rInterfaceShapeDerivativesNeg Negative side shape functions derivatives at the interface-gauss-points
     * @param rInterfaceWeightsNeg Negative side weights for the interface-gauss-points
     * @param rInterfaceNormalsNeg Negative side normal vectors for the interface-gauss-points
     * @param rContactShapeFunctionsNeg Negative side shape functions at the contact-line gauss-points: vector type for multiple C.L.
     * @param rContactShapeDerivativesNeg Negative side shape functions derivatives at the contact-line gauss-: vector type for multiple C.L.
     * @param rContactWeightsNeg Negative side weights for the contact-line gauss-points: vector type for multiple C.L.
     * @param rContactTangentialsNeg Negative side tangential vectors for the contact-line gauss-points: vector type for multiple C.L.
     * @param rHasContactLine Boolean: DEPRECATED
     */
    void ComputeSplitting(
		TElementData& rData,
		MatrixType& rShapeFunctionsPos,
        MatrixType& rShapeFunctionsNeg,
        MatrixType& rEnrichedShapeFunctionsPos,
        MatrixType& rEnrichedShapeFunctionsNeg,
        GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesNeg,
        GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesPos,
        GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesNeg,
        MatrixType& rInterfaceShapeFunctionNeg,
        MatrixType& rEnrInterfaceShapeFunctionPos,
        MatrixType& rEnrInterfaceShapeFunctionNeg,
        GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        Kratos::Vector& rInterfaceWeightsNeg,
        std::vector<Vector>& rInterfaceNormalsNeg,
        std::vector<MatrixType>& rContactShapeFunctionNeg,
        std::vector<GeometryType::ShapeFunctionsGradientsType>& rContactShapeDerivativesNeg,
        std::vector<Kratos::Vector>& rContactWeightsNeg,
        std::vector<Vector>& rContactTangentialsNeg);
        //MatrixType& rContactShapeFunctionNeg,
        //GeometryType::ShapeFunctionsGradientsType& rContactShapeDerivativesNeg,
        //Kratos::Vector& rContactWeightsNeg,
        //Vector& rContactTangentialsNeg,
        //bool& rHasContactLine);

    /**
     * @brief Calculates curvature at the gauss points of the interface.
     * @param rInterfaceCurvature Vector containing curvature values at the gauss points
     * @param rInterfaceShapeDerivativesNeg Negative side shape functions derivatives at the interface-gauss-points
     */
    void CalculateCurvature(
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        Kratos::Vector& rInterfaceCurvature);

    /**
     * @brief Calculates curvature at the gauss points of the interface.
     * @param rInterfaceCurvature Vector containing curvature values at the gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     */
    void CalculateCurvature(
        const Matrix& rIntShapeFunctions,
        Kratos::Vector& rInterfaceCurvature);

    /**
     * @brief Calculates curvature at the gauss points of the interface.
     * @param rInterfaceCurvature Vector containing curvature values at the gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rInterfaceShapeDerivativesNeg Negative side shape functions derivatives at the interface-gauss-points
     */
    void CalculateCurvature(
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        Kratos::Vector& rInterfaceCurvature);

    /**
     * @brief Impose pressure discontinuity at the interface due to the surface tension
     * A penalty method is NOT needed and integration is done on the interface
     * @param coefficient surface tension coefficient
     * @param rCurvature curvature calculated at the interface gauss points
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rIntEnrShapeFunctionsPos Enriched Shape functions calculated at the interface gauss points (positive side)
     * @param rIntEnrShapeFunctionsNeg Enriched Shape functions calculated at the interface gauss points (negative side)
     * @param rKeeTot Pressure enrichment contribution related to pressure enrichment DOFs will be modified by penalty method
     * @param rRHSeeTot Right Hand Side vector associated to the pressure enrichment DOFs will be modified by surface tension contribution
     */
	void PressureDiscontinuity(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const Matrix& rIntEnrShapeFunctionsPos,
        const Matrix& rIntEnrShapeFunctionsNeg,
		MatrixType& rKeeTot,
		VectorType& rRHSeeTot);

    /**
     * @brief Impose pressure discontinuity at the edges due to the surface tension
     * @param coefficient surface tension coefficient
     * @param rKeeTot Pressure enrichment contribution related to pressure enrichment DOFs will be modified by penalty method
     * @param rRHSeeTot Right Hand Side vector associated to the pressure enrichment DOFs will be modified by surface tension contribution
     */
	void PressureDiscontinuity(
        const double coefficient,
		MatrixType& rKeeTot,
		VectorType& rRHSeeTot);

    /**
     * @brief Computes the LHS terms associated with the pressure stabilizations at cut
     * This method is mainly proposed for Nitsche-XFEM
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     * @param rV Contribution related to the pressure enrichment DOFs in the N-S standard equations
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    void PressureGradientStabilization(
        TElementData& rData,
        MatrixType& rLHS,
        MatrixType& rV,
        VectorType& rRHS);

    /**
     * @brief Computes the LHS terms associated with the pressure stabilizations at the interface
     * This method is derived from what is done for Nitsche-XFEM
     * @param rIntWeights Weights associated with interface gauss points
     * @param rInterfaceShapeDerivativesNeg Shape functions derivatives at the interface-gauss-points
     * @param rKeeTot Pressure enrichment contribution related to pressure enrichment DOFs
     * @param rRHSeeTot Right Hand Side vector associated to the pressure enrichment DOFs
     */
    void PressureGradientStabilization(
        TElementData& rData,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntEnrShapeFunctionsPos,
        const Matrix& rIntEnrShapeFunctionsNeg,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        MatrixType& rKeeTot,
		VectorType& rRHSeeTot);

    /**
     * @brief Computes the LHS terms associated with the pressure stabilizations at the interface
     * This method is derived from what is done for Nitsche-XFEM
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntEnrShapeFunctionsPos (negative nodes) Enriched Shape functions calculated at the interface gauss points
     * @param rIntEnrShapeFunctionsNeg (positive nodes) Enriched Shape functions calculated at the interface gauss points
     * @param rInterfaceShapeDerivativesNeg Shape functions derivatives at the interface-gauss-points
     * @param rKeeTot Pressure enrichment contribution related to pressure enrichment DOFs
     */
    void PressureGradientStabilization(
        TElementData& rData,
        const Kratos::Vector& rIntWeights,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        MatrixType& rKeeTot);

    /**
     * @brief Computes the LHS terms associated with the pressure stabilizations
     * This method is derived from what is done in ghost penalty method
     * @param rWeights Weights associated with gauss points
     * @param rShapeDerivatives Shape functions derivatives at the gauss points
     * @param rKeeTot Pressure enrichment contribution related to pressure enrichment DOFs
     */
    void GhostPressureGradientStabilization(
        TElementData& rData,
        const Kratos::Vector& rWeights,
        const GeometryType::ShapeFunctionsGradientsType& rShapeDerivatives,
        MatrixType& rKeeTot);


    /**
     * @brief Impose pressure discontinuity at the interface and computes the surface tension
     * A penalty method is acquired and integration is done on the interface
     * @param coefficient surface tension coefficient
     * @param rCurvature curvature calculated at the interface gauss points
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntNormalsNeg Normal vectors (negative side) associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rIntEnrShapeFunctionsPos Enriched Shape functions calculated at the interface gauss points (positive side)
     * @param rIntEnrShapeFunctionsNeg Enriched Shape functions calculated at the interface gauss points (negative side)
     * @param rSurfaceTensionForce Surface tension force calculated by integrating \sigma \kappa n over the cut face
     * @param rKeeTot Pressure enrichment contribution related to pressure enrichment DOFs will be modified by penalty method
     * @param rRHSeeTot Right Hand Side vector associated to the pressure enrichment DOFs will be modified by surface tension contribution
     */
	void PressureDiscontinuityandSurfaceTension(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const std::vector<Vector>& rIntNormalsNeg,
        const Matrix& rIntShapeFunctions,
        const Matrix& rIntEnrShapeFunctionsPos,
        const Matrix& rIntEnrShapeFunctionsNeg,
        Vector& rSurfaceTensionForce,
		MatrixType& rKeeTot,
		VectorType& rRHSeeTot);    

    /**
     * @brief Computes the surface tension on the interface
     * @param coefficient surface tension coefficient
     * @param rCurvature curvature calculated at the interface gauss points
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntNormalsNeg Normal vectors (negative side) associated with interface gauss points
     * @param rSurfaceTensionForce Surface tension force calculated by integrating \sigma \kappa n over the cut face
     */
	void SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const std::vector<Vector>& rIntNormalsNeg,
        Vector& rSurfaceTensionForce);  

    /**
     * @brief Computes the surface tension on the interface and implement its effect on the RHS vector
     * @param coefficient surface tension coefficient
     * @param rCurvature curvature calculated at the interface gauss points
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rIntNormalsNeg Normal vectors (negative side) associated with interface gauss points
     * @param rRHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
	void SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const std::vector<Vector>& rIntNormalsNeg,
        VectorType& rRHS);  

    /**
     * @brief Computes the surface tension on the interface and implement its effect on the RHS vector
     * Added the effect of contact line for an open interface.
     * @param coefficient surface tension coefficient
     * @param rCurvature curvature calculated at the interface gauss points
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rIntNormalsNeg Normal vectors (negative side) associated with interface gauss points
     * @param rCLWeights Weights associated with contact line gauss points
     * @param rCLShapeFunctions Shape functions calculated at the contact line gauss points
     * @param rTangential Tangential vectors (according to negative side interfaces) associated with contact line gauss points
     * @param HasContactLine shows if there is a contact line
     * @param rRHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
    void SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const std::vector<Vector>& rIntNormalsNeg,
        const Kratos::Vector& rCLWeights,
        const Matrix& rCLShapeFunctions,
        const Vector& rTangential,
        bool HasContactLine,
        VectorType& rRHS);

   /**
     * @brief Computes the surface tension on the interface and implement its effect on the RHS vector
     * Added the effect of contact line for an open interface.
     * @param rData Element data container
     * @param coefficient surface tension coefficient
     * @param coefficientS solid surface contact net coefficient
     * @param zeta dissipative coefficient at the contact line
     * @param rCurvature curvature calculated at the interface gauss points
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rIntNormalsNeg Normal vectors (negative side) associated with interface gauss points
     * @param rCLWeights Weights associated with contact line gauss points
     * @param rCLShapeFunctions Shape functions calculated at the contact line gauss points
     * @param rTangential Tangential vectors (according to negative side interfaces) associated with contact line gauss points
     * @param HasContactLine shows if there is a contact : DEPRECATED
     * @param rLHS The contribution of contact line dissipative force to LHS
     * @param rRHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
    void SurfaceTension(
        const TElementData& rData,
        const double coefficient,
        const double coefficientS,
        const double zeta,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const std::vector<Vector>& rIntNormalsNeg,
        const std::vector<Kratos::Vector>& rCLWeights,
        const std::vector<Matrix>& rCLShapeFunctions,
        const std::vector<Vector>& rTangential,
        //bool HasContactLine,
        MatrixType& rLHS,
        VectorType& rRHS);

    /**
     * @brief Computes the surface tension on the interface and implement its effect on the RHS vector
     * Added the effect of contact line for an open interface.
     * @param rData Element data container
     * @param coefficient surface tension coefficient
     * @param coefficientS solid surface contact net coefficient
     * @param zeta dissipative coefficient at the contact line
     * @param micro_length_scale characteristic micro length-scale
     * @param rCurvature curvature calculated at the interface gauss points
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rIntNormalsNeg Normal vectors (negative side) associated with interface gauss points
     * @param rCLWeights Weights associated with contact line gauss points
     * @param rCLShapeFunctions Shape functions calculated at the contact line gauss points
     * @param rTangential Tangential vectors (according to negative side interfaces) associated with contact line gauss points
     * @param rLHS The contribution of contact line dissipative force to LHS
     * @param rRHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
    void SurfaceTension(
        const TElementData& rData,
        const double coefficient,
        const double coefficientS,
        const double zeta,
        const double micro_length_scale,
        const Kratos::Vector& rCurvature,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const std::vector<Vector>& rIntNormalsNeg,
        const std::vector<Kratos::Vector>& rCLWeights,
        const std::vector<Matrix>& rCLShapeFunctions,
        const std::vector<Vector>& rTangential,
        MatrixType& rLHS,
        VectorType& rRHS);

    /**
     * @brief Computes the surface tension on the interface and implement its effect on the RHS vector
     * Curvature is implicit in the formulation
     * @param coefficient surface tension coefficient
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rInterfaceShapeDerivativesNeg Negative side shape functions derivatives at the interface-gauss-points
     * @param rIntNormalsNeg Normal vectors (negative side) associated with interface gauss points
     * @param rRHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
	void SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        const std::vector<Vector>& rIntNormalsNeg,
        VectorType& rRHS); 

    /**
     * @brief Computes the surface tension on the interface and implement its effect on the RHS vector
     * Curvature is implicit in the formulation and normal vector is calculated from distance
     * @param coefficient surface tension coefficient
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rInterfaceShapeDerivativesNeg Negative side shape functions derivatives at the interface-gauss-points
     * @param rRHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
	void SurfaceTension(
        const double coefficient,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        VectorType& rRHS); 

    /**
     * @brief Computes the surface tension on the interface and implement its effect on the RHS vector
     * Curvature is implicit in the formulation and normal vector is calculated from distance
     * @param coefficient surface tension coefficient
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rInterfaceShapeDerivativesNeg Negative side shape functions derivatives at the interface-gauss-points
     * @param rRHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
	void SurfaceTensionMixed(
        const double coefficient,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        VectorType& rRHS); 

    /**
     * @brief Computes the pressure difference (surface tension) on the interface and implement its effect on the RHS vector
     * Curvature is implicit in the formulation and normal vector is calculated from distance
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rIntEnrShapeFunctionsPos Enriched Shape functions calculated at the interface gauss points (positive side)
     * @param rIntEnrShapeFunctionsNeg Enriched Shape functions calculated at the interface gauss points (negative side)
     * @param rVLHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
	void SurfaceTensionDP(
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const Matrix& rIntEnrShapeFunctionsPos,
        const Matrix& rIntEnrShapeFunctionsNeg,
        MatrixType& rVLHS); 

    /**
     * @brief Computes the surface tension on the interface and implement its effect on the RHS vector
     * Curvature is implicit in the formulation and normal vector is calculated from distance
     * LHS is also modified to include the contribution of surface tension
     * @param coefficient surface tension coefficient
     * @param rIntWeights Weights associated with interface gauss points
     * @param rIntShapeFunctions Shape functions calculated at the interface gauss points
     * @param rInterfaceShapeDerivativesNeg Negative side shape functions derivatives at the interface-gauss-points
     * @param rLHS The contribution of surface tension to LHS is included
     * @param rRHS The effect of pressure discontinuity is implemented as an interfacial integral on the RHS
     */
	void SurfaceTension(
        TElementData& rData,
        const double coefficient,
        const Kratos::Vector& rIntWeights,
        const Matrix& rIntShapeFunctions,
        const GeometryType::ShapeFunctionsGradientsType& rInterfaceShapeDerivativesNeg,
        MatrixType& rLHS,
        VectorType& rRHS); 

    /**
     * @brief Condense the enrichment
     * This method performs the static condensation of the enrichment terms, by adding
     * its local contributions to both the LHS and RHS elemental matrices.
     * @param rData Element data container
     * @param rLeftHandSideMatrix Reference to the element Left Hand Side matrix
     * @param rRightHandSideVector Reference to the element Right Hand Side vector
     * @param rVTot Common N-S equations term associated to pressure enrichment DOFs
     * @param rHTot Pressure enrichment contribution related to velocity and pressure DOFs
     * @param rKeeTot Pressure enrichment contribution related to pressure enrichment DOFs
     * @param rRHSeeTot Right Hand Side vector associated to the pressure enrichment DOFs
     */
	void CondenseEnrichment(
		const TElementData& rData,
		Matrix& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		const MatrixType& rVTot,
		const MatrixType& rHTot,
		MatrixType& rKeeTot,
		const VectorType& rRHSeeTot);

    /**
     * @brief Computes the contribution of the dissipation (Navier slip boundary condition) to the LHS
     * @param rData Element data container
     * @param betaIn Dissipation coefficient inside droplet at the contact surface
     * @param betaOut Dissipation coefficient outside droplet at the contact surface
     * @param betaContact Dissipation coefficient at the contact line
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     * @param rRHS RHS correction due to - LHS x {un}
     */
    void ContactSurfaceDissipation(
        const TElementData& rData,
        const double betaIn,
        const double betaOut,
        const double betaContact,
        MatrixType& rLHS,
        VectorType& rRHS);

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
    TwoFluidNavierStokes& operator=(TwoFluidNavierStokes const& rOther);

    /// Copy constructor.
    TwoFluidNavierStokes(TwoFluidNavierStokes const& rOther);

    ///@}

}; // Class TwoFluidNavierStokes
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< class TElementData >
inline std::istream& operator >> (std::istream& rIStream,
    TwoFluidNavierStokes<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
    const TwoFluidNavierStokes<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_TWO_FLUID_NAVIER_STOKES