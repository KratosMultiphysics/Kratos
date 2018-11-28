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
    KRATOS_CLASS_POINTER_DEFINITION(TwoFluidNavierStokes);

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
     * @brief Computes the material response
     * This method computes the material response taking into account if the 
     * element is air or fluid. The material response is saved in the data container.
     * @param rData Reference to the element data container
     */
    void CalculateMaterialResponse(
        TElementData& rData) const override;

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