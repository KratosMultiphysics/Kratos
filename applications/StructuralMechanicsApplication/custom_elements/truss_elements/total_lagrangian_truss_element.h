// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

using SizeType = std::size_t;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class TotalLagrangianTrussElement
 * @ingroup StructuralMechanicsApplication
 * @brief This is the Total Lagrangian  TRUSS element of 2 nodes for 2D and 3D.
 * O---------O --> x'
 *  0         1
 * @author Alejandro Cornejo
 * Reference: "Non-linear Finite Element Analysis of Solids and Structures" by Belytschko, Liu, Moran.
 */
template <SizeType TDimension>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TotalLagrangianTrussElement
    : public Element
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = Element;
    static constexpr SizeType NNodes      = 2;
    static constexpr SizeType Dimension   = TDimension;
    static constexpr SizeType DofsPerNode = TDimension;
    static constexpr SizeType SystemSize  = TDimension * 2;
    using SystemSizeBoundedArrayType      = array_1d<double, SystemSize>;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TotalLagrangianTrussElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    TotalLagrangianTrussElement()
    {
    }

    // Constructor using an array of nodes
    TotalLagrangianTrussElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry)
    {
        mThisIntegrationMethod = CalculateIntegrationMethod();
    }

    // Constructor using an array of nodes with properties
    TotalLagrangianTrussElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId,pGeometry,pProperties)
    {
        mThisIntegrationMethod = CalculateIntegrationMethod();
    }

    // Copy constructor
    TotalLagrangianTrussElement(TotalLagrangianTrussElement const& rOther)
        : BaseType(rOther),
        mThisIntegrationMethod(rOther.mThisIntegrationMethod),
        mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TotalLagrangianTrussElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override
    {
        return Kratos::make_intrusive<TotalLagrangianTrussElement>(NewId, pGeom, pProperties);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns the angle of the FE axis in the reference configuration
     * 2D calculations
     */
    double GetAngle() const
    {
        return StructuralMechanicsElementUtilities::GetReferenceRotationAngle2D2NBeam(GetGeometry());
    }

    /**
     * @brief This function builds the Frenet Serret matrix that rotates from global to local axes in the reference configuration
     * T = | <- t -> |  x, local
     *     | <- n -> |  y, local
     *     | <- m -> |  z, local
     * 3D calculations
    */
    BoundedMatrix<double, 3, 3> GetFrenetSerretMatrix() const;

    /**
     * @brief This method returns the integration method depending on the Number of Nodes
     */
    IntegrationMethod CalculateIntegrationMethod()
    {
        return GeometryData::IntegrationMethod::GI_GAUSS_1;
    }

    /**
     * @brief This method returns the base shape functions in a reduced size vector
     */
    Vector GetBaseShapeFunctions(const double xi) const;

    /**
     * @brief Returns a n component vector including the values of the DoFs
     * in LOCAL truss axes
     */
    void GetNodalValuesVector(SystemSizeBoundedArrayType& rNodalValue) const;

    /**
     * @brief Computes the reference (initial) length of the FE and returns it
     */
    double CalculateReferenceLength() const
    {
        if constexpr (Dimension == 2) {
            // Same implementation for 2N and 3N
            return StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
        } else {
            return StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        }
        return 0.0; // TODO
    }

    /**
     * @brief Computes the current length (updated) of the FE and returns it
     */
    double CalculateCurrentLength() const
    {
        if constexpr (Dimension == 2) {
            // Same implementation for 2N and 3N
            return StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
        } else {
            return StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
        }
        return 0.0; // TODO
    }

    /**
     * @brief Computes the Green-Lagrange strain of the element
     * @param CurrentLength The current length of the element
     * @param ReferenceLength The reference length of the element
     * @return The Green-Lagrange strain
     */
    double CalculateGreenLagrangeStrain(const double CurrentLength, const double ReferenceLength) const
    {
        return 0.5 * (std::pow(CurrentLength / ReferenceLength, 2) - 1.0);
    }

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The vector containing the dof of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Returns the used integration method
     * @return default integration method of the used Geometry
     */
    IntegrationMethod GetIntegrationMethod() const override
    {
        return mThisIntegrationMethod;
    }

    /**
    * element can be integrated using the GP provided by the geometry or custom ones
    * by default, the base element will use the standard integration provided by the geom
    * @return bool to select if use/not use GPs given by the geometry
    */
    bool UseGeometryIntegrationMethod() const
    {
        return true;
    }

    /**
     * @brief Returns the set of integration points
     */
    const GeometryType::IntegrationPointsArrayType IntegrationPoints() const
    {
        return GetGeometry().IntegrationPoints();
    }

    /**
     * @brief Returns the set of integration points
     */
    const GeometryType::IntegrationPointsArrayType IntegrationPoints(IntegrationMethod ThisMethod) const
    {
        return GetGeometry().IntegrationPoints(ThisMethod);
    }

    /**
     * @brief This function returns the 2 shape functions used for interpolating the displacements in x,y,z
     * Also the derivatives of u to compute the longitudinal strain
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetShapeFunctionsValues(SystemSizeBoundedArrayType& rN, const double Length, const double xi) const;
    void GetShapeFunctionsValuesY(SystemSizeBoundedArrayType& rN, const double Length, const double xi) const;
    void GetShapeFunctionsValuesZ(SystemSizeBoundedArrayType& rN, const double Length, const double xi) const;
    void GetFirstDerivativesShapeFunctionsValues(SystemSizeBoundedArrayType& rN, const double Length, const double xi) const;

    /**
     * @brief This function rotates the LHS from local to global coordinates
     * @param rLHS the left hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateLHS(
        MatrixType &rLHS);

    /**
     * @brief This function rotates the RHS from local to global coordinates
     * @param rRHS the right hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateRHS(
        VectorType &rRHS);

    /**
     * @brief This function rotates the LHS and RHS from local to global coordinates
     * @param rLHS the left hand side
     * @param rRHS the right hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateAll(
        MatrixType &rLHS,
        VectorType &rRHS);

    /**
     * @brief This function retrieves the body forces in local axes
     * @param rElement the element reference
     * @param rIntegrationPoints array of IP
     * @param PointNumber tthe IP to be evaluated
    */
    array_1d<double, 3> GetLocalAxesBodyForce(
        const Element &rElement,
        const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
        const IndexType PointNumber) const;

    /**
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
     * @param rRightHandSideVector container for the desired RHS output
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector the elemental right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
 * @brief Calculate a Vector Variable on the Element Constitutive Law
 * @param rVariable The variable we want to get
 * @param rOutput The values obtained in the integration points
 * @param rCurrentProcessInfo the current process info instance
 */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable, std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Get on rVariable Constitutive Law from the element
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Truss Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    IntegrationMethod mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1; /// Currently selected integration methods
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Sets the used integration method
     * @param ThisIntegrationMethod Integration method used
     */
    void SetIntegrationMethod(const IntegrationMethod& rThisIntegrationMethod)
    {
        mThisIntegrationMethod = rThisIntegrationMethod;
    }

    /**
     * @brief Sets the used constitutive laws
     * @param ThisConstitutiveLawVector Constitutive laws used
     */
    void SetConstitutiveLawVector(const std::vector<ConstitutiveLaw::Pointer>& rThisConstitutiveLawVector)
    {
        mConstitutiveLawVector = rThisConstitutiveLawVector;
    }

    /**
     * @brief It initializes the material
     */
    void InitializeMaterial();


    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer &rSerializer) const override;

    void load(Serializer &rSerializer) override;

}; // class TotalLagrangianTrussElement.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
