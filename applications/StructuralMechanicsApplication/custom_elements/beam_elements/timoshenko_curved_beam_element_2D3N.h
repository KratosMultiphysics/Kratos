// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
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
 * @class LinearTimoshenkoCurvedBeamElement2D3N
 * @ingroup StructuralMechanicsApplication
 * @brief This is the Timoshenko curved beam element of 3 nodes. Reference:
 * Connecting beams and continua: variational basis and mathematical analysis, Romero and Schenk, Meccanica, 2023
 * DOI: https://doi.org/10.1007/s11012-023-01702-0
 * 
 * Ordering of the nodes:      0 ------ 2 ------- 1
 * 
 * Quadratic interpolation of the curved geometry, displacements and rotation. Reduced integration to avoid shear locking
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearTimoshenkoCurvedBeamElement2D3N
    : public Element
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = Element;

    static constexpr SizeType SystemSize    = 9;
    static constexpr SizeType NumberOfNodes = 3;
    static constexpr SizeType DoFperNode    = 3;
    static constexpr SizeType StrainSize    = 3;

    using GlobalSizeVector = BoundedVector<double, 9>;
    using array_3 = array_1d<double, 3>;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LinearTimoshenkoCurvedBeamElement2D3N);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    LinearTimoshenkoCurvedBeamElement2D3N() = default;

    // Constructor using an array of nodes
    LinearTimoshenkoCurvedBeamElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry)
    {
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
    }

    // Constructor using an array of nodes with properties
    LinearTimoshenkoCurvedBeamElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId,pGeometry,pProperties)
    {
        // This is needed to prevent uninitialised integration method in inactive elements
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
    }

    // Copy constructor
    LinearTimoshenkoCurvedBeamElement2D3N(LinearTimoshenkoCurvedBeamElement2D3N const& rOther)
        : BaseType(rOther),
        mThisIntegrationMethod(rOther.mThisIntegrationMethod),
        mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LinearTimoshenkoCurvedBeamElement2D3N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LinearTimoshenkoCurvedBeamElement2D3N>(NewId, pGeom, pProperties);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns a 9 component vector including the values of the DoFs
     * in GLOBAL beam axes
     */
    void GetNodalValuesVector(GlobalSizeVector &rNodalValues) const;

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
     * @brief Returns a 3 component vector with the values of the shape
     * functions at each node
     * xi: isoparametric coordinate
     */
    array_3 GetShapeFunctionsValues(const double xi) const;

    /**
     * @brief Returns a 3 component vector with the values of the shape
     * functions derivatives in the real space at each node
     * xi: isoparametric coordinate
     * J: Jacobian
     */
    array_3 GetFirstDerivativesShapeFunctionsValues(const double xi, const double J) const;

    /**
     * @brief Returns a 3 component vector with the values of the shape
     * functions second derivatives in the real space at each node
     * xi: isoparametric coordinate
     * J: Jacobian
     */
    array_3 GetSecondDerivativesShapeFunctionsValues(const double xi, const double J) const;

    /**
     * @brief Returns a 3 component vector with the values of the shape
     * functions derivatives in the natural space at each node
     * xi: isoparametric coordinate
     */
    array_3 GetLocalFirstDerivativesShapeFunctionsValues(const double xi) const;

    /**
     * @brief Returns a 3 component vector with the values of the shape
     * functions second derivatives in the natural space at each node
     * xi: isoparametric coordinate
     */
    array_3 GetLocalSecondDerivativesShapeFunctionsValues(const double xi) const;

    /**
     * @brief This method fills the global vectors of shape functions
     * in shuch a way that:
     * u = Nu * U
     * v = N * U
     * theta = Ntheta * U
     */
    void GetShapeFunctionsValuesGlobalVectors(
        const array_3 &rShapeFunctions,
        GlobalSizeVector &rNshape,
        GlobalSizeVector &rNu,
        GlobalSizeVector &rNtheta) const;

    /**
     * @brief Returns the Jacobian of the isoparametric transformation of arc length s
     *     J = sqrt((dx)^2 + (dy)^2)
     */
    double GetJacobian(const double xi) const;

    /**
     * @brief This function returns tangent and transverse unit vectors of the beam at coordinate xi
    */
    void GetTangentandTransverseUnitVectors(
        const double xi,
        array_3 &rt,
        array_3 &rn) const;

    /**
     * @brief This function builds the Frenet Serret matrix that rotates from global to local axes
    */
    BoundedMatrix<double, 2, 2> GetFrenetSerretMatrix(
        const double xi,
        const array_3 &rt,
        const array_3 &rn) const;

    /**
     * @brief This function retrieves the body forces in local axes
     * @param rElement the element reference
     * @param rIntegrationPoints array of IP
     * @param PointNumber tthe IP to be evaluated
    */
    array_3 GetBodyForce(
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
        rOStream << "Timoshenko 3N curved Beam Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
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

    IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods

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

}; // class LinearTimoshenkoCurvedBeamElement2D3N.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
