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
 * @class LinearTrussElement2D
 * @ingroup StructuralMechanicsApplication
 * @brief This is the Linear TRUSS element of 2 nodes.
 * 2 and 3 noded elements
 * @author Alejandro Cornejo
 */
template <SizeType TNNodes>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearTrussElement2D
    : public Element
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = Element;
    static constexpr SizeType DofsPerNode = 2;
    static constexpr SizeType SystemSize = DofsPerNode * TNNodes;
    using SystemSizeBoundedArrayType = array_1d<double, SystemSize>;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LinearTrussElement2D);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    LinearTrussElement2D()
    {
    }

    // Constructor using an array of nodes
    LinearTrussElement2D(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry)
    {
        if constexpr (TNNodes == 2) {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
        } else {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
    }

    // Constructor using an array of nodes with properties
    LinearTrussElement2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId,pGeometry,pProperties)
    {
        // This is needed to prevent uninitialised integration method in inactive elements
        if constexpr (TNNodes == 2) {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
        } else {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
    }

    // Copy constructor
    LinearTrussElement2D(LinearTrussElement2D const& rOther)
        : BaseType(rOther),
        mThisIntegrationMethod(rOther.mThisIntegrationMethod),
        mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LinearTrussElement2D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override
    {
        return Kratos::make_intrusive<LinearTrussElement2D>(NewId, pGeom, pProperties);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Indicates the amount of DoFs per node (u0, v, theta)
     */
    IndexType GetDoFsPerNode() const
    {
        return 2;
    }

    /**
     * @brief This method returns the angle of the FE axis
     */
    double GetAngle() const
    {
        // return StructuralMechanicsElementUtilities::GetReferenceRotationAngle2D2NBeam(GetGeometry());
        return 1.0;
    }

    /**
     * @brief Returns a 6 component vector including the values of the DoFs
     * in LOCAL beam axes
     */
    void GetNodalValuesVector(VectorType& rNodalValue) const;

    /**
     * @brief Computes the length of the FE and returns it
     */
    double CalculateLength() const
    {
        return StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
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
     * @brief This function returns the 4 shape functions used for interpolating the transverse displacement v. (denoted as N)
     * Also its derivatives
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const;
    void GetFirstDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const;

    /**
     * @brief This function rotates the LHS from local to global coordinates
     * @param rLHS the left hand side
     * @param rGeometry the geometry of the FE
    */
    virtual void RotateLHS(
        MatrixType &rLHS,
        const GeometryType &rGeometry);

    /**
     * @brief This function rotates the RHS from local to global coordinates
     * @param rRHS the right hand side
     * @param rGeometry the geometry of the FE
    */
    virtual void RotateRHS(
        VectorType &rRHS,
        const GeometryType &rGeometry);

    /**
     * @brief This function rotates the LHS and RHS from local to global coordinates
     * @param rLHS the left hand side
     * @param rRHS the right hand side
     * @param rGeometry the geometry of the FE
    */
    virtual void RotateAll(
        MatrixType &rLHS,
        VectorType &rRHS,
        const GeometryType &rGeometry);

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
        rOStream << "Truss 2D Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
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

    IntegrationMethod mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2; /// Currently selected integration methods

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

}; // class LinearTrussElement2D.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
