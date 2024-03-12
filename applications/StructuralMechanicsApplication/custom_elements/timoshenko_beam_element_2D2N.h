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
// #include "includes/define.h"
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
// #include "utilities/integration_utilities.h"
// #include "custom_utilities/structural_mechanics_element_utilities.h"

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
 * @class TimoshenkoBeamElement2D2N
 * @ingroup StructuralMechanicsApplication
 * @brief This is the Timoshenko beam element of 2 nodes. Reference: Felippa and Oñate,
 * "Accurate Timoshenko Beam Elements For Linear Elastostatics and LPB Stability",
 * Archives of Comp. Methods in Eng. (2021) 28:2021-2080
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TimoshenkoBeamElement2D2N
    : public Element
{

public:

    ///@name Type Definitions
    ///@{

    ///Type definition for integration methods
    // using IntegrationMethod = GeometryData::IntegrationMethod;

    /// The base element type
    using BaseType = Element;

    /// The definition of the sizetype
    // using IndexType = std::size_t;
    // using SizeType  = std::size_t;

    // using VectorType = BaseType::VectorType;

    // using MatrixType = BaseType::MatrixType;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TimoshenkoBeamElement2D2N);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    TimoshenkoBeamElement2D2N()
    {
    }

    // Constructor using an array of nodes
    TimoshenkoBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry){};

    // Constructor using an array of nodes with properties
    TimoshenkoBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId,pGeometry,pProperties)
    {
        // This is needed to prevent uninitialised integration method in inactive elements
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    }

    // Copy constructor
    TimoshenkoBeamElement2D2N(TimoshenkoBeamElement2D2N const& rOther)
        : BaseType(rOther),
        mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_intrusive<TimoshenkoBeamElement2D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_intrusive<TimoshenkoBeamElement2D2N>(NewId, pGeom, pProperties);
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
    const IndexType GetDoFsPerNode() const
    {
        return 3;
    }

    /**
     * @brief Returns a 6 component vector including the values of the DoFs
     */
    void GetNodalValuesVector(VectorType& rNodalValue);

    /**
     * @brief Computes the axial strain (El), shear strain (gamma_xy) and bending curvature (kappa)
     */
    const double CalculateAxialStrain     (const double Length, const double Phi, const double xi, const VectorType& rNodalValues);
    const double CalculateShearStrain     (const double Length, const double Phi, const double xi, const VectorType& rNodalValues);
    const double CalculateBendingCurvature(const double Length, const double Phi, const double xi, const VectorType& rNodalValues);

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

    const GeometryType::IntegrationPointsArrayType IntegrationPoints() const 
    {
        return GetGeometry().IntegrationPoints();
    }

    const GeometryType::IntegrationPointsArrayType IntegrationPoints(IntegrationMethod ThisMethod) const
    {
        return GetGeometry().IntegrationPoints(ThisMethod);
    }

    /**
     * @brief This function computes the Phi parameter in Felippa et al.
    */
    const double CalculatePhi(ConstitutiveLaw::Parameters &rValues);

    /**
     * @brief This function returns the 4 shape functions used for interpolating the transverse displacement v. (denoted as N)
     * Also its derivatives
    */
    void GetShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi);
    void GetFirstDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi);
    void GetSecondDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi);
    void GetThirdDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi);

    /**
     * @brief This function returns the 4 shape functions used for interpolating the total rotation Theta (N_theta)
     * Also its derivative
    */
    void GetNThetaShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi);
    void GetFirstDerivativesNThetaShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi);

    /**
     * @brief This function returns the 2 shape functions used for interpolating the axial displacement u0
     * Also its derivative
    */
    void GetNu0ShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi);
    void GetFirstDerivativesNu0ShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi);

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
    int Check( const ProcessInfo& rCurrentProcessInfo ) const override;

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
        rOStream << "Timoshenko Beam Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
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
    virtual void InitializeMaterial();

    /**
     * @brief This functions computes the integration weight to consider
     * @param ThisIntegrationMethod The array containing the integration points
     * @param PointNumber The id of the integration point considered
     * @param detJ The determinant of the jacobian of the element
     */
    double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
        const SizeType PointNumber,
        const double detJ
        ) const;

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

}; // class TimoshenkoBeamElement2D2N.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
