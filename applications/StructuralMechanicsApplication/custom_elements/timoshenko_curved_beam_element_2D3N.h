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
 * @brief This is the Timoshenko curved beam element of 3 nodes. References:
 *  Felippa and Oñate,
 * "Accurate Timoshenko Beam Elements For Linear Elastostatics and LPB Stability",
 * Archives of Comp. Methods in Eng. (2021) 28:2021-2080
 * 
 * and
 * 
 * Hosseini et al.
 * Isogeometric analysis of free-form Timoshenko curved beams including the nonlinear effects of large deformations
 * Acta Mechanica Sinica/Lixue Xuebao
 * 
 * Ordering of the nodes:      0 ------ 2 ------- 1
 * 
 * Quadratic interpolation of the curved geometry and longitudinal displacement, u
 * Quintic interpolation of deflection v
 * Quartic interpolation of the total rotation Theta
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

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LinearTimoshenkoCurvedBeamElement2D3N);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    LinearTimoshenkoCurvedBeamElement2D3N()
    {
    }

    // Constructor using an array of nodes
    LinearTimoshenkoCurvedBeamElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry)
    {
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
    }

    // Constructor using an array of nodes with properties
    LinearTimoshenkoCurvedBeamElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId,pGeometry,pProperties)
    {
        // This is needed to prevent uninitialised integration method in inactive elements
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
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
     * @brief Returns a 6 component vector including the values of the DoFs
     * in GLOBAL beam axes
     */
    void GetNodalValuesVector(
        GlobalSizeVector &rNodalValues,
        const double angle1,
        const double angle2,
        const double angle3);

    /**
     * @brief Computes the axial strain (El), shear strain (gamma_xy) and bending curvature (kappa)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
     * @param rNodalValues The vector containing the nodal values in local axes
     */
    double CalculateAxialStrain     (const double J, const double xi, const GlobalSizeVector& rNodalValues);
    double CalculateShearStrain     (const double J, const double xi, const GlobalSizeVector& rNodalValues);
    double CalculateBendingCurvature(const double J, const double xi, const GlobalSizeVector& rNodalValues);

    /**
     * @brief Computes the axial strain (El), shear strain (gamma_xy) and bending curvature (kappa) and builds the strain vector
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
     * @param rNodalValues The vector containing the nodal values in local axes
     */
    void CalculateGeneralizedStrainsVector(VectorType& rStrain, const double J, const double xi, const GlobalSizeVector &rNodalValues);

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
     * @brief Returns the Jacobian of the isoparametric transformation
     *     J = sqrt((dx)^2 + (dy)^2)
     */
    const double GetJacobian(const double xi);

    /**
     * @brief Returns the bending/shear ratio stiffness
     */
    const double GetBendingShearStiffnessRatio();

    /**
     * @brief Returns the curvature of the geometry
     */
    const double GetGeometryCurvature(
        const double J,
        const double xi);

    /**
     * @brief This function returns the 4 shape functions used for interpolating the transverse displacement v. (denoted as N)
     * Also its derivatives
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetShapeFunctionsValues                 (GlobalSizeVector& rN, const double J, const double xi);
    void GetFirstDerivativesShapeFunctionsValues (GlobalSizeVector& rN, const double J, const double xi);
    void GetSecondDerivativesShapeFunctionsValues(GlobalSizeVector& rN, const double J, const double xi);
    void GetThirdDerivativesShapeFunctionsValues (GlobalSizeVector& rN, const double J, const double xi);
    void GetFourthDerivativesShapeFunctionsValues(GlobalSizeVector& rN, const double J, const double xi);



    BoundedMatrix<double, 9, 9> GetGlobalSizeRotationMatrixGlobalToLocalAxes(const double xi)
    {
        GlobalSizeVector dN_dxi, d2N_dxi2;
        GetLocalFirstDerivativesNu0ShapeFunctionsValues (dN_dxi,   xi);
        GetLocalSecondDerivativesNu0ShapeFunctionsValues(d2N_dxi2, xi);
        const auto& r_geom = GetGeometry();

        double dx_dxi = 0.0;
        double dy_dxi = 0.0;

        double d2x_dxi2 = 0.0;
        double d2y_dxi2 = 0.0;

        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            const IndexType u_coord = DoFperNode * i;
            const auto &r_coords_node = r_geom[i].GetInitialPosition();
            dx_dxi += r_coords_node[0] * dN_dxi[u_coord];
            dy_dxi += r_coords_node[1] * dN_dxi[u_coord];

            d2x_dxi2 += r_coords_node[0] * d2N_dxi2[u_coord];
            d2y_dxi2 += r_coords_node[1] * d2N_dxi2[u_coord];
        }

        VectorType x_prime(3), x_prime2(3), t(3), n(3), b(3);
        x_prime.clear();
        // x_prime2.clear();
        t.clear();
        n.clear();
        b.clear();
        b[2] = -1.0;

        x_prime[0] = dx_dxi;
        x_prime[1] = dy_dxi;
        x_prime[2] = 0.0;

        // x_prime2[0] = d2x_dxi2;
        // x_prime2[1] = d2y_dxi2;
        // x_prime2[2] = 0.0;

        noalias(t) = x_prime / norm_2(x_prime);
        noalias(n) = MathUtils<double>::CrossProduct(t, b);

        BoundedMatrix<double, 3, 3> T;
        T.clear();

        T(0, 0) = t[0];
        T(0, 1) = t[1];

        T(1, 0) = n[0];
        T(1, 1) = n[1];

        T(2, 2) = 1.0;

        BoundedMatrix<double, 9, 9> global_size_T;
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);
        return global_size_T;
    }

    void RotateLHS(
        MatrixType &rLHS,
        const double angle1,
        const double angle2,
        const double angle3);

    double GetAngle(const double xi);

    void RotateRHS(
        VectorType &rRHS,
        const double angle1,
        const double angle2,
        const double angle3);

    void RotateAll(
        MatrixType &rLHS,
        VectorType &rRHS,
        const double angle1,
        const double angle2,
        const double angle3);

    /**
     * @brief This function returns the 4 shape functions used for interpolating the total rotation Theta (N_theta)
     * Also its derivative
     * @param rN reference to the shape functions (or derivatives)
     * @param J The jacobian of the beam element at xi
     * @param ShearFactor The shear slenderness parameter
     * @param k0 The curvature of the geometry
     * @param xi The coordinate in the natural axes
    */
    void GetNThetaShapeFunctionsValues                (GlobalSizeVector& rN, const double J, const double xi);
    void GetFirstDerivativesNThetaShapeFunctionsValues(GlobalSizeVector& rN, const double J, const double xi);

    /**
     * @brief This function returns the 2 shape functions used for interpolating the axial displacement u0
     * Also its derivatives
     * @param rN reference to the shape functions (or derivatives)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
    */
    void GetNu0ShapeFunctionsValues                 (GlobalSizeVector& rNu,                const double xi);
    void GetFirstDerivativesNu0ShapeFunctionsValues (GlobalSizeVector& rNu, const double J, const double xi);
    void GetSecondDerivativesNu0ShapeFunctionsValues(GlobalSizeVector& rNu, const double J, const double xi);
    void GetLocalFirstDerivativesNu0ShapeFunctionsValues (GlobalSizeVector& rNu, const double xi);
    void GetLocalSecondDerivativesNu0ShapeFunctionsValues(GlobalSizeVector& rNu, const double xi);


    /**
     * @brief This function retrieves the body forces in local axes
     * @param rElement the element reference
     * @param rIntegrationPoints array of IP
     * @param PointNumber tthe IP to be evaluated
    */
    array_1d<double, 3> GetLocalAxesBodyForce(
        const Element &rElement,
        const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
        const IndexType PointNumber,
        const double angle);

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
