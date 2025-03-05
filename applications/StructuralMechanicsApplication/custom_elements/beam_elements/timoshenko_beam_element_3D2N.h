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

#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "timoshenko_beam_element_2D2N.h"

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
 * @class LinearTimoshenkoBeamElement3D2N
 * @ingroup StructuralMechanicsApplication
 * @brief This is the 3D Timoshenko beam element of 2 nodes.
 * This element employs 3rd order locking-free Hermitic polynomials for the vertical deflection (v and w) that satisfy
 * the kinematic constraints of:
 *     theta_y = w' - E·Iy/GAsz * w'''
 *     theta_z = v' + E·Iz/GAsy * v'''
 * The longitudinal displacements (u) and torsional rotation (theta_x) are discretized wit standard linear shape functions
 * The elemental nodal unknown vector is:
 *     a = [u0, v0, w0, thetax_0, thetay_0, thetaz_0,   u1, v1, w1, thetax_1, thetay_1, thetaz_1]
 * The element is formulated in natural coordinates xi ranging [-1, 1]
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearTimoshenkoBeamElement3D2N
    : public LinearTimoshenkoBeamElement2D2N
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = LinearTimoshenkoBeamElement2D2N;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LinearTimoshenkoBeamElement3D2N);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    LinearTimoshenkoBeamElement3D2N()
    {
    }

    // Constructor using an array of nodes
    LinearTimoshenkoBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
    }

    // Constructor using an array of nodes with properties
    LinearTimoshenkoBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId,pGeometry,pProperties)
    {
        // This is needed to prevent uninitialised integration method in inactive elements
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
    }

    // Copy constructor
    LinearTimoshenkoBeamElement3D2N(LinearTimoshenkoBeamElement3D2N const& rOther)
        : BaseType(rOther)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LinearTimoshenkoBeamElement3D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override
    {
        return Kratos::make_intrusive<LinearTimoshenkoBeamElement3D2N>(NewId, pGeom, pProperties);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Indicates the amount of DoFs per node (u, v, w, theta_x, theta_y, theta_z)
     */
    IndexType GetDoFsPerNode() const override
    {
        return 6;
    }

    /**
     * @brief Returns a 12 component vector including the values of the DoFs in LOCAL beam axes
     */
    void GetNodalValuesVector(VectorType& rNodalValue) const override;

    /**
     * @brief Computes:
     * Axial strain:        E_l = du/dx
     * Torsional curvature: k_x = d theta_x / dx
     * Bending curvature y: k_y = d theta_y / dx
     * Bending curvature z: k_y = d theta_z / dx
     * Shear angle xy:      phi_y = dv/dx - theta_z
     * Shear angle xz:      phi_z = dv/dx + theta_z
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
     * @param rNodalValues The vector containing the nodal values in local axes
     */
    double         CalculateAxialStrain(  const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const override;
    virtual double CalculateShearStrainXY(const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const;
    virtual double CalculateShearStrainXZ(const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const;
    virtual double CalculateBendingCurvatureX(const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const;
    virtual double CalculateBendingCurvatureY(const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const;
    virtual double CalculateBendingCurvatureZ(const double Length, const double Phi, const double xi, const VectorType& rNodalValues) const;

    /**
     * @brief Computes
     *      Axial strain:        E_l = du/dx
     *      Torsional curvature: k_x = d theta_x / dx
     *      Bending curvature y: k_y = d theta_y / dx
     *      Bending curvature z: k_z = d theta_z / dx
     *      Shear angle xy:      phi_y = dv/dx - theta_z
     *      Shear angle xz:      phi_z = dv/dx + theta_z
     * and stores them in rStrain
     * @param rStrain The strain vector (6 components)
     * @param Length The size of the beam element
     * @param Phi The shear slenderness parameter
     * @param xi The coordinate in the natural axes
     * @param rNodalValues The vector containing the nodal values in local axes
     */
    void CalculateGeneralizedStrainsVector(VectorType& rStrain, const double Length, const double Phi, const double xi, const VectorType &rNodalValues) const override;

    /**
     * @brief Computes the length of the FE and returns it
     */
    double CalculateLength() const override
    {
        return StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    }

    /**
     * @brief Modifies a vector to include the components of a local size vector to the global size
     * @param rGlobalSizeVector The global size vector multiplying v and theta_z components
     * @param rLocalSizeVector The 4 local components of Nv
     */
    virtual void GlobalSizeVectorTransversalY(VectorType& rGlobalSizeVector, const VectorType& rLocalSizeVector) const
    {
        rGlobalSizeVector.clear();
        rGlobalSizeVector[1]  = rLocalSizeVector[0];
        rGlobalSizeVector[5]  = rLocalSizeVector[1];
        rGlobalSizeVector[7]  = rLocalSizeVector[2];
        rGlobalSizeVector[11] = rLocalSizeVector[3];
    }

    /**
     * @brief Modifies a vector to include the components of a local size vector to the global size
     * @param rGlobalSizeVector The global size vector multiplying w and theta_y components
     * @param rLocalSizeVector The 4 local components of Nw
     */
    virtual void GlobalSizeVectorTransversalZ(VectorType& rGlobalSizeVector, const VectorType& rLocalSizeVector) const
    {
        rGlobalSizeVector.clear();
        rGlobalSizeVector[2]  = rLocalSizeVector[0];
        rGlobalSizeVector[4]  = rLocalSizeVector[1];
        rGlobalSizeVector[8]  = rLocalSizeVector[2];
        rGlobalSizeVector[10] = rLocalSizeVector[3];
    }

    /**
     * @brief Modifies a vector to include the components of a local size vector to the global size
     * @param rGlobalSizeVector The global size vector including only the axial u terms
     * @param rLocalSizeVector The 2 local components of u
     */
    void GlobalSizeAxialVector(VectorType& rGlobalSizeVector, const VectorType& rLocalSizeVector) override
    {
        rGlobalSizeVector.clear();
        rGlobalSizeVector[0] = rLocalSizeVector[0];
        rGlobalSizeVector[6] = rLocalSizeVector[1];
    }

    /**
     * @brief Modifies a vector to include the components of a local size vector to the global size
     * @param rGlobalSizeVector The global size vector including only the axial theta_x terms
     * @param rLocalSizeVector The 2 local components of theta_x
     */
    virtual void GlobalSizeVectorAxialRotation(
        VectorType& rGlobalSizeVector,
        const VectorType& rLocalSizeVector)
    {
        rGlobalSizeVector.clear();
        rGlobalSizeVector[3] = rLocalSizeVector[0];
        rGlobalSizeVector[9] = rLocalSizeVector[1];
    }

    /**
     * @brief Assembles the 3 dimension rotation matrix to a global element size one
     * @param rT The local 3 dimension rotation matrix
     * @param rGlobalT The global 12 dimension rotation matrix
     */
    virtual void AssembleGlobalRotationMatrix(
        const BoundedMatrix<double, 3, 3>& rT,
        BoundedMatrix<double, 12, 12>& rGlobalT);

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
     * @brief This function rotates the LHS from local to global coordinates
     * @param rLHS the left hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateLHS(
        MatrixType &rLHS,
        const GeometryType &rGeometry) override;

    /**
     * @brief This function rotates the RHS from local to global coordinates
     * @param rRHS the right hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateRHS(
        VectorType &rRHS,
        const GeometryType &rGeometry) override;

    /**
     * @brief This function rotates the LHS and RHS from local to global coordinates
     * @param rLHS the left hand side
     * @param rRHS the right hand side
     * @param rGeometry the geometry of the FE
    */
    void RotateAll(
        MatrixType &rLHS,
        VectorType &rRHS,
        const GeometryType &rGeometry) override;

    /**
     * @brief This function retrieves the body forces in local axes
     * @param rElement the element reference
     * @param rIntegrationPoints array of IP
     * @param PointNumber tthe IP to be evaluated
    */
    array_1d<double, 3> GetLocalAxesBodyForce(
        const Element &rElement,
        const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
        const IndexType PointNumber) const override;

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
    // void CalculateOnIntegrationPoints(
    //     const Variable<double>& rVariable,
    //     std::vector<double>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo
    //     ) override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    // void CalculateOnIntegrationPoints(
    //     const Variable<Vector>& rVariable,
    //     std::vector<Vector>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo
    //     ) override;

    /**
     * @brief Get on rVariable Constitutive Law from the element
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    // void CalculateOnIntegrationPoints(
    //     const Variable<ConstitutiveLaw::Pointer>& rVariable,
    //     std::vector<ConstitutiveLaw::Pointer>& rValues,
    //     const ProcessInfo& rCurrentProcessInfo
    //     ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    // int Check(const ProcessInfo &rCurrentProcessInfo) const override;

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
        rOStream << "3D Timoshenko 2N Beam Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


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

}; // class LinearTimoshenkoBeamElement3D2N.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
