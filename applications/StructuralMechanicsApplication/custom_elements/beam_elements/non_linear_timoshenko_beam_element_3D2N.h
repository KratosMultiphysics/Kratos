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
#pragma once

#include "includes/element.h"

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
 * @class NonLinearTimoshenkoBeamElement3D2N
 * @ingroup StructuralMechanicsApplication
 * @brief This is the 3D Timoshenko beam element of 2 nodes. 
 * Reference:  I. Romero and F. Armero, "An objective finite element approximation o the kinetics of geometrically exact rods
 * and its use in the formulation of an energy-momentum conservative scheme in dynamics", IJNME, 2002.
 * DOI: 10.1002/nme.486 
 * @author Alejandro Cornejo
 */
class NonLinearTimoshenkoBeamElement3D2N : public Element
{
public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NonLinearTimoshenkoBeamElement3D2N);

    ///@name Type Definitions
    ///@{
    using BaseType = Element;

    using array3 = array_1d<double, 3>;

    NonLinearTimoshenkoBeamElement3D2N() {}

    NonLinearTimoshenkoBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {}

    NonLinearTimoshenkoBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {}

    NonLinearTimoshenkoBeamElement3D2N(NonLinearTimoshenkoBeamElement3D2N const& rOther)
        : BaseType(rOther)
    {}

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearTimoshenkoBeamElement3D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearTimoshenkoBeamElement3D2N>(NewId, pGeom, pProperties);
    }

    Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override;

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
     * @brief Calculate local system
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate left hand side
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate right hand side
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief General calculation method with flags
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLHS,
        const bool ComputeRHS
        );

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief It initializes the material
     */
    void InitializeMaterial();

    /**
     * @brief Indicates the amount of DoFs per node (u, v, w, theta_x, theta_y, theta_z)
     */
    IndexType GetDoFsPerNode() const
    {
        return 6; // 3 displacements and 3 rotations
    }

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

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
     * @brief Sets the used constitutive laws
     * @param ThisConstitutiveLawVector Constitutive laws used
     */
    void SetRotationOperators(const std::vector<BoundedMatrix<double, 3, 3>>& rThisRotationOperators)
    {
        mRotationOperators = rThisRotationOperators;
    }

    /**
     * @brief This method computes the DoF mapping operator defined in Eq. (37) from the ref.
     * NOTE the ref is incorrect. It has been corrected in here.
     * This matrix maps the (du, dd1, dd2, dd3) -> (du, dtheta), being "d" the variation
     */
    BoundedMatrix<double, 12, 6> CalculateDoFMappingMatrix(
        const Vector &rD1,
        const Vector &rD2,
        const Vector &rD3);

    /**
     * @brief This method computes the B matrix that relates the variation of the DoFs and the variation of the
     * strains, its size is a 6 by 12.
     */
    BoundedMatrix<double, 6, 12> CalculateB(
        const double N1,
        const double N2,
        const double dN1,
        const double dN2);

    /**
     * @brief This method computes the element length in the reference configuration
     */
    double CalculateReferenceLength() const;

    /**
     * @brief This method computes the generalized strain vector (6 components)
     * The ordering is of the components is the Kratos one:
     * generalized_strain = [Gamma, Omega], see Timoshenko beam constitutive law.
     */
    Vector CalculateStrainVector(
        const double N1,
        const double N2,
        const double dN1,
        const double dN2);

    /**
     * @brief This method computes the generalized stress vector (6 components) and/or the generalized constitutive matrix (6x6)
     * The ordering of the strain is supposed to be the Kratos one, however the stress and constitutive are returned as
     * Romero and Armero for coherence
     */
    void CalculateGeneralizedResponse(
        const IndexType IntegrationPoint,
        ConstitutiveLaw::Parameters rValues);

    /**
     * @brief This method builds the rotation operator in the reference configuration
     * This corresponds to the "Lambda" operator in Romero and Armero, Eq. (3)
     * Operator = [d1, d2, d3] as column vectors
     */
    BoundedMatrix<double, 3, 3> CalculateInitialRotationOperator(const bool UseCurrentConfiguration = false);

    /**
     * @brief this is called for non-linear analysis at the end of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

private:

    /* The rotation operators are built with [d1, d2, d3] as col vectors */
    std::vector<BoundedMatrix<double, 3, 3>> mRotationOperators; // The two rotation matrices, one per each IP.

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; // The vector containing the beam constitutive laws, one per each IP
    IntegrationMethod mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_LOBATTO_1; // By default the quadrature points are located at the nodes of the beam


};

} // namespace Kratos
