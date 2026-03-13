// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Marcelo Raschi
//                   Manuel Caicedo
//                   Javier Mroginski

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "base_solid_element.h"
#include "includes/variables.h"

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
 * @class SmallDisplacementBbar
 * @ingroup StructuralMechanicsApplication
 * @brief Infinitesimal strain definition with mixed B-bar formulation
 * @details Implements an infinitesimal strain definition with mixed B-bar
formulation, that avoids volumetric locking in elastic-plastic response
in nearly incompressible deformations (Hughes TJR, "The Finite Element Method:
Linear Static and Dynamic Finite Element Analysis" 1st Edition, section 4.5.2,
page 232)
 * @author Marcelo Raschi
 * @author Manuel Caicedo
 * @author Javier Mroginski
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallDisplacementBbar
        : public BaseSolidElement
{
protected:
/**
* Internal variables used in the kinematic calculations
*/
struct KinematicVariablesBbar
    : public KinematicVariables
{
    Vector Bh;

    /**
    * The default constructor
    * @param StrainSize The size of the strain vector in Voigt notation
    * @param Dimension The problem dimension: 2D or 3D
    * @param NumberOfNodes The number of nodes in the element
    */
    KinematicVariablesBbar(
        const SizeType StrainSize,
        const SizeType Dimension,
        const SizeType NumberOfNodes
        ) : KinematicVariables(StrainSize, Dimension, NumberOfNodes)
    {
        Bh = ZeroVector(Dimension * NumberOfNodes);
    }
};

public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// The base element type
    typedef BaseSolidElement BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of SmallDisplacementStrElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementBbar);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SmallDisplacementBbar(IndexType NewId, GeometryType::Pointer pGeometry);
    SmallDisplacementBbar(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    SmallDisplacementBbar(SmallDisplacementBbar const& rOther)
        :BaseType(rOther)
    {};


    /// Destructor.
    ~SmallDisplacementBbar() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param pGeom The pointer to the geometry of the element
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param ThisNodes The array containing nodes
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
    * Calculate a Matrix Variable on the Element Constitutive Law
    * @param rVariable The variable we want to get
    * @param rOutput The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateOnIntegrationPoints(
            const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rOutput,
            const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * Calculate a Vector Variable on the Element Constitutive Law
    * @param rVariable The variable we want to get
    * @param rOutput The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateOnIntegrationPoints(
            const Variable<Vector>& rVariable,
            std::vector<Vector>& rOutput,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

    /**
    * Calculate a double Variable on the Element Constitutive Law
    * @param rVariable The variable we want to get
    * @param rOutput The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateOnIntegrationPoints(
            const Variable<double>& rVariable,
            std::vector<double>& rOutput,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

    /**
    * Called at the end of each solution step
    * @param rCurrentProcessInfo the current process info instance
    */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    SmallDisplacementBbar() : BaseSolidElement()
    {
    }

    /**
    * @brief This method returns if the element provides the strain
    */
    bool UseElementProvidedStrain() const override;

    /**
    * This functions updates the kinematics variables
    * @param rThisKinematicVariables The kinematic variables to be calculated
    * @param PointNumber The integration point considered
    */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        ) override;

    /**
    * Calculation of the RHS
    */
    void CalculateAndAddResidualVector(
        VectorType& rRightHandSideVector,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const array_1d<double, 3>& rBodyForce,
        const Vector& rStressVector,
        const double IntegrationWeight
        ) const override;

    /**
     * This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     * @param ThisStressMeasure The stress measure considered
     */
    void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        ) override;


    /**
    * This functions calculates both the RHS and the LHS
    * @param rLeftHandSideMatrix The LHS
    * @param rRightHandSideVector The RHS
    * @param rCurrentProcessInfo The current process info instance
    * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
    * @param CalculateResidualVectorFlag The flag to set if compute the RHS
    */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    /**
    * Calculation of the Deformation Matrix B
    * @param B The deformation matrix
    * @param DN_DX The derivatives of the shape functions
    */
    void CalculateB(
        Matrix& rB,
        const Matrix& DN_DX
    );

    Matrix ComputeEquivalentF(const Vector& rStrainTensor);

    /**
        * This functions updates the kinematics variables
        * @param rThisKinematicVariables The kinematic variables to be calculated
        * @param PointNumber The integration point considered
        */
    void CalculateKinematicVariablesBbar(
        KinematicVariablesBbar& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        );

    /**
    * Calculation of the Deformation Matrix Bbar
    * @param B The deformation matrix
    * @param DN_DX The derivatives of the shape functions
    */
    void CalculateBbar(
        Matrix &rB,
        Vector &rBh,
        const Matrix &DN_DX,
        const GeometryType::IntegrationPointsArrayType &IntegrationPoints,
        const IndexType PointNumber
    );

    // Compute Bbar components
    /**
    * This functions updates the kinematics variables
    * @param rThisKinematicVariables The kinematic variables to be calculated
    */
    void CalculateHydrostaticDeformationMatrix(KinematicVariablesBbar& rThisKinematicVariables);

private:

    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //SmallDisplacementStrElement& operator=(const SmallDisplacementStrElement& rOther);
    /// Copy constructor.
    //SmallDisplacementStrElement(const SmallDisplacementStrElement& rOther);
    ///@}

}; // Class SmallDisplacementBbar

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
