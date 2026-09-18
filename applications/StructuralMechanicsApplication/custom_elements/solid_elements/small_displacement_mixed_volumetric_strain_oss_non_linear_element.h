// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"

// Application includes
#include "structural_mechanics_application_variables.h"
#include "small_displacement_mixed_volumetric_strain_oss_element.h"

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
 * @class SmallDisplacementMixedVolumetricStrainOssNonLinearElement
 * @ingroup StructuralMechanicsApplication
 * @brief Small displacement with strain based mixed formulation element
 * @details This implements a small displacements element formulation with an extra volumetric strain nodal DOF and OSS stabilization
 * As a difference to its parent element, this one introduces the OSS projections as unknowns or, in other words, these are not linearized
 * This makes it possible to consider the subscale projections in the modal analysis as these are included into the eigensystem (otherwise they are not)
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallDisplacementMixedVolumetricStrainOssNonLinearElement
    : public SmallDisplacementMixedVolumetricStrainOssElement
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = SmallDisplacementMixedVolumetricStrainOssElement;

    /// The definition of the index type
    using IndexType = BaseType::IndexType;

    /// The definition of the sizetype
    using SizeType = BaseType::SizeType;

    /// This is the definition of the node
    using NodeType = typename BaseType::NodeType;

    // Counted pointer of SmallDisplacementMixedVolumetricStrainOssNonLinearElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementMixedVolumetricStrainOssNonLinearElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    SmallDisplacementMixedVolumetricStrainOssNonLinearElement()
    {};

    // Constructor using an array of nodes
    SmallDisplacementMixedVolumetricStrainOssNonLinearElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BaseType(
            NewId,
            pGeometry)
    {};

    // Constructor using an array of nodes with properties
    SmallDisplacementMixedVolumetricStrainOssNonLinearElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : BaseType(
            NewId,
            pGeometry,
            pProperties)
    {};

    // Copy constructor
    SmallDisplacementMixedVolumetricStrainOssNonLinearElement(SmallDisplacementMixedVolumetricStrainOssNonLinearElement const& rOther)
        : BaseType(rOther)
    {};

    // Destructor
    ~SmallDisplacementMixedVolumetricStrainOssNonLinearElement() override
    {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;

    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

    void GetSecondDerivativesVector(
        Vector &rValues,
        int Step = 0) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    const Parameters GetSpecifications() const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SmallDisplacementMixedVolumetricStrainOssNonLinearElement #" << Id();
        if (!BaseType::mConstitutiveLawVector.empty()
            && BaseType::mConstitutiveLawVector[0] != nullptr) {
          buffer << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        } else {
          buffer << " (no constitutive law)";
        }
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SmallDisplacementMixedVolumetricStrainOssNonLinearElement #" << Id();
        if (!BaseType::mConstitutiveLawVector.empty()
            && BaseType::mConstitutiveLawVector[0] != nullptr) {
          rOStream << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        } else {
          rOStream << " (no constitutive law)";
        }
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
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

    /**
     * @brief Calculates the finite element residual projection operator
     * @param rOrthogonalSubScalesOperator Reference to the residual projection operator to be filled
     * @param rProcessInfo Reference to the process info container
     */
    void CalculateOrthogonalSubScalesOperator(
        MatrixType &rOrthogonalSubScalesOperator,
        const ProcessInfo &rProcessInfo) const;

    /**
     * @brief Calculates the stabilization operator that applies the projections to the balance equations
     * @param rOrthogonalSubScalesStabilizationOperator Reference to the stabilization operator to be filled
     * @param rProcessInfo Reference to the process info container
     */
    void CalculateOrthogonalSubScalesStabilizationOperator(
        MatrixType &rOrthogonalSubScalesStabilizationOperator,
        const ProcessInfo &rProcessInfo) const;

    /**
     * @brief Calculates the residual lumped projection operator
     * @param rOrthogonalSubScalesLumpedProjectionOperator Reference to the lumped projection operator to be filled
     * @param rProcessInfo Reference to the process info container
     */
    void CalculateOrthogonalSubScalesLumpedProjectionOperator(
        MatrixType &rOrthogonalSubScalesLumpedProjectionOperator,
        const ProcessInfo &rProcessInfo) const;

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

    void save( Serializer& rSerializer ) const override;

    void load( Serializer& rSerializer ) override;

    ///@}
}; // class SmallDisplacementMixedVolumetricStrainOssNonLinearElement.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
