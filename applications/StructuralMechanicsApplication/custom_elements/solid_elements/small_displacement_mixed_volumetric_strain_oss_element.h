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
#include "small_displacement_mixed_volumetric_strain_element.h"

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
 * @class SmallDisplacementMixedVolumetricStrainOssElement
 * @ingroup StructuralMechanicsApplication
 * @brief Small displacement with strain based mixed formulation element
 * @details This implements a small displacements element formulation with an extra volumetric strain nodal DOF and OSS stabilization
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallDisplacementMixedVolumetricStrainOssElement
    : public SmallDisplacementMixedVolumetricStrainElement
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = SmallDisplacementMixedVolumetricStrainElement;

    /// The definition of the index type
    using IndexType = BaseType::IndexType;

    /// The definition of the sizetype
    using SizeType = BaseType::SizeType;

    // Counted pointer of SmallDisplacementMixedVolumetricStrainOssElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SmallDisplacementMixedVolumetricStrainOssElement );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    SmallDisplacementMixedVolumetricStrainOssElement()
    {};

    // Constructor using an array of nodes
    SmallDisplacementMixedVolumetricStrainOssElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BaseType(
            NewId,
            pGeometry)
    {};

    // Constructor using an array of nodes with properties
    SmallDisplacementMixedVolumetricStrainOssElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : BaseType(
            NewId,
            pGeometry,
            pProperties)
    {};

    // Copy constructor
    SmallDisplacementMixedVolumetricStrainOssElement(SmallDisplacementMixedVolumetricStrainOssElement const &rOther)
        : BaseType(rOther)
    {};

    // Destructor
    ~SmallDisplacementMixedVolumetricStrainOssElement() override
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

    void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

    void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double,3>>& rVariable,
        std::vector<array_1d<double,3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

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
        buffer << "SmallDisplacementMixedVolumetricStrainOssElement #" << Id();
        if (!BaseType::mConstitutiveLawVector.empty()) {
          buffer << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        } else {
          buffer << "\nNo constitutive law.";
        }
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SmallDisplacementMixedVolumetricStrainOssElement #" << Id();
        if (!BaseType::mConstitutiveLawVector.empty()) {
          rOStream << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        } else {
          rOStream << "\nNo constitutive law.";
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
     * @brief Calculates the operator to be applied to the OSS projections
     * @param rOrthogonalSubScalesOperator Matrix to store the OSS projection operator
     * @param rThisKinematicVariables Reference to the kinematics container
     * @param rThisGaussPointAuxiliaryVariables Reference to the Gauss point data container
     */
    void CalculateOssStabilizationOperatorGaussPointContribution(
        MatrixType &rOrthogonalSubScalesOperator,
        const KinematicVariables &rThisKinematicVariables,
        const GaussPointAuxiliaryVariables& rThisGaussPointAuxiliaryVariables) const;

    void UpdateGaussPointDisplacementSubscaleHistory(
        const KinematicVariables& rThisKinematicVariables,
        const ConstitutiveVariables& rThisConstitutiveVariables,
        const ProcessInfo& rProcessInfo,
        const IndexType PointIndex) override;

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

    void save( Serializer& rSerializer ) const override;

    void load( Serializer& rSerializer ) override;

    ///@}
}; // class SmallDisplacementMixedVolumetricStrainOssElement.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
