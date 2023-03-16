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

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

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
 * @class SmallDisplacementMixedVolumetricStrainModalAnalysisElement
 * @ingroup StructuralMechanicsApplication
 * @brief Small displacement with strain based mixed formulation element
 * @details This implements a small displacements element formulation with an extra volumetric strain nodal DOF
 * Note that this element is specifically designed to include the Orthogonal SubScales projections into the modal analysis
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallDisplacementMixedVolumetricStrainModalAnalysisElement
    : public SmallDisplacementMixedVolumetricStrainElement
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = SmallDisplacementMixedVolumetricStrainElement;

    ///Reference type definition for constitutive laws
    using ConstitutiveLawType = typename BaseType::ConstitutiveLawType;

    ///Pointer type for constitutive laws
    using ConstitutiveLawPointerType = typename BaseType::ConstitutiveLawPointerType;

    ///StressMeasure from constitutive laws
    using StressMeasureType = typename BaseType::StressMeasureType;

    ///Type definition for integration methods
    using IntegrationMethod = typename BaseType::IntegrationMethod;

    /// This is the definition of the node
    using NodeType = typename BaseType::NodeType;

    // Counted pointer of SmallDisplacementMixedVolumetricStrainModalAnalysisElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementMixedVolumetricStrainModalAnalysisElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    SmallDisplacementMixedVolumetricStrainModalAnalysisElement()
    {};

    // Constructor using an array of nodes
    SmallDisplacementMixedVolumetricStrainModalAnalysisElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BaseType(
            NewId,
            pGeometry)
    {};

    // Constructor using an array of nodes with properties
    SmallDisplacementMixedVolumetricStrainModalAnalysisElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : BaseType(
            NewId,
            pGeometry,
            pProperties)
    {};

    // Copy constructor
    SmallDisplacementMixedVolumetricStrainModalAnalysisElement(SmallDisplacementMixedVolumetricStrainModalAnalysisElement const& rOther)
        : BaseType(rOther)
    {};

    // Destructor
    ~SmallDisplacementMixedVolumetricStrainModalAnalysisElement() override
    {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // /**
    //  * @brief Called to initialize the element.
    //  * @warning Must be called before any calculation is done
    //  */
    // void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    // /**
    //  * @brief Called at the beginning of each solution step
    //  * @param rCurrentProcessInfo the current process info instance
    //  */
    // void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    // /**
    //  * @brief Called at the end of eahc solution step
    //  * @param rCurrentProcessInfo the current process info instance
    //  */
    // void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

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

    // void Calculate(
    //     const Variable<double>& rVariable,
    //     double& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // void Calculate(
    //     const Variable<array_1d<double, 3>>& rVariable,
    //     array_1d<double, 3>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    // /**
    //  * @brief Calculate a double Variable on the element Gauss points
    //  * @param rVariable The variable we want to get
    //  * @param rOutput The values obtained int the integration points
    //  * @param rCurrentProcessInfo the current process info instance
    //  */
    // void CalculateOnIntegrationPoints(
    //     const Variable<double>& rVariable,
    //     std::vector<double>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // /**
    //  * @brief Calculate a Vector Variable on the element Gauss points
    //  * @param rVariable The variable we want to get
    //  * @param rOutput The values obtained int the integration points
    //  * @param rCurrentProcessInfo the current process info instance
    //  */
    // void CalculateOnIntegrationPoints(
    //     const Variable<Vector>& rVariable,
    //     std::vector<Vector>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // /**
    //  * @brief Calculate an array Variable on the element Gauss points
    //  * @param rVariable The variable we want to get
    //  * @param rOutput The values obtained int the integration points
    //  * @param rCurrentProcessInfo the current process info instance
    //  */
    // void CalculateOnIntegrationPoints(
    //     const Variable<array_1d<double,3>>& rVariable,
    //     std::vector<array_1d<double,3>>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

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
        buffer << "SmallDisplacementMixedVolumetricStrainModalAnalysisElement #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SmallDisplacementMixedVolumetricStrainModalAnalysisElement #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
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

}; // class SmallDisplacementMixedVolumetricStrainModalAnalysisElement.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
