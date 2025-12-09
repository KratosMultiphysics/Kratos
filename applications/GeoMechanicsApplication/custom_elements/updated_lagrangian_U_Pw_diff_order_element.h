// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "stress_state_policy.h"

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
 * @class UpdatedLagrangianUPwDiffOrderElement
 * @brief Updated Lagrangian element for 2D and 3D geometries.
 * @details Implements an Updated Lagrangian definition for different order U-P elements. This works for arbitrary geometries in 2D and 3D
 * @author Vahid Galavi (Geomechanics)
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) UpdatedLagrangianUPwDiffOrderElement : public SmallStrainUPwDiffOrderElement
{
public:
    ///@name Type Definitions
    ///@{
    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType     = Vector;
    using MatrixType     = Matrix;

    /// Type definition for integration methods
    using IntegrationMethod = GeometryData::IntegrationMethod;

    /// The definition of the sizetype
    using SizeType = std::size_t;
    using SmallStrainUPwDiffOrderElement::CalculateDerivativesOnInitialConfiguration;
    using SmallStrainUPwDiffOrderElement::CalculateGreenLagrangeStrain;
    using SmallStrainUPwDiffOrderElement::mConstitutiveLawVector;
    using SmallStrainUPwDiffOrderElement::mStateVariablesFinalized;
    using SmallStrainUPwDiffOrderElement::mStressVector;

    using ElementVariables = typename SmallStrainUPwDiffOrderElement::ElementVariables;

    /// Counted pointer of UpdatedLagrangianUPwDiffOrderElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UpdatedLagrangianUPwDiffOrderElement);

    /// Default Constructor
    UpdatedLagrangianUPwDiffOrderElement() : SmallStrainUPwDiffOrderElement() {}

    /// Constructor using Geometry
    UpdatedLagrangianUPwDiffOrderElement(IndexType                          NewId,
                                         GeometryType::Pointer              pGeometry,
                                         std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                         std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : SmallStrainUPwDiffOrderElement(
              NewId, pGeometry, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
    }

    /// Constructor using Properties
    UpdatedLagrangianUPwDiffOrderElement(IndexType                          NewId,
                                         GeometryType::Pointer              pGeometry,
                                         PropertiesType::Pointer            pProperties,
                                         std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                         std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : SmallStrainUPwDiffOrderElement(
              NewId, pGeometry, pProperties, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
    }

    /// Destructor
    ~UpdatedLagrangianUPwDiffOrderElement() override = default;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param rNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType NewId, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    /**
     * @brief Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;
    using SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      bool               CalculateStiffnessMatrixFlag,
                      bool               CalculateResidualVectorFlag) override;

    std::vector<double> GetOptionalPermeabilityUpdateFactors(const std::vector<Vector>&) const override;

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
    void CalculateAndAddGeometricStiffnessMatrix(MatrixType&   rLeftHandSideMatrix,
                                                 const Vector& rStressVector,
                                                 const Matrix& rDNuDx,
                                                 const double  IntegrationCoefficient) const;

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    // Copy constructor
    UpdatedLagrangianUPwDiffOrderElement(UpdatedLagrangianUPwDiffOrderElement const& rOther);

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
}; // Class UpdatedLagrangianUPwDiffOrderElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
