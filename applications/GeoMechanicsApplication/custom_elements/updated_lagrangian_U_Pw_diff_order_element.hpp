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
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/math_utilities.h"
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
template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UpdatedLagrangianUPwDiffOrderElement
    : public SmallStrainUPwDiffOrderElement<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
    using BaseType       = SmallStrainUPwDiffOrderElement<TDim, TNumNodes>;
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
    using BaseType::CalculateDerivativesOnInitialConfiguration;
    using BaseType::CalculateGreenLagrangeStrain;
    using BaseType::mConstitutiveLawVector;
    using BaseType::mStateVariablesFinalized;
    using BaseType::mStressVector;

    using ElementVariables = typename BaseType::ElementVariables;

    /// Counted pointer of UpdatedLagrangianUPwDiffOrderElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UpdatedLagrangianUPwDiffOrderElement);

    /// Default Constructor
    UpdatedLagrangianUPwDiffOrderElement() : SmallStrainUPwDiffOrderElement<TDim, TNumNodes>() {}

    /// Constructor using Geometry
    UpdatedLagrangianUPwDiffOrderElement(IndexType                          NewId,
                                         GeometryType::Pointer              pGeometry,
                                         std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                         std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : SmallStrainUPwDiffOrderElement<TDim, TNumNodes>(
              NewId, pGeometry, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
    }

    /// Constructor using Properties
    UpdatedLagrangianUPwDiffOrderElement(IndexType                          NewId,
                                         GeometryType::Pointer              pGeometry,
                                         PropertiesType::Pointer            pProperties,
                                         std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                         std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : SmallStrainUPwDiffOrderElement<TDim, TNumNodes>(
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
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<UpdatedLagrangianUPwDiffOrderElement>(
            NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
            this->CloneIntegrationCoefficientModifier());
    }

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param rNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType NewId, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const override
    {
        return Create(NewId, this->GetGeometry().Create(rNodes), pProperties);
    }

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        if (rVariable == REFERENCE_DEFORMATION_GRADIENT_DETERMINANT) {
            rOutput = GeoMechanicsMathUtilities::CalculateDeterminants(this->CalculateDeformationGradients());
        } else {
            BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    /**
     * @brief Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        rOutput.resize(mConstitutiveLawVector.size());

        if (rVariable == REFERENCE_DEFORMATION_GRADIENT) {
            rOutput = this->CalculateDeformationGradients();
        } else {
            BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        KRATOS_TRY

        rOutput.resize(this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()));

        if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
            const auto deformation_gradients = this->CalculateDeformationGradients();

            std::ranges::transform(deformation_gradients, rOutput.begin(), [this](const Matrix& rDeformationGradient) {
                return this->CalculateGreenLagrangeStrain(rDeformationGradient);
            });
        } else {
            BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

    using BaseType::CalculateOnIntegrationPoints;

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
    std::string Info() const override
    {
        const std::string constitutive_info =
            !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";
        return "Updated Lagrangian U-Pw different order Element #" + std::to_string(this->Id()) +
               "\nConstitutive law: " + constitutive_info;
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->pGetGeometry()->PrintData(rOStream);
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
                      bool               CalculateResidualVectorFlag) override
    {
        KRATOS_TRY

        BaseType::CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                               CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

        ElementVariables variables;
        this->InitializeElementVariables(variables, rCurrentProcessInfo);

        if (CalculateStiffnessMatrixFlag && variables.ConsiderGeometricStiffness) {
            const auto& r_integration_points =
                this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
            const auto integration_coefficients =
                this->CalculateIntegrationCoefficients(r_integration_points, variables.detJuContainer);

            for (IndexType g_point = 0; g_point < r_integration_points.size(); ++g_point) {
                this->CalculateAndAddGeometricStiffnessMatrix(rLeftHandSideMatrix, mStressVector[g_point],
                                                              variables.DNu_DXContainer[g_point],
                                                              integration_coefficients[g_point]);
            }
        }

        KRATOS_CATCH("")
    }

    std::vector<double> GetOptionalPermeabilityUpdateFactors(const std::vector<Vector>&) const override
    {
        return {};
    }

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
                                                 const double  IntegrationCoefficient) const
    {
        KRATOS_TRY

        const Matrix reduced_Kg_matrix =
            prod(rDNuDx, Matrix(prod(MathUtils<double>::StressVectorToTensor(rStressVector), trans(rDNuDx)))) *
            IntegrationCoefficient;

        Matrix geometric_stiffness_matrix = ZeroMatrix(TNumNodes * TDim, TNumNodes * TDim);
        MathUtils<double>::ExpandAndAddReducedMatrix(geometric_stiffness_matrix, reduced_Kg_matrix, TDim);

        GeoElementUtilities::AssembleUUBlockMatrix(rLeftHandSideMatrix, geometric_stiffness_matrix);

        KRATOS_CATCH("")
    }

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    }
}; // Class UpdatedLagrangianUPwDiffOrderElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
