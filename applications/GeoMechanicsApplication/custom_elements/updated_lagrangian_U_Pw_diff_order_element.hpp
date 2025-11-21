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

#if !defined(KRATOS_GEO_UPDATED_LAGRANGIAN_U_PW_DIFF_ORDER_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_UPDATED_LAGRANGIAN_U_PW_DIFF_ORDER_ELEMENT_H_INCLUDED

// Project includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_utilities/element_utilities.hpp"
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
#endif // KRATOS_GEO_UPDATED_LAGRANGIAN_U_PW_DIFF_ORDER_ELEMENT_H_INCLUDED defined
