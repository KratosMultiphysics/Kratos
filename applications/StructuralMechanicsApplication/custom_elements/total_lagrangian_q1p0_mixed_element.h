// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//


#if !defined(KRATOS_TOTAL_LAGRANGIAN_Q1P0_MIXED_ELEMENT_H_INCLUDED )
#define  KRATOS_TOTAL_LAGRANGIAN_Q1P0_MIXED_ELEMENT_H_INCLUDED


// System includes


// External include

// Project includes
#include "includes/define.h"
#include "custom_elements/total_lagrangian.h"
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
 * @class TotalLagrangianQ1P0MixedElement
 * @ingroup StructuralMechanicsApplication
 * @brief Total Lagrangian mixed u-p element (Q1P0) for 2D and 3D geometries.
 * @details Implements a mixed u-p total Lagrangian element for structural analysis, especially when incompressibility is used at constitutive level. This works for arbitrary geometries in 2D and 3D
 * @details Reference: Nonlinear Finite Element Methods, prof. Peter Wriggers, Springer, section 10.2.1. "Mixed Q1-P0 Element"
 * @author Alejandro Cornejo
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TotalLagrangianQ1P0MixedElement
    : public TotalLagrangian
{
public:
    ///@name Type Definitions
    ///@{

    /// The base element type
    typedef TotalLagrangian BaseType;

    /// Counted pointer of TotalLagrangianQ1P0MixedElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TotalLagrangianQ1P0MixedElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TotalLagrangianQ1P0MixedElement(IndexType NewId, GeometryType::Pointer pGeometry);
    TotalLagrangianQ1P0MixedElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    TotalLagrangianQ1P0MixedElement(TotalLagrangianQ1P0MixedElement const& rOther)
        :BaseType(rOther)
    {};

    /// Destructor.
    ~TotalLagrangianQ1P0MixedElement() override;

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
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

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
        std::stringstream buffer;
        buffer << "TotalLagrangianQ1P0MixedElement #" << Id() << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TotalLagrangianQ1P0MixedElement #" << Id() << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
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

    TotalLagrangianQ1P0MixedElement() : TotalLagrangian()
    {
    }

    /**
     * @brief This functions calculates both the RHS and the LHS
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
     * @brief It computes the stress vector and tangent constitutive tensor by
     * automatic differentiation
     */
    void CalculateNeoHookeanStressAndTangent(const Matrix &rC, const double Pressure, const double LameMu, Vector &rStress, Matrix &rTangentTensor);

    /**
     * @brief It computes the current volume of the element
     */
    // double GetCurrentVolume() const;

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
    double mPressure = 0.0;

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

    ///@}
    ///@name Serialization
    ///@{
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
    //TotalLagrangian& operator=(const TotalLagrangian& rOther);
    /// Copy constructor.
    //TotalLagrangian(const TotalLagrangian& rOther);
    ///@}

}; // Class TotalLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_LAGRANGIAN_Q1P0_MIXED_ELEMENT_H_INCLUDED  defined
