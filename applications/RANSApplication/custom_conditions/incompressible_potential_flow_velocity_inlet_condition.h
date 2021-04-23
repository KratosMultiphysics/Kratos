//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_VELOCITY_INLET_CONDITION_H_INCLUDED)
#define KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_VELOCITY_INLET_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "includes/condition.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{
///@name Kratos Classes
///@{

/**
 * @brief
 *
 * @tparam TDim       Dimensionality of the condition (2D or 3D)
 * @tparam TNumNodes  Number of nodes in the condition
 */
template <
    unsigned int TDim,
    unsigned int TNumNodes = TDim>
class IncompressiblePotentialFlowVelocityInletCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IncompressiblePotentialFlowVelocityInletCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(IncompressiblePotentialFlowVelocityInletCondition);

    using BaseType = Condition;
    using NodeType = Node<3>;
    using PropertiesType = Properties;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    using IndexType = std::size_t;
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    explicit IncompressiblePotentialFlowVelocityInletCondition(
        IndexType NewId = 0)
    : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    IncompressiblePotentialFlowVelocityInletCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    IncompressiblePotentialFlowVelocityInletCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    IncompressiblePotentialFlowVelocityInletCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    IncompressiblePotentialFlowVelocityInletCondition(
        IncompressiblePotentialFlowVelocityInletCondition const& rOther)
    : BaseType(rOther)
    {
    }

    /// Destructor.
    ~IncompressiblePotentialFlowVelocityInletCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    IncompressiblePotentialFlowVelocityInletCondition& operator=(
        IncompressiblePotentialFlowVelocityInletCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new IncompressiblePotentialFlowVelocityInletCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<IncompressiblePotentialFlowVelocityInletCondition>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<IncompressiblePotentialFlowVelocityInletCondition>(
            NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Condition::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer p_new_condition = Create(
            NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

        p_new_condition->SetData(this->GetData());
        p_new_condition->SetFlags(this->GetFlags());

        return p_new_condition;
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType& ConditionDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    void GetValuesVector(
        VectorType& rValues,
        int Step = 0) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() override;

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

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class IncompressiblePotentialFlowVelocityInletCondition

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(
    std::istream& rIStream,
    IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const IncompressiblePotentialFlowVelocityInletCondition<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_VELOCITY_INLET_CONDITION_H_INCLUDED
