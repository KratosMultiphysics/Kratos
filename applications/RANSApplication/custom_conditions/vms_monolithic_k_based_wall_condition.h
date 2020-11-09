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

#if !defined(KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_WALL_H_INCLUDED)
#define KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_WALL_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "custom_conditions/monolithic_wall_condition.h"

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
template <unsigned int TDim, unsigned int TNumNodes = TDim>
class VMSMonolithicKBasedWallCondition : public MonolithicWallCondition<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VMSMonolithicKBasedWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(VMSMonolithicKBasedWallCondition);

    using BaseType = MonolithicWallCondition<TDim, TNumNodes>;
    using NodeType = Node<3>;
    using PropertiesType = Properties;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    explicit VMSMonolithicKBasedWallCondition(
        IndexType NewId = 0)
    : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    VMSMonolithicKBasedWallCondition(
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
    VMSMonolithicKBasedWallCondition(
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
    VMSMonolithicKBasedWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    VMSMonolithicKBasedWallCondition(
        VMSMonolithicKBasedWallCondition const& rOther)
    : BaseType(rOther)
    {
    }

    /// Destructor.
    ~VMSMonolithicKBasedWallCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    VMSMonolithicKBasedWallCondition& operator=(
        VMSMonolithicKBasedWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new VMSMonolithicKBasedWallCondition object.
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
        return Kratos::make_intrusive<VMSMonolithicKBasedWallCondition>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<VMSMonolithicKBasedWallCondition>(
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

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

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

protected:
    ///@name Protected Operations
    ///@{

    /**
     * @brief Applies rans based wall law
     *
     * This method calculate left hand side matrix and right hand side vector for following equation
     *
     * \[
     *      u_\tau = max\left(C_\mu^0.25 \sqrt{k}, \frac{||\underline{u}||}{\frac{1}{\kappa}ln(y^+)+\beta}\right)
     * \]
     *
     * integration point value = \rho \frac{u_\tau^2}{||\underline{u}||}\underline{u}
     *
     * @param rLocalMatrix         Left hand side matrix
     * @param rLocalVector         Right hand side vector
     * @param rCurrentProcessInfo  Current process info from model part
     */
    void ApplyWallLaw(
        MatrixType& rLocalMatrix,
        VectorType& rLocalVector,
        ProcessInfo& rCurrentProcessInfo) override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    double mWallHeight;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class VMSMonolithicKBasedWallCondition

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(
    std::istream& rIStream,
    VMSMonolithicKBasedWallCondition<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const VMSMonolithicKBasedWallCondition<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_WALL_H_INCLUDED
