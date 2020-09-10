/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if !defined(KRATOS_IGA_BASE_ELEMENT_H_INCLUDED)
#define KRATOS_IGA_BASE_ELEMENT_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/element.h"

// External includes

// Project includes

namespace Kratos
{

/** Base element for isogeometric elements.
 */
template <std::size_t TDofsPerNode>
class IgaBaseElement
    : public Element
{
public:
    using IgaBaseElementType = IgaBaseElement<TDofsPerNode>;

    // Alias for a three-dimensional vector - used a lot.
    using Vector3 = BoundedVector<double, 3>;

    // Inherit constructors of the Kratos element.
    using Element::Element;

    /** Gets the number of dofs per node.
     *
     * @return Number of dofs per node.
     */
    static constexpr inline std::size_t DofsPerNode()
    {
        return TDofsPerNode;
    }

    /** Gets the number of nodes.
     *
     * @return Number of nodes.
     */
    std::size_t inline NumberOfNodes() const
    {
        return GetGeometry().size();
    }

    /** Gets the number of degrees of freedom.
     *
     * @return Number of degrees of freedom.
     */
    std::size_t inline NumberOfDofs() const
    {
        return NumberOfNodes() * DofsPerNode();
    }

    /** Calculates the elemental left- and right-hand-side
     *
     * @param rLeftHandSideMatrix Elemental left-hand-side matrix
     * @param rRightHandSideVector Elemental right-hand-side
     * @param rCurrentProcessInfo Current process info
     *
     * @note Child-classes should implement CalculateAll.
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        if (rLeftHandSideMatrix.size1() != number_of_dofs
            || rLeftHandSideMatrix.size2() != number_of_dofs) {
            rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
        }

        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(number_of_dofs);
        }

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    /** Calculates the elemental left-hand-side
     *
     * @param rLeftHandSideMatrix Elemental left-hand-side matrix
     * @param rCurrentProcessInfo Current process info
     *
     * @note Child-classes should implement CalculateAll.
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        if (rLeftHandSideMatrix.size1() != number_of_dofs
            || rLeftHandSideMatrix.size2() != number_of_dofs) {
            rLeftHandSideMatrix.resize(number_of_dofs, number_of_dofs);
        }

        VectorType right_hand_side_vector = Vector(0);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    /** Calculates the elemental right-hand-side
     *
     * @param rRightHandSideVector Elemental right-hand-side vector
     * @param rCurrentProcessInfo Current process info
     *
     * @note Child-classes should implement CalculateAll.
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override
    {
        const std::size_t number_of_dofs = NumberOfDofs();

        MatrixType left_hand_side_matrix = Matrix(0, 0);

        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(number_of_dofs);
        }

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    /** Calculates the elemental left- and right-hand-side
     *
     * @note This function should be implemented by the child-classes to
     *       calculate left- and right-hand-side
     *
     * @param rLeftHandSideMatrix Elemental left-hand-side matrix.
     * @param rRightHandSideVector Elemental right-hand-side vector.
     * @param rCurrentProcessInfo Current process info.
     * @param ComputeLeftHandSide True whether the left-hand-side matrix
     *                            should be calculated.
     * @param ComputeRightHandSide True whether the right-hand-side vector
     *                             should be calculated.
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide) = 0;

    /** Get the geometry information as a string.
     *
     * @return The geometry information as a string.
     */
    std::string Info() const override
    {
        std::stringstream buffer;
        PrintInfo(buffer);
        return buffer.str();
    }

    /** Write the geometry info to a stream.
     *
     * @param rOStream Output stream.
     */
    void PrintData(
        std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

protected:

    /** Helper method for setting-up the elemental list of degrees of freedom.
     *
     * @param rElementalDofList Elemental list of degrees of freedom.
     * @param NodeIndex Index of the node.
     * @param DofTypeIndex Index of the degree of freedom type.
     * @param Variable Variable of the degree of freedom.
     */
    template <typename TVariable>
    void inline SetElementDof(
        DofsVectorType& rElementalDofList,
        const std::size_t NodeIndex,
        const std::size_t DofTypeIndex,
        const TVariable& Variable)
    {
        Node<3>& node = GetGeometry()[NodeIndex];

        rElementalDofList[NodeIndex * DofsPerNode() + DofTypeIndex] =
            node.pGetDof(Variable);
    }

    /** Helper method for setting-up the elemental list of equation ids.
     *
     * @param rResult Elemental list of equation ids.
     * @param NodeIndex Index of the node.
     * @param DofTypeIndex Index of the degree of freedom type.
     * @param Variable Variable of the degree of freedom.
     */
    template <typename TVariable>
    void inline SetElementEquationId(
        EquationIdVectorType& rResult,
        const std::size_t NodeIndex,
        const std::size_t DofTypeIndex,
        const TVariable& Variable)
    {
        Node<3>& node = GetGeometry()[NodeIndex];

        rResult[NodeIndex * DofsPerNode() + DofTypeIndex] =
            node.GetDof(Variable).EquationId();
    }

    /** Helper method for getting the index of the degree of freedom type.
     *
     * @param DofIndex Index of the degree of freedom.
     *
     * @return The index of the degree of freedom type.
     */
    static inline std::size_t GetDofTypeIndex(
        std::size_t DofIndex)
    {
        return DofIndex % DofsPerNode();
    }

    /** Helper method for getting the index of the shape function.
     *
     * @param DofIndex Index of the degree of freedom.
     *
     * @return The index of the shape function.
     */
    static inline std::size_t GetShapeIndex(
        std::size_t DofIndex)
    {
        return DofIndex / DofsPerNode();
    }
};

} // namespace Kratos

#endif // !defined(KRATOS_IGA_BASE_ELEMENT_H_INCLUDED)
