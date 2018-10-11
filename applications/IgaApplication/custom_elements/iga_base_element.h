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

template <std::size_t TDofsPerNode>
class IgaBaseElement
    : public Element
{
public:
    using IgaBaseElementType = IgaBaseElement<TDofsPerNode>;

    using Vector3 = BoundedVector<double, 3>;

    using Element::Element;

    static constexpr inline std::size_t DofsPerNode()
    {
        return TDofsPerNode;
    }

    std::size_t inline NumberOfNodes() const
    {
        return GetGeometry().size();
    }

    std::size_t inline NumberOfDofs() const
    {
        return NumberOfNodes() * DofsPerNode();
    }

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

    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide) = 0;

    std::string Info() const override
    {
        std::stringstream buffer;
        PrintInfo(buffer);
        return buffer.str();
    }

    void PrintData(std::ostream& rOStream) const
    {
        pGetGeometry()->PrintData(rOStream);
    }

protected:

    template <typename TVariable>
    void inline SetElementDof(
        DofsVectorType& rElementalDofList,
        const std::size_t NodeIndex,
        const std::size_t DofIndex,
        const TVariable& variable)
    {
        Node<3>& node = GetGeometry()[NodeIndex];

        rElementalDofList[NodeIndex * DofsPerNode() + DofIndex] =
            node.pGetDof(variable);
    }

    template <typename TVariable>
    void inline SetEquationId(
        EquationIdVectorType& rResult,
        const std::size_t NodeIndex,
        const std::size_t DofIndex,
        const TVariable& variable)
    {
        Node<3>& node = GetGeometry()[NodeIndex];

        rResult[NodeIndex * DofsPerNode() + DofIndex] =
            node.GetDof(variable).EquationId();
    }

    static inline std::size_t GetDofType(
        std::size_t DofIndex)
    {
        return DofIndex % DofsPerNode();
    }

    static inline std::size_t GetShapeIndex(
        std::size_t DofIndex)
    {
        return DofIndex / DofsPerNode();
    }
};

} // namespace Kratos

#endif // !defined(KRATOS_IGA_BASE_ELEMENT_H_INCLUDED)
