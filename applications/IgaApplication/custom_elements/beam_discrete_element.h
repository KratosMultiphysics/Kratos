/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if !defined(KRATOS_BEAM_DISCRETE_ELEMENT_H_INCLUDED)
#define KRATOS_BEAM_DISCRETE_ELEMENT_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

// External includes

// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"


namespace Kratos
{

template <int TDofsPerNode>
class IgaBaseElement
    : public Element
{
public:
    using IgaBaseElementType = IgaBaseElement<TDofsPerNode>;

    using Vector3D = BoundedVector<double, 3>;

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

    template <typename TVariable>
    void inline SetDof(
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
};

class BeamDiscreteElement
    : public IgaBaseElement<3>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( BeamDiscreteElement );

    using IgaBaseElementType::IgaBaseElementType;

    ~BeamDiscreteElement() override
    {
    };

    Vector3D mReferenceBaseVector;

    Element::Pointer BeamDiscreteElement::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        auto geometry = GetGeometry().Create(ThisNodes);

        return Kratos::make_shared<BeamDiscreteElement>(NewId, geometry,
            pProperties);
    }

    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rElementalDofList.resize(NumberOfDofs());

        for (std::size_t i = 0; i < NumberOfNodes(); i++) {
            SetDof(rElementalDofList, i, 0, DISPLACEMENT_X);
            SetDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
            SetDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
        }

        KRATOS_CATCH("")
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rResult.resize(NumberOfDofs());

        for (std::size_t i = 0; i < NumberOfNodes(); i++) {
            SetEquationId(rResult, i, 0, DISPLACEMENT_X);
            SetEquationId(rResult, i, 1, DISPLACEMENT_Y);
            SetEquationId(rResult, i, 2, DISPLACEMENT_Z);
        }

        KRATOS_CATCH("")
    }

    void Initialize() override
    {
        mReferenceBaseVector = GetActualBaseVector();
    }

    Vector3D GetActualBaseVector() const
    {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        Vector3D actual_base_vector = ZeroVector(3);

        for (std::size_t i = 0; i < NumberOfNodes(); i++)
        {
            actual_base_vector[0] += DN_De(0, i) * GetGeometry()[i].X();
            actual_base_vector[1] += DN_De(0, i) * GetGeometry()[i].Y();
            actual_base_vector[2] += DN_De(0, i) * GetGeometry()[i].Z();
        }

        return actual_base_vector;
    }

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide) override
    {
        KRATOS_TRY;

        // get integration data
        
        const double& integration_weight = GetValue(INTEGRATION_WEIGHT);
        Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        // get properties

        const auto& properties = GetProperties();

        const double E = properties[YOUNG_MODULUS];
        const double A = properties[CROSS_AREA];
        const double prestress = properties[PRESTRESS_CAUCHY];

        // compute base vectors

        const Vector3D actual_base_vector = GetActualBaseVector();

        const double reference_a = norm_2(mReferenceBaseVector);
        const double actual_a = norm_2(actual_base_vector);

        const double actual_aa = actual_a * actual_a;
        const double reference_aa = reference_a * reference_a;

        // green-lagrange strain

        const double e11_membrane = 0.5 * (actual_aa - reference_aa);

        // normal force

        const double s11_membrane = prestress * A + e11_membrane * A * E /
            reference_aa;

        for (std::size_t r = 0; r < NumberOfDofs(); r++) {
            const std::size_t dof_type_r = r % DofsPerNode();
            const std::size_t shape_index_r = r / DofsPerNode();

            const double epsilon_var_r = actual_base_vector[dof_type_r] *
                shape_derivatives(shape_index_r, 0) / reference_aa;

            if (ComputeLeftHandSide) {
                for (std::size_t s = 0; s < NumberOfDofs(); s++) {
                    const std::size_t dof_type_s = s % DofsPerNode();
                    const std::size_t shape_index_s = s / DofsPerNode();

                    const double epsilon_var_s =
                        actual_base_vector[dof_type_s] *
                        shape_derivatives(shape_index_s, 0) / reference_aa;

                    rLeftHandSideMatrix(r, s) = E * A * epsilon_var_r *
                        epsilon_var_s;

                    if (dof_type_r == dof_type_s) {
                        const double epsilon_var_rs =
                            shape_derivatives(shape_index_r, 0) *
                            shape_derivatives(shape_index_s, 0) / reference_aa;

                        rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs;
                    }
                }
            }

            if (ComputeRightHandSide) {
                rRightHandSideVector[r] = -s11_membrane * epsilon_var_r;
            }
        }

        if (ComputeLeftHandSide) {
            rLeftHandSideMatrix *= reference_a * integration_weight;
        }

        if (ComputeRightHandSide) {
            rRightHandSideVector *= reference_a * integration_weight;
        }

        KRATOS_CATCH("")
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"BeamDiscreteElement\" #" << Id();
    }
};

} // namespace Kratos

#endif // !defined(KRATOS_BEAM_DISCRETE_ELEMENT_H_INCLUDED)
