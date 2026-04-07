//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_STABILIZATION_PENALTY_CONDITION_3P_H_INCLUDED)
#define KRATOS_STABILIZATION_PENALTY_CONDITION_3P_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/condition.h"

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/iga_flags.h"

namespace Kratos
{
    /// 3P stabilization penalty condition: translation only
    class StabilizationPenaltyCondition3P
        : public Condition
    {
    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StabilizationPenaltyCondition3P);

        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        /// Constructor with Id and geometry
        StabilizationPenaltyCondition3P(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {}

        /// Constructor with Id, geometry and property
        StabilizationPenaltyCondition3P(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Condition(NewId, pGeometry, pProperties)
        {}

        /// Default constructor
        StabilizationPenaltyCondition3P() : Condition()
        {}

        /// Destructor
        ~StabilizationPenaltyCondition3P() override
        {}

        Condition::Pointer Create(
            IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive<StabilizationPenaltyCondition3P>(
                NewId, pGeom, pProperties);
        }

        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive<StabilizationPenaltyCondition3P>(
                NewId, GetGeometry().Create(ThisNodes), pProperties);
        }

        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetGeometry().size() * 3;

            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);

            noalias(rRightHandSideVector) = ZeroVector(mat_size);

            MatrixType left_hand_side_matrix;
            CalculateAll(
                left_hand_side_matrix,
                rRightHandSideVector,
                rCurrentProcessInfo,
                false,
                true
            );
        }

        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetGeometry().size() * 3;

            VectorType right_hand_side_vector;

            if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

            CalculateAll(
                rLeftHandSideMatrix,
                right_hand_side_vector,
                rCurrentProcessInfo,
                true,
                false
            );
        }

        void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetGeometry().size() * 3;

            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            noalias(rRightHandSideVector) = ZeroVector(mat_size);

            if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

            CalculateAll(
                rLeftHandSideMatrix,
                rRightHandSideVector,
                rCurrentProcessInfo,
                true,
                true
            );
        }

        void AddExplicitContribution(
            const VectorType& rRHS,
            const Variable<VectorType>& rRHSVariable,
            const Variable<array_1d<double,3> >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo
        ) override;

        void EquationIdVector(
            EquationIdVectorType& rResult,
            const ProcessInfo& rCurrentProcessInfo
        ) const override;

        void GetDofList(
            DofsVectorType& rElementalDofList,
            const ProcessInfo& rCurrentProcessInfo
        ) const override;

        void CalculateAll(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo,
            const bool CalculateStiffnessMatrixFlag,
            const bool CalculateResidualVectorFlag
        );

        int Check(const ProcessInfo& rCurrentProcessInfo) const override;

        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "\"StabilizationPenaltyCondition3P\" #" << Id();
            return buffer.str();
        }

        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"StabilizationPenaltyCondition3P\" #" << Id();
        }

        void PrintData(std::ostream& rOStream) const override
        {
            pGetGeometry()->PrintData(rOStream);
        }

    private:
        friend class Serializer;

        void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        }
    };

} // namespace Kratos

#endif // KRATOS_STABILIZATION_PENALTY_CONDITION_3P_H_INCLUDED