//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_SUPPORT_PENALTY_CONDITION_H_INCLUDED )
#define  KRATOS_SUPPORT_PENALTY_CONDITION_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"

namespace Kratos
{
    /// Condition for penalty support condition
    class SupportPenaltyCondition
        : public Condition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer definition of SupportPenaltyCondition
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SupportPenaltyCondition);

        /// Size types
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor with Id and geometry
        SupportPenaltyCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {};

        /// Constructor with Id, geometry and property
        SupportPenaltyCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Condition(NewId, pGeometry, pProperties)
        {};

        /// Default constructor
        SupportPenaltyCondition() : Condition()
        {};

        /// Destructor
        virtual ~SupportPenaltyCondition() override
        {};

        ///@}
        ///@name Life Cycle
        ///@{

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive<SupportPenaltyCondition>(
                NewId, pGeom, pProperties);
        };

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive<SupportPenaltyCondition>(
                NewId, GetGeometry().Create(ThisNodes), pProperties);
        };

        ///@}
        ///@name Operations
        ///@{

        /**
        * @brief This is called during the assembling process in order
        *        to calculate the condition right hand side matrix
        * @param rLeftHandSideMatrix the condition right hand side matrix
        * @param rCurrentProcessInfo the current process info
        */
        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetGeometry().size() * 3;

            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            noalias(rRightHandSideVector) = ZeroVector(mat_size);

            MatrixType left_hand_side_matrix;

            CalculateAll(left_hand_side_matrix, rRightHandSideVector,
                rCurrentProcessInfo, false, true);
        }

        /**
        * @brief This is called during the assembling process in order
        *        to calculate the condition left hand side matrix
        * @param rLeftHandSideMatrix the condition left hand side matrix
        * @param rCurrentProcessInfo the current process info
        */
        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetGeometry().size() * 3;

            VectorType right_hand_side_vector;

            if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2())
                rLeftHandSideMatrix.resize(mat_size, mat_size);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

            CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
                rCurrentProcessInfo, true, false);
        }

        /**
         * @brief This function provides a more general interface to the element.
         * @details It is designed so that rLHSvariables and rRHSvariables are
         *          passed to the element thus telling what is the desired output
         * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
         * @param rRightHandSideVector container for the desired RHS output
         * @param rCurrentProcessInfo the current process info instance
         */
        void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetGeometry().size() * 3;

            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            noalias(rRightHandSideVector) = ZeroVector(mat_size);

            if (rLeftHandSideMatrix.size1() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

            CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                rCurrentProcessInfo, true, true);
        }

        /**
        * @brief Sets on rResult the ID's of the element degrees of freedom
        * @param rResult The vector containing the equation id
        * @param rCurrentProcessInfo The current process info instance
        */
        void EquationIdVector(
            EquationIdVectorType& rResult,
            const ProcessInfo& rCurrentProcessInfo
        ) const override;

        /**
        * @brief Sets on rConditionDofList the degrees of freedom of the considered element geometry
        * @param rElementalDofList The vector containing the dof of the element
        * @param rCurrentProcessInfo The current process info instance
        */
        void GetDofList(
            DofsVectorType& rElementalDofList,
            const ProcessInfo& rCurrentProcessInfo
        ) const override;

        /// Calculates left (K) and right (u) hand sides, according to the flags
        void CalculateAll(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo,
            const bool CalculateStiffnessMatrixFlag,
            const bool CalculateResidualVectorFlag
        );

        void DeterminantOfJacobianInitial(
            const GeometryType& rGeometry,
            Vector& rDeterminantOfJacobian);

        ///@}
        ///@name Check
        ///@{

        /// Performs check if Penalty factor is provided.
        int Check(const ProcessInfo& rCurrentProcessInfo) const override;

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "\"SupportPenaltyCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SupportPenaltyCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override
        {
            pGetGeometry()->PrintData(rOStream);
        }

        ///@}

    private:
        ///@name Serialization
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        }

        ///@}

    }; // Class SupportPenaltyCondition

}  // namespace Kratos.

#endif // KRATOS_SUPPORT_PENALTY_CONDITION_H_INCLUDED  defined
