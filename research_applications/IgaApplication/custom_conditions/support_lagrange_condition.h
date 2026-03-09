//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_SUPPORT_LAGRANGE_CONDITION_H_INCLUDED )
#define  KRATOS_SUPPORT_LAGRANGE_CONDITION_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"


namespace Kratos
{

    /// Lagrange Multiplier based support condition.
    /** This condition can be used to apply support conditions
    *   to all types of geometry discretizations with the
    *   Lagrange Multiplier approach.
    *
    *   The aproach is described in https://doi.org/10.1186/s40323-018-0109-4
    */
    class KRATOS_API(IGA_APPLICATION) SupportLagrangeCondition
        : public Condition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of SupportLagrangeCondition
        KRATOS_CLASS_POINTER_DEFINITION(SupportLagrangeCondition);

        /// Size types
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor with Id and geometry
        SupportLagrangeCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {};

        /// Constructor with Id, geometry and property
        SupportLagrangeCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Condition(NewId, pGeometry, pProperties)
        {};

        /// Default constructor
        SupportLagrangeCondition()
            : Condition()
        {};

        /// Destructor.
        virtual ~SupportLagrangeCondition() = default;

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
            return Kratos::make_intrusive<SupportLagrangeCondition>(
                NewId, pGeom, pProperties);
        };

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties
            ) const override
        {
            return Kratos::make_intrusive< SupportLagrangeCondition >(
                NewId, GetGeometry().Create(ThisNodes), pProperties);
        };

        ///@}
        ///@name Operations
        ///@{

        /// Calculates RightHandSide F-vector only
        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetNumberOfNonZeroNodes() * 6;

            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            noalias(rRightHandSideVector) = ZeroVector(mat_size);

            MatrixType left_hand_side_matrix;

            CalculateAll(left_hand_side_matrix, rRightHandSideVector,
                rCurrentProcessInfo, false, true);
        }

        /// Calculates LeftHandSide K-Matrix only
        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetNumberOfNonZeroNodes() * 6;

            VectorType right_hand_side_vector;

            if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2())
                rLeftHandSideMatrix.resize(mat_size, mat_size);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

            CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
                rCurrentProcessInfo, true, false);
        }

        /// Calculates LeftHandSide K-Matrix and f-vector
        void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            const SizeType mat_size = GetNumberOfNonZeroNodes() * 6;

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
        * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
        * @param rElementalDofList The vector containing the dof of the element
        * @param rCurrentProcessInfo The current process info instance
        */
        void GetDofList(
            DofsVectorType& rElementalDofList,
            const ProcessInfo& rCurrentProcessInfo
            ) const override;

        /**
        * This functions calculates both the RHS and the LHS
        * @param rLeftHandSideMatrix: The LHS
        * @param rRightHandSideVector: The RHS
        * @param rCurrentProcessInfo: The current process info instance
        * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
        * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
        */
        void CalculateAll(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo,
            const bool CalculateStiffnessMatrixFlag,
            const bool CalculateResidualVectorFlag
            );

        /* @brief This function computes the determinants of jacobian for
         *         the initial configuration. Required for conservative integration.
         * @param rGeometry: corresponding geometry
         * @param rDeterminantOfJacobian: output determinants
         */
        void DeterminantOfJacobianInitial(
            const GeometryType& rGeometry,
            Vector& rDeterminantOfJacobian);

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "\"SupportLagrangeCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SupportLagrangeCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override {
            pGetGeometry()->PrintData(rOStream);
        }

        ///@}

    private:
        ///@name private variables
        ///@{

        const double shape_function_tolerance = 1e-6;

        ///@}
        ///@name private operations
        ///@{

        /* @brief checks all shape functions and
         *        returns the number of non zero nodes.
         */
        SizeType GetNumberOfNonZeroNodes() const;

        ///@}
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

    }; // Class SupportLagrangeCondition

}  // namespace Kratos.

#endif // KRATOS_SUPPORT_LAGRANGE_CONDITION_H_INCLUDED  defined