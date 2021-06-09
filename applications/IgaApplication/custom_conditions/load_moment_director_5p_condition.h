//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_MOMENT_LOAD_DIRECTOR_5p_CONDITION_H_INCLUDED )
#define  KRATOS_MOMENT_LOAD_DIRECTOR_5p_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "iga_application_variables.h"

#include "includes/condition.h"

namespace Kratos
{
    /// Condition for moment loads for the 5p shell based on directors
    class LoadMomentDirector5pCondition final
        : public Condition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer definition of LoadCondition
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LoadMomentDirector5pCondition);

        /// Size types
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef Geometry<Node<3>> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor with Id and geometry
        LoadMomentDirector5pCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {};

        /// Constructor with Id, geometry and property
        LoadMomentDirector5pCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Condition(NewId, pGeometry, pProperties)
        {};

        /// Default constructor
        LoadMomentDirector5pCondition() : Condition()
        {};

        /// Destructor
        ~LoadMomentDirector5pCondition() = default;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties
        ) const final
        {
            return Kratos::make_intrusive<LoadMomentDirector5pCondition>(
                NewId, pGeom, pProperties);
        };

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties
        ) const final
        {
            return Kratos::make_intrusive<LoadMomentDirector5pCondition>(
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
            const ProcessInfo& rCurrentProcessInfo) final
        {
            MatrixType left_hand_side_matrix = Matrix(0, 0);

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
            const ProcessInfo& rCurrentProcessInfo) final
        {
            VectorType right_hand_side_vector = Vector(0);

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
            const ProcessInfo& rCurrentProcessInfo) final
        {
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
        ) const final;

        /**
        * @brief Sets on rConditionDofList the degrees of freedom of the considered element geometry
        * @param rElementalDofList The vector containing the dof of the element
        * @param rCurrentProcessInfo The current process info instance
        */
        void GetDofList(
            DofsVectorType& rElementalDofList,
           const ProcessInfo& rCurrentProcessInfo
        ) const final;

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

        void DeterminantOfJacobianInitial(
            const GeometryType& rGeometry,
            Vector& rDeterminantOfJacobian);

        array_1d<double, 3> calculateMomentLoadTimesDirectorTestFunction(
            const GeometryType& rGeometry,
            const Matrix& r_N,
            const IndexType& point_number,
            const array_1d<double, 3>& momentload);

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const final
        {
            std::stringstream buffer;
            buffer << "\"LoadCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const final
        {
            rOStream << "\"LoadCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const final
        {
            pGetGeometry()->PrintData(rOStream);
        }

        ///@}

    private:
        ///@name Serialization
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const final
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        }

        virtual void load(Serializer& rSerializer) final
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        }

        ///@}

    }; // Class LoadMomentDirector5pCondition

}  // namespace Kratos.

#endif // KRATOS_MOMENT_LOAD_DIRECTOR_5p_CONDITION_H_INCLUDED  defined