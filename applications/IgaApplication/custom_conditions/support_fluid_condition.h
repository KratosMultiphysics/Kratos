//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

# pragma once

// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
    /// Condition for penalty support condition
    class SupportFluidCondition
        : public Condition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer definition of SupportFluidCondition
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SupportFluidCondition);

        /// Size types
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        /// Type for shape function derivatives container
        typedef Kratos::Matrix ShapeDerivativesType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor with Id and geometry
        SupportFluidCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {};

        /// Constructor with Id, geometry and property
        SupportFluidCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Condition(NewId, pGeometry, pProperties)
        {};

        /// Default constructor
        SupportFluidCondition() : Condition()
        {};

        /// Destructor
        virtual ~SupportFluidCondition() override
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
            return Kratos::make_intrusive<SupportFluidCondition>(
                NewId, pGeom, pProperties);
        };

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& rThisNodes,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive<SupportFluidCondition>(
                NewId, GetGeometry().Create(rThisNodes), pProperties);
        };

        void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

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
            VectorType right_hand_side_vector;

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
            buffer << "\"SupportFluidCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SupportFluidCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override
        {
            pGetGeometry()->PrintData(rOStream);
        }

        void GetSolutionCoefficientVector(Vector& rValues) const;

        ///@}

    protected:

        void InitializeMaterial();

        ConstitutiveLaw::Pointer mpConstitutiveLaw; /// The pointer containing the constitutive laws

        /**
         * Internal variables used in the kinematic calculations
         */
        struct ConstitutiveVariables
        {
            ConstitutiveLaw::StrainVectorType StrainVector;
            ConstitutiveLaw::StressVectorType StressVector;
            ConstitutiveLaw::VoigtSizeMatrixType D;

            /**
             * The default constructor
             * @param StrainSize The size of the strain vector in Voigt notation
             */
            ConstitutiveVariables(const SizeType StrainSize)
            {
                if (StrainVector.size() != StrainSize)
                    StrainVector.resize(StrainSize);

                if (StressVector.size() != StrainSize)
                    StressVector.resize(StrainSize);

                if (D.size1() != StrainSize || D.size2() != StrainSize)
                    D.resize(StrainSize, StrainSize);

                noalias(StrainVector) = ZeroVector(StrainSize);
                noalias(StressVector) = ZeroVector(StrainSize);
                noalias(D)            = ZeroMatrix(StrainSize, StrainSize);
            }
        };    

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

        /**
         * @brief Initialize member variables
         * 
         * This function initializes the member variables of the condition.
         * It computes the basis function order and penalty factor based on the geometry.
         */
        void InitializeMemberVariables();

        void CalculateB(
            Matrix& rB,
            const ShapeDerivativesType& r_DN_DX) const;

        void ApplyConstitutiveLaw(
            const Matrix& rB, 
            ConstitutiveLaw::Parameters& rValues,
            ConstitutiveVariables& rConstitutiveVariables) const;

        // member variables
        unsigned int mDim;
        double mPenalty;
        IndexType mBasisFunctionsOrder;

        ///@}

    }; // Class SupportFluidCondition

}  // namespace Kratos.