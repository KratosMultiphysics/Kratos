//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

#if !defined(KRATOS_SUPPORT_NITSCHE_NON_SYMMETRIC_CONDITION_H_INCLUDED )
#define  KRATOS_SUPPORT_NITSCHE_NON_SYMMETRIC_CONDITION_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/constitutive_law.h"

// External includes

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/iga_flags.h"

#include "geometries/coupling_geometry.h"

namespace Kratos
{

/// Nitsche-type based support condition.

    class SupportNonSymmetricNitscheCondition
        : public Condition
    {
    protected:

        /// Internal variables used for metric transformation
        struct KinematicVariables
        {
            // covariant metric
            array_1d<double, 3> a_ab_covariant;

            //base vector 1
            array_1d<double, 3> a1;
            //base vector 2
            array_1d<double, 3> a2;
            //base vector 3 normalized
            array_1d<double, 3> a3;
            //not-normalized base vector 3
            array_1d<double, 3> a3_tilde;
        
            //differential area
            double dA;

            //the tangent to the surface vector 
            array_1d<double, 3> t;
            //the normal to the surface boundary vector
            array_1d<double, 3> n;
            //the normal to the surface boundary vector in contravariant basis
            array_1d<double, 2> n_contravariant;

            /**
            * The default constructor
            * @param Dimension: The size of working space dimension
            */
            KinematicVariables(SizeType Dimension)
            {
                noalias(a_ab_covariant) = ZeroVector(Dimension);

                noalias(a1) = ZeroVector(Dimension);
                noalias(a2) = ZeroVector(Dimension);
                noalias(a3) = ZeroVector(Dimension);

                noalias(a3_tilde) = ZeroVector(Dimension);

                noalias(n) = ZeroVector(Dimension);
                noalias(n_contravariant) = ZeroVector(2);
                noalias(t) = ZeroVector(Dimension);

                dA = 1.0;
            }
        };

        /**
        * Internal variables used in the constitutive equations
        */
        struct ConstitutiveVariables
        {
            Vector StrainVector;
            Vector StressVector;
            Matrix ConstitutiveMatrix;

            /**
            * @param StrainSize: The size of the strain vector in Voigt notation
            */
            ConstitutiveVariables(SizeType StrainSize)
            {
                StrainVector = ZeroVector(StrainSize);
                StressVector = ZeroVector(StrainSize);
                ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
            }
        };
        
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of SupportNonSymmetricNitscheCondition
        KRATOS_CLASS_POINTER_DEFINITION(SupportNonSymmetricNitscheCondition);

        /// Size types
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor with Id and geometry
        SupportNonSymmetricNitscheCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {};

        /// Constructor with Id, geometry and property
        SupportNonSymmetricNitscheCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Condition(NewId, pGeometry, pProperties)
        {};

        /// Default constructor
        SupportNonSymmetricNitscheCondition()
            : Condition()
        {};

        /// Destructor.
        virtual ~SupportNonSymmetricNitscheCondition() = default;

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
            return Kratos::make_intrusive<SupportNonSymmetricNitscheCondition>(
                NewId, pGeom, pProperties);
        };

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive< SupportNonSymmetricNitscheCondition >(
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
            const ProcessInfo& rCurrentProcessInfo) override
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
            const ProcessInfo& rCurrentProcessInfo) override
        {
            if (rCurrentProcessInfo[BUILD_LEVEL] == 2)
            {
                CalculateNitscheStabilizationMatrix(rLeftHandSideMatrix, rRightHandSideVector,
                    rCurrentProcessInfo);
            }
            else
            {
                CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                    rCurrentProcessInfo, true, true);
            }
            
        }

        /**
         * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable rDestinationVariable.
         * @param rRHSVector input variable containing the RHS vector to be assembled
         * @param rRHSVariable variable describing the type of the RHS vector to be assembled
         * @param rDestinationVariable variable in the database to which the rRHSvector will be assembled
         * @param rCurrentProcessInfo The current process info instance
         */
        void AddExplicitContribution(const VectorType& rRHS,
            const Variable<VectorType>& rRHSVariable,
            const Variable<array_1d<double,3> >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo
        ) override;

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

        /// Calculates left (K) and right (u) hand sides
        void CalculateAll(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo,
            const bool CalculateStiffnessMatrixFlag,
            const bool CalculateResidualVectorFlag
        );

        void CalculateNitscheStabilizationMatrix(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo
        );

        void DeterminantOfJacobianInitial(
            const GeometryType& rGeometry,
            Vector& rDeterminantOfJacobian);

        ///@}
        ///@name Base Class Operations
        ///@{

        void GetValuesVector(
            Vector& rValues,
            int Step = 0) const override;

        enum class ConfigurationType {
            Current,
            Reference
        };

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "\"SupportNonSymmetricNitscheCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SupportNonSymmetricNitscheCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override {
            pGetGeometry()->PrintData(rOStream);
        }

        ///@}

    private:

        ///@name Member Variables
        ///@{

        // Components of the metric coefficient tensor on the contravariant basis
        std::vector<array_1d<double, 3>> m_A_ab_covariant_vector;

        // Determinant of the geometrical Jacobian.
        Vector m_dA_vector;

        /* Transformation the strain tensor from the curvilinear system
        *  to the local cartesian in voigt notation including a 2 in the
        *  shear part. */
        std::vector<Matrix> m_T_vector;

        /* Transformation the stress tensor from the local cartesian 
        *  to the curvilinear system in voigt notation. */
        std::vector<Matrix> m_T_hat_vector;

        // The reference contravariant basis vector
        std::vector<array_1d< array_1d<double, 3>,2>> m_reference_contravariant_base;

        // The normal to the boundary vector
        std::vector<array_1d<double, 2>> m_n_contravariant_vector;

        void CalculateKinematics(
            IndexType IntegrationPointIndex,
            KinematicVariables& rKinematicVariables,
            const Matrix& rShapeFunctionGradientValues, const ConfigurationType& rConfiguration);

        // Computes transformation
        void CalculateTransformation(
            const KinematicVariables& rKinematicVariables,
            Matrix& rT, Matrix& rT_hat, array_1d<array_1d<double, 3>,2>& rReferenceContraVariantBase);
        
        /// Traction-related functions
        // Computes the traction
        void CalculateTraction(
            IndexType IntegrationPointIndex,
            array_1d<double, 3>& rTraction,
            const KinematicVariables& rActualKinematic,
            ConstitutiveVariables& rThisConstitutiveVariablesMembrane);
        
        // Computes the first derivatives of the stresses at the covariant bases
        void CalculateFirstVariationStressCovariant(
            IndexType IntegrationPointIndex,
            Matrix& rFirstVariationStressCovariant,
            const KinematicVariables& rActualKinematic,
            ConstitutiveVariables& rThisConstitutiveVariablesMembrane);

        // Computes the first derivatives of traction
        void CalculateFirstVariationTraction(
            IndexType IntegrationPointIndex,
            Matrix& rFirstVariationTraction,
            Matrix& rFirstVariationStressCovariant,
            const KinematicVariables& rActualKinematic,
            ConstitutiveVariables& rThisConstitutiveVariablesMembrane);

        // Computes the necessary products needed for second derivatives of traction
        void CalculateSecondVariationTractionProduct(
            IndexType IntegrationPointIndex,
            Matrix& rPi,
            const KinematicVariables& rActualKinematic,
            ConstitutiveVariables& rThisConstitutiveVariablesMembrane);

        // Computes the second derivatives of traction
        void CalculateSecondVariationTraction(
            IndexType IntegrationPointIndex,
            Matrix& rSecondVariationTraction,
            const KinematicVariables& rActualKinematic,
            Matrix& rFirstVariationStressCovariant, 
            array_1d<double, 3>& rDisplacement,
            array_1d<double, 3>& rSecondVariationTractionProduct);
        
        /**
        * This functions updates the constitutive variables
        * @param rActualMetric: The actual metric
        * @param rThisConstitutiveVariablesMembrane: The constitutive variables to be calculated
        * @param rValues: The CL parameters
        * @param ThisStressMeasure: The stress measure considered
        */
        void CalculateConstitutiveVariables(
            IndexType IntegrationPointIndex,
            KinematicVariables& rActualMetric,
            ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
            ConstitutiveLaw::Parameters& rValues,
            const ConstitutiveLaw::StressMeasure ThisStressMeasure
        );

        /**
        * This functions calculates the matrix for prestress transformation
        * @param rActualKinematic: The actual metric
        * @param rPrestresstransVariables: 
        */
        void CalculateTransformationPrestress(
            Matrix& rTransformationPrestress,
            const KinematicVariables& rActualKinematic
        );

        ///@}
        ///@name Serialization
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
            rSerializer.save("A_ab_covariant_vector", m_A_ab_covariant_vector);
            rSerializer.save("dA_vector", m_dA_vector);
            rSerializer.save("T_vector", m_T_vector);
            rSerializer.save("reference_contravariant_base", m_reference_contravariant_base);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
            rSerializer.load("A_ab_covariant_vector", m_A_ab_covariant_vector);
            rSerializer.load("dA_vector", m_dA_vector);
            rSerializer.load("T_vector", m_T_vector);
            rSerializer.load("reference_contravariant_base", m_reference_contravariant_base);
        }

        ///@}

    }; // Class SupportNonSymmetricNitscheCondition

}  // namespace Kratos.

#endif // KRATOS_SUPPORT_NITSCHE_CONDITION_H_INCLUDED  defined 
