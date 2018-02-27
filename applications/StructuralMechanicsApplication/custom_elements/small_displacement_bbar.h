// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_SMALL_DISPLACEMENT_BBAR_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_BBAR_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Total Lagrangian element for 2D and 3D geometries.

/**
 * Implements a total Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 2D and 3D
 */

    class SmallDisplacementBbar
            : public SmallDisplacement
    {
    public:
        ///@name Type Definitions
        ///@{
        ///Reference type definition for constitutive laws
        typedef ConstitutiveLaw ConstitutiveLawType;
        ///Pointer type for constitutive laws
        typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
        ///Type definition for integration methods
        typedef GeometryData::IntegrationMethod IntegrationMethod;

        /// Counted pointer of SmallDisplacementStrElement
        KRATOS_CLASS_POINTER_DEFINITION(SmallDisplacementBbar);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        SmallDisplacementBbar(IndexType NewId, GeometryType::Pointer pGeometry);
        SmallDisplacementBbar(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /// Destructor.
        virtual ~SmallDisplacementBbar();

        ///@}
        ///@name Operators
        ///@{
        ///@}
        ///@name Operations
        ///@{
        /**
         * Returns the currently selected integration method
         * @return current integration method selected
         */
        //TODO: ADD THE OTHER CREATE FUNCTION
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /**
         * Calculate a Matrix Variable on the Element Constitutive Law
         * @param rVariable: The variable we want to get
         * @param rOutput: The values obtained int the integration points
         * @param rCurrentProcessInfo: the current process info instance
         */
        void CalculateOnIntegrationPoints(
                const Variable<Matrix>& rVariable,
                std::vector<Matrix>& rOutput,
                const ProcessInfo& rCurrentProcessInfo
        ) override;

        /**
         * Calculate a Vector Variable on the Element Constitutive Law
         * @param rVariable: The variable we want to get
         * @param rOutput: The values obtained int the integration points
         * @param rCurrentProcessInfo: the current process info instance
         */
        void CalculateOnIntegrationPoints(
                const Variable<Vector>& rVariable,
                std::vector<Vector>& rOutput,
                const ProcessInfo& rCurrentProcessInfo
        ) override;

        /**
         * Calculate a double Variable on the Element Constitutive Law
         * @param rVariable: The variable we want to get
         * @param rOutput: The values obtained int the integration points
         * @param rCurrentProcessInfo: the current process info instance
         */
        void CalculateOnIntegrationPoints(
                const Variable<double>& rVariable,
                std::vector<double>& rOutput,
                const ProcessInfo& rCurrentProcessInfo
        ) override;

        /**
         * Called at the end of eahc solution step
         * @param rCurrentProcessInfo: the current process info instance
         */
        void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    protected:

        /**
         * Internal variables used in the kinematic calculations
         */
        struct KinematicVariables
        {
            Vector  N;
            Matrix  B;
            Vector Bh;
            double  detF;
            Matrix  F;
            double  detJ0;
            Matrix  J0;
            Matrix  InvJ0;
            Matrix  DN_DX;

            /**
             * The default constructor
             * @param StrainSize: The size of the strain vector in Voigt notation
             * @param Dimension: The size of the strain vector in Voigt notation
             * @param NumberOfNodes: The size of the strain vector in Voigt notation
             */
            KinematicVariables(
                    const unsigned int& StrainSize,
                    const unsigned int& Dimension,
                    const unsigned int& NumberOfNodes
            )
            {
                detF = 1.0;
                detJ0 = 1.0;
                N = ZeroVector(NumberOfNodes);
                B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
                Bh = ZeroVector(Dimension * NumberOfNodes);
                F = IdentityMatrix(Dimension);
                DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
                J0 = ZeroMatrix(Dimension, Dimension);
                InvJ0 = ZeroMatrix(Dimension, Dimension);
            }
        };

        ///@name Protected static Member Variables
        ///@{
        ///@}
        ///@name Protected member Variables
        ///@{

        ///@}
        ///@name Protected Operators
        ///@{
        SmallDisplacementBbar() : SmallDisplacement()
        {
        }

        /**
         * This functions updates the kinematics variables
         * @param rThisKinematicVariables: The kinematic variables to be calculated
         * @param PointNumber: The integration point considered
         */

        /**
         * Calculation of the RHS
         */
        void CalculateAndAddResidualVector(
                VectorType& rRightHandSideVector,
                const KinematicVariables& rThisKinematicVariables,
                const ProcessInfo& CurrentProcessInfo,
                const Vector& BodyForce,
                const Vector& StressVector,
                const double IntegrationWeight
        );

        /**
         * Calculation of B standard matrices and other parameters
         */
        void CalculateKinematicVariables(
                KinematicVariables& rThisKinematicVariables,
                const unsigned int PointNumber,
                const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        );

        /**
        * This functions updates the constitutive variables
        * @param rThisKinematicVariables: The kinematic variables to be calculated
        * @param rThisConstitutiveVariables: The constitutive variables
        * @param rValues: The CL parameters
        * @param PointNumber: The integration point considered
        * @param IntegrationPoints: The list of integration points
        * @param ThisStressMeasure: The stress measure considered
        * @param Displacements: The displacements vector
        */
        void CalculateConstitutiveVariables(
                KinematicVariables& rThisKinematicVariables,
                ConstitutiveVariables& rThisConstitutiveVariables,
                ConstitutiveLaw::Parameters& rValues,
                const unsigned int PointNumber,
                const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                const ConstitutiveLaw::StressMeasure ThisStressMeasure,
                const Vector Displacements
        );

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
                ProcessInfo& rCurrentProcessInfo,
                const bool CalculateStiffnessMatrixFlag,
                const bool CalculateResidualVectorFlag
        ) override;

        /**
         * Calculation of the Deformation Matrix B
         * @param B: The deformation matrix
         * @param DN_DX: The derivatives of the shape functions
         */
        virtual void CalculateB(
                Matrix& rB,
                const Matrix& DN_DX,
                const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                const unsigned int PointNumber
        );


        /**
         * This functions updates the kinematics variables
         * @param rThisKinematicVariables: The kinematic variables to be calculated
         * @param PointNumber: The integration point considered
         */
        void CalculateKinematicVariablesBbar(
                KinematicVariables& rThisKinematicVariables,
                const unsigned int PointNumber,
                const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        );


        /**
         * Calculation of the Deformation Matrix Bbar
         * @param B: The deformation matrix
         * @param DN_DX: The derivatives of the shape functions
         */
        virtual void CalculateMatrixBbar(
                Matrix& rB,
                Vector& rBh,
                const Matrix& DN_DX,
                const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                const unsigned int PointNumber
        );

        // Compute Bbar components
        /**
         * This functions updates the kinematics variables
         * @param rThisKinematicVariables: The kinematic variables to be calculated
         */
        void CalculateHydrostaticDeformationMatrix(KinematicVariables& rThisKinematicVariables);

    private:

        friend class Serializer;

        // A private default constructor necessary for serialization

        virtual void save(Serializer& rSerializer) const override;

        virtual void load(Serializer& rSerializer) override;

        ///@name Private Inquiry
        ///@{
        ///@}
        ///@name Un accessible methods
        ///@{
        /// Assignment operator.
        //SmallDisplacementStrElement& operator=(const SmallDisplacementStrElement& rOther);
        /// Copy constructor.
        //SmallDisplacementStrElement(const SmallDisplacementStrElement& rOther);
        ///@}

    }; // Class SmallDisplacementBbar

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_SMALL_DISPLACEMENT_BBAR_H_INCLUDED defined
