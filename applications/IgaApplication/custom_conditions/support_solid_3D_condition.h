//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_SUPPORT_SOLID_3D_CONDITION_H_INCLUDED )
#define  KRATOS_SUPPORT_SOLID_3D_CONDITION_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/condition.h"
#include "utilities/math_utils.h"
#include "includes/variables.h"


// External includes

// Project includes
#include "iga_application_variables.h"

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
    /// Condition for penalty support condition
    class SupportSolid3DCondition
        : public Condition
    {
    protected:

    /**
     * Internal variables used in the constitutive calculations
     */
        struct KinematicVariables
            {
                Vector  N;
                Matrix  B;
                double  detF;
                Matrix  F;
                double  detJ0;
                Matrix  J0;
                Matrix  InvJ0;
                Matrix  DN_DX;
                Vector Displacements;

                /**
                 * The default constructor
                 * @param StrainSize The size of the strain vector in Voigt notation
                 * @param Dimension The problem dimension: 2D or 3D
                 * @param NumberOfNodes The number of nodes in the element
                 */
                KinematicVariables(
                    const SizeType StrainSize,
                    const SizeType Dimension,
                    const SizeType NumberOfNodes
                    )
                {
                    detF = 1.0;
                    detJ0 = 1.0;
                    N = ZeroVector(NumberOfNodes);
                    B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
                    F = IdentityMatrix(Dimension);
                    DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
                    J0 = ZeroMatrix(Dimension, Dimension);
                    InvJ0 = ZeroMatrix(Dimension, Dimension);
                    Displacements = ZeroVector(Dimension * NumberOfNodes);
                }
            };
    IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods

    //std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws
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

    ///@name Protected static Member Variables
    ///@{
    void InitializeMaterial();


/**
     * @brief Gives the StressMeasure used
     */
    ConstitutiveLaw::StressMeasure GetStressMeasure() const;

    //@}
    ///@name Protected member Variables
    ///@{
    ConstitutiveLaw::Pointer mpConstitutiveLaw; /// The pointer containing the constitutive law

    ///@}

    /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     */
    void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        );



    /**
     * @brief This method returns if the element provides the strain
     */
    virtual bool UseElementProvidedStrain() const;    
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer definition of SupportSolid3DCondition
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SupportSolid3DCondition);

        /// Size types
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        ///@}
        ///@name Life Cycle
        ///@{
        
        void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

        /// Constructor with Id and geometry
        SupportSolid3DCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {};

        /// Constructor with Id, geometry and property
        SupportSolid3DCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Condition(NewId, pGeometry, pProperties)
        {};

        /// Default constructor
        SupportSolid3DCondition() : Condition()
        {};

        /// Destructor
        virtual ~SupportSolid3DCondition() override
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
            return Kratos::make_intrusive<SupportSolid3DCondition>(
                NewId, pGeom, pProperties);
        };

        const Matrix& ShapeFunctionsValues(IntegrationMethod ThisMethod) const
        {
            return GetGeometry().ShapeFunctionsValues(ThisMethod);
        }

        const GeometryType::IntegrationPointsArrayType  IntegrationPoints() const 
        {
            return GetGeometry().IntegrationPoints();
        }

        const GeometryType::IntegrationPointsArrayType  IntegrationPoints(IntegrationMethod ThisMethod) const
        {
            return GetGeometry().IntegrationPoints(ThisMethod);
        }
    // /**
    //  * @brief Returns the used integration method
    //  * @return default integration method of the used Geometry
    //  */
    // IntegrationMethod GetIntegrationMethod() const override
    // {
    //     return mThisIntegrationMethod;
    // }

    // /**
    // * element can be integrated using the GP provided by the geometry or custom ones
    // * by default, the base element will use the standard integration provided by the geom
    // * @return bool to select if use/not use GPs given by the geometry
    // */
    // bool UseGeometryIntegrationMethod() const
    // {
    //     return true;
    // }

    // const GeometryType::IntegrationPointsArrayType  IntegrationPoints() const 
    // {
    //     return GetGeometry().IntegrationPoints();
    // }

    // const GeometryType::IntegrationPointsArrayType  IntegrationPoints(IntegrationMethod ThisMethod) const
    // {
    //     return GetGeometry().IntegrationPoints(ThisMethod);
    // }

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive<SupportSolid3DCondition>(
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

            MatrixType left_hand_side_matrix = ZeroMatrix(mat_size, mat_size);

            CalculateAll(left_hand_side_matrix, rRightHandSideVector,
                rCurrentProcessInfo, false, true);
        }

        // /**
        // * @brief This is called during the assembling process in order
        // *        to calculate the condition left hand side matrix
        // * @param rLeftHandSideMatrix the condition left hand side matrix
        // * @param rCurrentProcessInfo the current process info
        // */
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

        void GetValuesVector(Vector& rValues) const;
        /* @brief This functions updates the kinematics variables
        * @param rThisKinematicVariables The kinematic variables to be calculated
        * @param PointNumber The integration point considered
        */
        void CalculateKinematicVariables(
            KinematicVariables& rThisKinematicVariables,
            const IndexType PointNumber,
            const GeometryType::IntegrationMethod& rIntegrationMethod
            );

        void CalculateB(Matrix& rB, Matrix& r_DN_DX) const;

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
            buffer << "\"SupportSolid3DCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SupportSolid3DCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override
        {
            pGetGeometry()->PrintData(rOStream);
        }

        ///@}
    
    /**
     * @brief Called at the end of eahc solution step
     * @param rCurrentProcessInfo the current process info instance
     */
        void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

        void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    


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

    }; // Class SupportPenaltyLaplacianCondition

}  // namespace Kratos.

#endif // SupportSolid3DCondition  defined
