//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Minas Apostolakis
//                   

#if !defined(KRATOS_COUPLING_SOLIDSHELL3P_NITSCHE_CONDITION_H_INCLUDED )
#define  KRATOS_COUPLING_SOLIDSHELL3P_NITSCHE_CONDITION_H_INCLUDED

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

/// Nitsche factor based coupling condition.
/** This condition can be used to apply continuity between different
*   discretizations with the nitsche approach.
*

*/
class CouplingSolidShell3pNitscheCondition
    : public Condition
{
protected:

    /// Internal variables used for metric transformation
    struct KinematicVariablesShell
    {
        // covariant metric
        array_1d<double, 3> a_ab_covariant;
        array_1d<double, 3> b_ab_covariant;

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

        // Hessian 
        Matrix Hessian = ZeroMatrix(3, 3);

        /**
        * The default constructor
        * @param Dimension: The size of working space dimension
        */
        KinematicVariablesShell(SizeType Dimension)
        {
            noalias(a_ab_covariant) = ZeroVector(Dimension);
            noalias(b_ab_covariant) = ZeroVector(Dimension);

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

    /// Internal variables used for metric transformation
    struct KinematicVariablesSolid
    {
        // covariant metric [g11, g22, g33, g12, g13, g23]
        array_1d<double, 6> a_ab_covariant;

        //base vector 1
        array_1d<double, 3> a1;
        //base vector 2
        array_1d<double, 3> a2;
        //base vector 3 normalized
        array_1d<double, 3> a3;

        ////differential area // TODO: Check if differential area is needed
        //double dA;

        //the tangent to the surface vector  // TODO: Check if differential area is needed
        //array_1d<double, 3> t;
        //the normal to the surface boundary vector
        array_1d<double, 3> n;
        //the normal to the surface boundary vector in contravariant basis 
        array_1d<double, 3> n_contravariant;

        /**
        * The default constructor
        * @param Dimension: The size of working space dimension
        */
        KinematicVariablesSolid(SizeType Dimension)
        {
            noalias(a_ab_covariant) = ZeroVector(Dimension);

            noalias(a1) = ZeroVector(Dimension);
            noalias(a2) = ZeroVector(Dimension);
            noalias(a3) = ZeroVector(Dimension);

            noalias(n) = ZeroVector(Dimension);
            noalias(n_contravariant) = ZeroVector(3);

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

    /// Counted pointer of CouplingSolidShell3pNitscheCondition
    KRATOS_CLASS_POINTER_DEFINITION(CouplingSolidShell3pNitscheCondition);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry
    CouplingSolidShell3pNitscheCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {};

    /// Constructor with Id, geometry and property
    CouplingSolidShell3pNitscheCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {};

    /// Default constructor
    CouplingSolidShell3pNitscheCondition()
        : Condition()
    {};

    /// Destructor.
    virtual ~CouplingSolidShell3pNitscheCondition() = default;

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
        return Kratos::make_intrusive<CouplingSolidShell3pNitscheCondition>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< CouplingSolidShell3pNitscheCondition >(
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
            /*CalculateNitscheStabilizationMatrix(rLeftHandSideMatrix, rRightHandSideVector,
                rCurrentProcessInfo);*/
        }
        else
        {
            CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                rCurrentProcessInfo, true, true);
        }
        
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

    /// Calculates left (K) and right (u) hand sides
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    //void CalculateNitscheStabilizationMatrix(
    //    MatrixType& rLeftHandSideMatrix,
    //    VectorType& rRightHandSideVector,
    //    const ProcessInfo& rCurrentProcessInfo
    //);

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
        buffer << "\"CouplingSolidShell3pNitscheCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"CouplingSolidShell3pNitscheCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

private:

    ///@name Member Variables
    ///@{
    
    // Shell base vactors and theta3
    std::vector<array_1d<double, 3>> _A1, _A2, _A3;
    std::vector<double> _theta3;
    // Components of the metric coefficient tensor on the contravariant basis
    std::vector<array_1d<double, 6>> m_A_ab_covariant_vector_master; // Solid has 6.
    std::vector<array_1d<double, 3>> m_A_ab_covariant_vector_slave;
    std::vector<array_1d<double, 3>> m_B_ab_covariant_vector_slave;

    // Determinant of the geometrical Jacobian.
    // Vector m_dA_vector_master;
    Vector m_dA_vector_slave;

    /* Transformation the strain tensor from the curvilinear system
    *  to the local cartesian in voigt notation including a 2 in the
    *  shear part. */
    // std::vector<Matrix> m_T_vector_master; for master Solid it is not needed
    std::vector<Matrix> m_T_vector_slave;

    /* Transformation the stress tensor from the local cartesian 
    *  to the curvilinear system in voigt notation. */
    //std::vector<Matrix> m_T_hat_vector_master; for master Solid it is not needed
    std::vector<Matrix> m_T_hat_vector_slave;

    std::vector<array_1d< array_1d<double, 3>,3>> m_reference_contravariant_base_master;
    std::vector<array_1d< array_1d<double, 3>,2>> m_reference_contravariant_base_slave;

    // The normal to the boundary vector
    std::vector<array_1d<double, 3>> m_n_contravariant_vector_master; // Solid -> ti = Sij nj : j = 1,2,3
    std::vector<array_1d<double, 2>> m_n_contravariant_vector_slave;

    double ComputeTheta3Shell(
        IndexType IntegrationPointIndex,
        const KinematicVariablesShell& rCurrentConfiguraitonKinematicVariablesShell);

    void OutOfPlaneDeformationFirstVariation(
        Matrix& OutOfPlaneDeformationWholeMatrix,
        const size_t& mat_size,
        const double theta,
        const array_1d<double, 3>& A1,
        const array_1d<double, 3>& A2,
        const Matrix& shape_functions_gradients_slave);

    array_1d<double, 3> Calculate_Phi_r_cross_A3(
        const array_1d<double, 3>& N_theta1,
        const array_1d<double, 3>& N_theta2,
        const array_1d<double, 3>& A1,
        const array_1d<double, 3>& A2);

    void CalculateKinematicsShell(
        IndexType IntegrationPointIndex,
        KinematicVariablesShell& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues, 
        const ConfigurationType& rConfiguration);

    void CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe) const;

    void CalculateTransposeInverseMatrix3x3(Matrix& TransposedInverted, const Matrix Input);

    void CalculateKinematicsSolid(
        IndexType IntegrationPointIndex,
        KinematicVariablesSolid& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues, const ConfigurationType& rConfiguration);

    // Computes transformation
    void CalculateTransformationShell(
        const KinematicVariablesShell& rKinematicVariables,
        Matrix& rT, Matrix& rT_hat, array_1d<array_1d<double, 3>,2>& rReferenceContraVariantBase); 

    // Traction-related functions
    void CalculateTractionShell(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariablesShell& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesShell);

    void CalculateTractionSolid(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariablesSolid& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesSolid);

    void CalculateFirstVariationStressCovariantShell(
        IndexType IntegrationPointIndex,
        const double& theta3,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariablesShell& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesShell
    );

    void CalculateFirstVariationStressCovariantSolid(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariablesSolid& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesSolid);
        
    void CalculateFirstVariationTractionShell(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationTraction,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariablesShell& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesShell);

    void CalculateFirstVariationTractionSolid(
        IndexType IntegrationPointIndex,
        Matrix& rFirstVariationTraction,
        Matrix& rFirstVariationStressCovariant,
        const KinematicVariablesSolid& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesSolid);

    void CalculateSecondVariationTractionSolid(
        IndexType IntegrationPointIndex,
        Matrix& rSecondVariationsTraction,
        const Matrix& rFirstVariationStressCovariant,
        const array_1d<double, 3>& rDisplacementMaster,
        const array_1d<double, 3>& rDisplacementSlave,
        const KinematicVariablesSolid& rActualKinematic,
        const ConstitutiveVariables& rThisConstitutiveVariablesSolid);

    void CalculateSecondVariationTractionShell(
        IndexType IntegrationPointIndex,
        Matrix& rProductSecondVariationTraction_Displacement,
        const double& theta3,
        const Matrix& rFirstVariationStressCovariant, 
        const array_1d<double, 3>& rDisplacementMaster,
        const array_1d<double, 3>& rDisplacementSlave,
        const KinematicVariablesShell& rActualKinematic,
        const ConstitutiveVariables& rThisConstitutiveVariablesShell);
    
    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariablesShell(
        IndexType IntegrationPointIndex,
        const double& theta3,
        KinematicVariablesShell& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );

    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
        void CalculateConstitutiveVariablesSolid(
            IndexType IntegrationPointIndex,
            KinematicVariablesSolid& rActualMetric,
            ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
            ConstitutiveLaw::Parameters& rValues,
            const ConstitutiveLaw::StressMeasure ThisStressMeasure
        );


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        rSerializer.save("A_ab_covariant_vector_master", m_A_ab_covariant_vector_master);
        rSerializer.save("A_ab_covariant_vector_slave", m_A_ab_covariant_vector_slave);
        //rSerializer.save("dA_vector_master", m_dA_vector_master);
        rSerializer.save("dA_vector_slave", m_dA_vector_slave);
        //rSerializer.save("T_vector_master", m_T_vector_master);
        rSerializer.save("T_vector_slave", m_T_vector_slave);
        //rSerializer.save("reference_contravariant_base_master", m_reference_contravariant_base_master);
        rSerializer.save("reference_contravariant_base_slave", m_reference_contravariant_base_slave);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        rSerializer.load("A_ab_covariant_vector_master", m_A_ab_covariant_vector_master);
        rSerializer.load("A_ab_covariant_vector_slave", m_A_ab_covariant_vector_slave);
        //rSerializer.load("dA_vector_master", m_dA_vector_master);
        rSerializer.load("dA_vector_slave", m_dA_vector_slave);
        //rSerializer.load("T_vector_master", m_T_vector_master);
        rSerializer.load("T_vector_slave", m_T_vector_slave);
        //rSerializer.load("reference_contravariant_base_master", m_reference_contravariant_base_master);
        rSerializer.load("reference_contravariant_base_slave", m_reference_contravariant_base_slave);
    }

    ///@}

}; // Class CouplingSolidShell3pNitscheCondition

}  // namespace Kratos.

#endif // KRATOS_COUPLING_NITSCHE_CONDITION_H_INCLUDED  defined 



