//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Michael Breitenberger
//                   Riccardo Rossi
//

#if !defined(KRATOS_NITSCHE_COUPLING_CONDITION_H_INCLUDED )
#define  KRATOS_NITSCHE_COUPLING_CONDITION_H_INCLUDED

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
class NitscheCouplingCondition
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
        //the normal to the surface boundary vector
        array_1d<double, 3> n;
        //the normal to the surface boundary vector in contravariant basis
        array_1d<double, 2> n_contravariant;
        //the tangent to the surface vector 
        array_1d<double, 3> t;

        //differential area
        double dA;

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


    struct PrestresstransVariables
    {
        Matrix Tpre;
        
        PrestresstransVariables(SizeType Dimension)
        {
           Matrix Tpre = ZeroMatrix(3, 3);
        }
    };
    
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of NitscheCouplingCondition
    KRATOS_CLASS_POINTER_DEFINITION(NitscheCouplingCondition);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry
    NitscheCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {};

    /// Constructor with Id, geometry and property
    NitscheCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {};

    /// Default constructor
    NitscheCouplingCondition()
        : Condition()
    {};

    /// Destructor.
    virtual ~NitscheCouplingCondition() = default;

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
        return Kratos::make_intrusive<NitscheCouplingCondition>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< NitscheCouplingCondition >(
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
        ProcessInfo& rCurrentProcessInfo) override
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
        ProcessInfo& rCurrentProcessInfo) override
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
        ProcessInfo& rCurrentProcessInfo) override
    {
        if (rCurrentProcessInfo[BUILD_LEVEL] == 5)
        {
            CalculateQMatrix(rLeftHandSideMatrix, rRightHandSideVector,
                rCurrentProcessInfo);
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
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList The vector containing the dof of the element
    * @param rCurrentProcessInfo The current process info instance
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    ) override;

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
    );

    void CalculateQMatrix(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    );

    /**
    * @brief This is called during the assembling process in order to calculate the elemental damping matrix
    * @param rDampingMatrix The elemental damping matrix
    * @param rCurrentProcessInfo The current process info instance
    */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    ///@}
    ///@name Base Class Operations
    ///@{

    void Initialize() override;

    void GetValuesVector(
        Vector& rValues,
        int Step = 0) const override;
    
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    enum class ConfigurationType {
      Current,
      Reference
    };

    enum class PatchType {
      Master,
      Slave
    };

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"NitscheCouplingCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"NitscheCouplingCondition\" #" << Id();
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
    std::vector<array_1d<double, 3>> m_A_ab_covariant_vector_master;
    std::vector<array_1d<double, 3>> m_A_ab_covariant_vector_slave;

    // Determinant of the geometrical Jacobian.
    Vector m_dA_vector_master;
    Vector m_dA_vector_slave;

    /* Transformation the strain tensor from the curvilinear system
    *  to the local cartesian in voigt notation including a 2 in the
    *  shear part. */
    std::vector<Matrix> m_T_vector_master;
    std::vector<Matrix> m_T_vector_slave;

    /* Transformation the stress tensor from the local cartesian 
    *  to the curvilinear system in voigt notation. */
    std::vector<Matrix> m_T_hat_vector_master;
    std::vector<Matrix> m_T_hat_vector_slave;

    /* */
    std::vector<array_1d< array_1d<double, 3>,2>> m_reference_contravariant_base_master;
    std::vector<array_1d< array_1d<double, 3>,2>> m_reference_contravariant_base_slave;

    /// The vector containing the constitutive laws for all integration points.
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector_master;
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector_slave;

    // The normal to the boundary vector
    std::vector<array_1d<double, 2>> m_n_contravariant_vector_master;
    std::vector<array_1d<double, 2>> m_n_contravariant_vector_slave;


    /// Initialize Operations
    void InitializeMaterial();

    void CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const Matrix& rShapeFunctionGradientValues, const ConfigurationType& rConfiguration, const PatchType& rPatch);

    // Computes transformation
    void CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& rT, Matrix& rT_hat, array_1d<array_1d<double, 3>,2>& rReferenceContraVariantBase);
    
    void CalculateTraction(
        IndexType IntegrationPointIndex,
        array_1d<double, 3>& rTraction,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch);
    
    void CalculateDTraction(
        IndexType IntegrationPointIndex,
        Matrix& rDTraction,
        Matrix& rNCurvilinear,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch);

    void CalculateNCurvilinear(
        IndexType IntegrationPointIndex,
        Matrix& rNCurvilinear,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane, 
        const PatchType& rPatch);

    void CalculateDDTraction(
        IndexType IntegrationPointIndex,
        Matrix& rDDTraction,
        const KinematicVariables& rActualKinematic,
        Matrix& rNCurvilinear, 
        array_1d<double, 3>& rDispMaster,
        array_1d<double, 3>& rDispSlave,
        array_1d<double, 3>& rDDTractionProduct,
        array_1d<double, 3>& rDDTractionProductMasterSlave,
        const PatchType& rPatch);

    void CalculateDDTractionProduct(
        IndexType IntegrationPointIndex,
        Matrix& rPi,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        const PatchType& rPatch);
    
    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
        KinematicVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure,
        const PatchType& rPatch
    );

    /**
    * This functions calculates the matrix for prestress transformation
    * @param rActualKinematic: The actual metric
    * @param rPrestresstransVariables: 
    */
    void CalculateTransformationmatrixPrestress(
        const KinematicVariables& rActualKinematic,
        PrestresstransVariables& rPrestresstransVariables
        );

    ///@}
    ///@name Geometrical Functions
    ///@{

    void CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe);

    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        rSerializer.save("A_ab_covariant_vector_master", m_A_ab_covariant_vector_master);
        rSerializer.save("A_ab_covariant_vector_slave", m_A_ab_covariant_vector_slave);
        rSerializer.save("dA_vector_master", m_dA_vector_master);
        rSerializer.save("dA_vector_slave", m_dA_vector_slave);
        rSerializer.save("T_vector_master", m_T_vector_master);
        rSerializer.save("T_vector_slave", m_T_vector_slave);
        rSerializer.save("constitutive_law_vector_master", mConstitutiveLawVector_master);
        rSerializer.save("constitutive_law_vector_slave", mConstitutiveLawVector_slave);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        rSerializer.load("A_ab_covariant_vector_master", m_A_ab_covariant_vector_master);
        rSerializer.load("A_ab_covariant_vector_slave", m_A_ab_covariant_vector_slave);
        rSerializer.load("dA_vector_master", m_dA_vector_master);
        rSerializer.load("dA_vector_slave", m_dA_vector_slave);
        rSerializer.load("T_vector_master", m_T_vector_master);
        rSerializer.load("T_vector_slave", m_T_vector_slave);
        rSerializer.load("constitutive_law_vector_master", mConstitutiveLawVector_master);
        rSerializer.load("constitutive_law_vector_slave", mConstitutiveLawVector_slave);
    }

    ///@}

}; // Class NitscheCouplingCondition

}  // namespace Kratos.

#endif // KRATOS_NITSCHE_COUPLING_CONDITION_H_INCLUDED  defined 


