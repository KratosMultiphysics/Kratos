//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Ilaria Iaconeta
//                  Bodhinanda Chandra
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
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

/// Large Displacement Lagrangian Element for 3D and 2D geometries. (base class)

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */

class UpdatedLagrangianElement
    : public Element
{
public:

    ///@name Type Definitions
    ///@{

    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;

    /// Counted pointer of LargeDisplacementElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UpdatedLagrangianElement );

    ///@}

protected:

    /**
     * Parameters to be used in the Element
     * for optimized memory allocation.
     */
    struct KinematicVariables
    {
    public:
        // DeterminantF
        double  detF;
        // DeterminantFT
        double  detFT;
        // Deformation Gradient
        Matrix  F;
        // Shape Function Derivatives in global space
        Matrix  DN_DX;

        KinematicVariables(
            SizeType WorkingSpaceDimension,
            SizeType NumberOfNodes)
        {
            // The size of the strain vector in Voigt notation
            SizeType StrainSize = 3;
            if (WorkingSpaceDimension == 3)
                StrainSize = 6;

            detF = 1.0;
            detFT = 1.0;

            F = ZeroMatrix(WorkingSpaceDimension, WorkingSpaceDimension);

            DN_DX = ZeroMatrix(NumberOfNodes, WorkingSpaceDimension);
        }
    };


    /*
    * Constitutive variables hold the strain vector, the stress vector
    * and the constitutive matrix with according memory allocation.
    */
    struct ConstitutiveVariables
    {
    public:
        Vector  StrainVector;
        Vector  StressVector;

        Matrix  ConstitutiveMatrix;

        /**
        * The default constructor to allocate memory.
        * @param WorkingSpaceDimension of the simulation
        */
        ConstitutiveVariables(const SizeType& WorkingSpaceDimension)
        {
            // The size of the strain vector in Voigt notation
            SizeType StrainSize = 3;
            if (WorkingSpaceDimension == 3)
                StrainSize = 6;

            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);

            ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
        }

    };

public:
    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangianElement() {
        mDeterminantF0 = 1.0;
        mFinalizedStep = false;
    };

    /// Default constructors
    UpdatedLagrangianElement(
        IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
        mDeterminantF0 = 1.0;
        mFinalizedStep = false;
    };

    /// Default constructors
    UpdatedLagrangianElement(
        IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
        mDeterminantF0 = 1.0;
        mFinalizedStep = false;
    }

    ///Copy constructor
    UpdatedLagrangianElement(UpdatedLagrangianElement const& rOther)
        : Element(rOther)
        , mDeformationGradientF0(rOther.mDeformationGradientF0)
        , mDeterminantF0(rOther.mDeterminantF0)
        , mpConstitutiveLaw(rOther.mpConstitutiveLaw)
        , mFinalizedStep(rOther.mFinalizedStep)
    {
    }

    /// Destructor.
    ~UpdatedLagrangianElement() override {};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianElement& operator=(UpdatedLagrangianElement const& rOther)
    {
        Element::operator=(rOther);

        mDeformationGradientF0.clear();
        mDeformationGradientF0 = rOther.mDeformationGradientF0;

        mDeterminantF0 = rOther.mDeterminantF0;
        mpConstitutiveLaw = rOther.mpConstitutiveLaw;

        return *this;
    }

    ///@}
    ///@name Create and Clone
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<UpdatedLagrangianElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Element::Pointer(new UpdatedLagrangianElement(
            NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Clones the selected element variables, creating a new one
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        UpdatedLagrangianElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

        NewElement.mpConstitutiveLaw = mpConstitutiveLaw->Clone();

        NewElement.mDeformationGradientF0 = mDeformationGradientF0;

        NewElement.mDeterminantF0 = mDeterminantF0;

        return Element::Pointer(new UpdatedLagrangianElement(NewElement));
    }

    ///@}
    ///@name Dofs
    ///@{

    /// Sets on rElementalDofList the degrees of freedom of the element
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    /// Sets on rResult the ID's of the element degrees of freedom
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Life Time Operations
    ///@{

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Allocation of Contribution Matrices
    ///@{

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();

        const SizeType matrix_size = number_of_nodes * dimension;   //number of degrees of freedom

        if (rLeftHandSideMatrix.size1() != matrix_size)
            rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(matrix_size, matrix_size); //resetting LHS

        if (rRightHandSideVector.size() != matrix_size)
            rRightHandSideVector.resize(matrix_size, false);

        rRightHandSideVector = ZeroVector(matrix_size);

        CalculateAll(
            rLeftHandSideMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            true, true);
    }

    /// Check if necessary at certain points!
    ///**
    // * this function provides a more general interface to the element.
    // * it is designed so that rLHSvariables and rRHSvariables are passed TO the element
    // * thus telling what is the desired output
    // * @param rLeftHandSideMatrices: container with the output left hand side matrices
    // * @param rLHSVariables: paramter describing the expected LHSs
    // * @param rRightHandSideVectors: container for the desired RHS output
    // * @param rRHSVariables: parameter describing the expected RHSs
    // */
    //void CalculateLocalSystem(std::vector< MatrixType >& rLeftHandSideMatrices,
    //                          const std::vector< Variable< MatrixType > >& rLHSVariables,
    //                          std::vector< VectorType >& rRightHandSideVectors,
    //                          const std::vector< Variable< VectorType > >& rRHSVariables,
    //                          ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental right hand side vector only
      * @param rRightHandSideVector: the elemental right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();

        //number of degrees of freedom
        const SizeType matrix_size = number_of_nodes * dimension;

        Matrix LeftHandSideMatrix; //resetting LHS

        if (rRightHandSideVector.size() != matrix_size)
            rRightHandSideVector.resize(matrix_size, false);

        rRightHandSideVector = ZeroVector(matrix_size);

        CalculateAll(
            LeftHandSideMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            false, true);
    }

    /// Check if necessary at certain points!
    ///**
    // * this function provides a more general interface to the element.
    // * it is designed so that rRHSvariables are passed TO the element
    // * thus telling what is the desired output
    // * @param rRightHandSideVectors: container for the desired RHS output
    // * @param rRHSVariables: parameter describing the expected RHSs
    // */
    //void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
    //                            const std::vector< Variable< VectorType > >& rRHSVariables,
    //                            ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side vector only
     * @param rLeftHandSideVector: the elemental left hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();

        //number of degrees of freedom
        const SizeType matrix_size = number_of_nodes * dimension;

        if (rLeftHandSideMatrix.size1() != matrix_size)
            rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(matrix_size, matrix_size);

        Vector RightHandSideVector;

        CalculateAll(
            rLeftHandSideMatrix,
            RightHandSideVector,
            rCurrentProcessInfo,
            true, false);
    }

    ///@}
    ///@name Operations for Dynamic
    ///@{

    /**
      * Computes the mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(MatrixType& rMassMatrix,
                             ProcessInfo& rCurrentProcessInfo) override;

    /**
      * Computes the damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& rCurrentProcessInfo) override;

    /// Sets nodal displacements on rValues
    void GetValuesVector(Vector& rValues, int Step = 0) override;

    /// Sets nodal velocitites on rValues
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    /// Sets nodal accelerations on rValues
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

    ///@}
    ///@name Access
    ///@{

    /// check completeness of input.
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MPM Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MPM Element #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        GetGeometry().PrintData(rOStream);
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{
    
    /// Container for historical total elastic deformation measure F0 = dx/dX
    Matrix mDeformationGradientF0;

    /// Container for the total deformation gradient determinants
    double mDeterminantF0;

    /// Container for constitutive law instances on each integration point
    ConstitutiveLaw::Pointer mpConstitutiveLaw;

    /// Finalize and Initialize label
    bool mFinalizedStep;

    ///@}
    ///@name Protected Operators
    ///@{

    /// Assembles required components for operation.
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag);

    ///@}
    ///@name Protected Stiffness and Force Contribution
    ///@{

    /// Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
    void UpdatedLagrangianElement::CalculateAndAddKuum(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rB,
        const Matrix& rConstitutiveMatrix,
        const double& rIntegrationWeight) const;

    /// Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
    void CalculateAndAddKuug(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rDN_DX,
        const Vector& rStressVector,
        const double& rIntegrationWeight) const;

    /// Calculation of the External Forces Vector. Fe = N * t + N * b
    virtual void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
            Vector& rVolumeForce,
            const double& rIntegrationWeight) const;

    /// Calculation of the Internal Forces Vector. Fi = B * sigma
    virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        const Matrix& rB,
        const Vector& rStressVector,
        const double& rIntegrationWeight);

    ///@}
    ///@name Protected Kinematics and Constitutive behavior
    ///@{

    /// Calculates Element Kinematics
    void CalculateKinematics(
        KinematicVariables& rKinematicVariables,
        ProcessInfo& rCurrentProcessInfo) const;

    /// Calculation of the Deformation Matrix B
    virtual void CalculateBMatrix(
        Matrix& rB,
        const Matrix& rDN_DX) const;

    /// Calculates Constitutive response
    void CalculateConstitutiveVariables(
        KinematicVariables& rDeformationGradient,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure) const;

    ///@}
    ///@name Protected Constitutive Law
    ///@{

    /// Initialize Material Properties on the Constitutive Law
    void InitializeMaterial();

    /// Reset the Constitutive Law Parameters
    void ResetConstitutiveLaw() override;

    ///@}
    ///@name Protected Integration Parameters
    ///@{

    /// Calculation of the Integration Weight
    virtual double CalculateIntegrationWeight(
        const KinematicVariables& rKinematicVariables);

    /// Calculation of the Volume Change of the Element
    virtual double& CalculateVolumeChange(double& rVolumeChange, KinematicVariables& rVariables);

    /// Calculation of the Volume Force of the Element
    virtual Vector& CalculateVolumeForce(Vector& rVolumeForce);

    ///@}
    ///@name Protected MPM Life Time Variables
    ///@{

    /// Calculation of the Current Displacement
    Matrix& SetCurrentDisplacement(
        Matrix& rCurrentDisp) const;

    /// Finalize Element Internal Variables
    virtual void FinalizeStepVariables(
        KinematicVariables & rVariables,
        ConstitutiveVariables& rConstitutiveVariables,
        const ProcessInfo& rCurrentProcessInfo);

    /// Update the position of the MP or Gauss point when Finalize Element Internal Variables is called
    virtual void UpdateGaussPoint(
        const ProcessInfo& rCurrentProcessInfo);

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
        rSerializer.save("ConstitutiveLawVector", mpConstitutiveLaw);
        rSerializer.save("DeformationGradientF0", mDeformationGradientF0);
        rSerializer.save("DeterminantF0", mDeterminantF0);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
        rSerializer.load("ConstitutiveLawVector", mpConstitutiveLaw);
        rSerializer.load("DeformationGradientF0", mDeformationGradientF0);
        rSerializer.load("DeterminantF0", mDeterminantF0);
    }

    ///@}

}; // Class UpdatedLagrangianElement

///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED defined
