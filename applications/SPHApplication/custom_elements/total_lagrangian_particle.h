// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#include "includes/element.h"
#include "sph_application_variables.h"
#include "custom_utilities/custom_kernels/kernel_factory.h"
#include "custom_utilities/compute_kernel_correction_utilities.h"
#include "custom_utilities/sph_element_utilities.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

using SizeType = std::size_t;

template<class TKernelType>
class KRATOS_API(SPH_APPLICATION) TotalLagrangianDisplacementParticle: public Element
{
protected:
    /**
     * Internal variables used in the kinematic calculations
     */
    struct KinematicVariables
    {
        VectorType W;
        MatrixType DW_DX;
        MatrixType B;
        double detF;
        MatrixType F;
        VectorType Displacement;

        /**
         * @brief Default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         * @param DomainSize The size of the problem domain
         * @param NumberOfNeighbours The number of neighbours of the particle
         */
        KinematicVariables(
            const SizeType StrainSize,
            const SizeType DomainSize,
            const SizeType NumberOfNeighbours
        )
        {
            W = ZeroVector(NumberOfNeighbours);
            DW_DX = ZeroMatrix(NumberOfNeighbours, DomainSize);
            B = ZeroMatrix(StrainSize, DomainSize * NumberOfNeighbours);
            detF = 1.0;
            F = IdentityMatrix(DomainSize);
            Displacement = ZeroVector(DomainSize * NumberOfNeighbours);
        }
    };
    
    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector;
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType C;

        /**
         * @brief Default constructor
         */
        ConstitutiveVariables(const SizeType StrainSize)
        {
            if (StrainVector.size() != StrainSize) StrainVector.resize(StrainSize);
            if (StressVector.size() != StrainSize) StressVector.resize(StrainSize);
            if (C.size1() != StrainSize || C.size2() != StrainSize) C.resize(StrainSize, StrainSize);

            noalias(StrainVector) = ZeroVector(StrainSize);
            noalias(StressVector) = ZeroVector(StrainSize);
            noalias(C) = ZeroMatrix(StrainSize, StrainSize);
        }
    };
public: 

    using BaseType = Element;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TotalLagrangianDisplacementParticle);

    // Constructor void 
    TotalLagrangianDisplacementParticle()
    {
    }

    // Constructor using an array of nodes 
    TotalLagrangianDisplacementParticle(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    // Constructor using an array of nodes with properties 
    TotalLagrangianDisplacementParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    { 
    }

    // Copy constructor
    TotalLagrangianDisplacementParticle(TotalLagrangianDisplacementParticle const& rOther)
        : BaseType(rOther),
        mThisConstitutiveLaw(rOther.mThisConstitutiveLaw)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TotalLagrangianDisplacementParticle>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TotalLagrangianDisplacementParticle>(NewId, pGeom, pProperties);
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     */
    Element::Pointer Clone( IndexType NewId, NodesArrayType const& rThisNodes) const override;
    
    /**
     * @brief Called to initialize the element
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Returns a vector that includes the values of the DoFs
     */
    virtual void GetNodalValuesVector(VectorType& rNodalValue) const;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     */
    void EquationIdVector(
        EquationIdVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;
    
    /**
     *  @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(VectorType& rValues, int step = 0) const override;

    /**
     * @brief Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(VectorType& rValues, int step = 0) const override;

    /**
     * @brief
     */
    void GlobalSizeVector(
        VectorType& rLocalVector,
        VectorType& rNodalValue,
        IndexType i
    )
    {
        IndexType dofs_per_node = rLocalVector.size();
        rLocalVector.clear();
        rLocalVector[0] = rNodalValue[dofs_per_node * i];
        rLocalVector[1] = rNodalValue[dofs_per_node * i + 1];
    };

    /**
     * @brief This is called during the assembling process in order to calculate the local system
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS matrix
     * @param rRightHandSideVector The RHS vector
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );
    
    /**
     * @brief This function is called to initialize the kinematic variables
     */
    virtual void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables, 
        const ProcessInfo& rCurrentProcessInfo
    );

    /**
     * @brief This function computes the deformation gradient of the particle
     */
    virtual void CalculateDeformationGradient(
        MatrixType& rF,
        MatrixType& rDW_DX,
        VectorType rW,
        const ProcessInfo& rProcessInfo
    );

    /**
     * @brief This function computes the deformation matrix B
     */
    void CalculateB(
        MatrixType& rB,
        const MatrixType& rF,
        const MatrixType& rDW_DX
    );

    /**
     * @brief This function is called to initialize the constitutive variables
     */
    virtual void CalculateConstitutiveVariables(
        ConstitutiveVariables& rThisConstitutiveVariables,
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );

    /**
     * @brief This function is called to set the variables for the constitutive law
     */
    virtual void SetConstitutiveLawVariables(
        ConstitutiveVariables& rThisConstitutiveVariables,
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveLaw::Parameters& rValues
    );
    
    /**
      * @brief Calculation of the Geometric Stiffness Matrix. Kg = dB * S
      * @param 
      * @param 
      */
    void CalculateAndAddKg(
        MatrixType& rLeftHandSideMatrix,
        const MatrixType& DW_DX,
        const VectorType& StressVector,
        const double IntegrationWeight
    ) const;

    /**
      * @brief Calculation of the Material Stiffness Matrix. Km = B^T * D *B
      * @param rLeftHandSideMatrix The local LHS of the element
      * @param B The deformationmmatrix
      * @param D The constitutive matrix
      * @param IntegrationWeight The integration weight of the corresponding Gauss point
      */
    void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& rConstitutiveMatrix,
        const double IntegrationWeight
    ) const; 

    /**
      * @brief Calculation of the right hand side
      * @param rRightHandSideVector The local component of the RHS due to external forces
      * @param rThisKinematicVariables The kinematic variables
      * @param rStressVector The vector containing the stress components
      */
    void CalculateAndAddResidualVector(
        VectorType& rRightHandSideVector,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const VectorType& rBodyForce,
        const Vector& rStressVector,
        const double IntegrationWeight
    ) const; 

    /**
     * @brief This function calculates the contribution of the external forces to the RHS
     */
    void CalculateAndAddExternelForcesContribution(
        const VectorType& rW,
        const ProcessInfo& rProcessInfo,
        const VectorType& rBodyForce,
        VectorType& rRHS,
        const double weight
    ) const;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix The elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental damping matrix
      * @param rDampingMatrix The elemental damping matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief These functions calculates the values of variables in the integrations points.
     * In SPH case coincide with the neighbouring particles 
     * @details These functions expect a std::vector of values for the specified variable type
     * @param rVariable This parameter selects the output 
     * @param SPH_KERNEL The function computes the kernel values in the neighbours  
     * @param SPH_KERNEL_GRADIENT The function computes the kernel gradient values in the neighbours  
     */

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

protected:

    ConstitutiveLaw::Pointer mThisConstitutiveLaw;

    /**
     * @brief This function sets the used constitutive laws
     */
    void SetConstitutiveLaw(const ConstitutiveLaw::Pointer rThisConstitutiveLaw)
    {
        mThisConstitutiveLaw = rThisConstitutiveLaw;
    }

    /**
     * @brief It initializes the material
     */
    void InitializeMaterial();

private:

};

}