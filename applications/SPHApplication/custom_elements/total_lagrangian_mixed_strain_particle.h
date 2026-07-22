// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#include "includes/element.h"
#include "sph_application_variables.h"
#include "custom_utilities/custom_kernels/kernel_factory.h"
#include "custom_utilities/compute_kernel_correction_utilities.h"
#include "custom_utilities/sph_element_utilities.h"

namespace Kratos
{
    
using SizeType = std::size_t;

template<class TKernelType>
class KRATOS_API(SPH_APPLICATION) TotalLagrangianMixedStrainParticle: public Element
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

    typedef GeometryData::IntegrationMethod IntegrationMethod; 

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TotalLagrangianMixedStrainParticle);

    // Constructor void 
    TotalLagrangianMixedStrainParticle()
    {
    }

    // Constructor using an array of nodes 
    TotalLagrangianMixedStrainParticle(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    // Constructor using an array of nodes with properties 
    TotalLagrangianMixedStrainParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod =  GeometryData::IntegrationMethod::GI_GAUSS_1;  ///
    }

    // Copy constructor
    TotalLagrangianMixedStrainParticle(TotalLagrangianMixedStrainParticle const& rOther)
        : BaseType(rOther),
        mThisConstitutiveLaw(rOther.mThisConstitutiveLaw),
        mThisIntegrationMethod(rOther.mThisIntegrationMethod)  ////
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TotalLagrangianMixedStrainParticle>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TotalLagrangianMixedStrainParticle>(NewId, pGeom, pProperties);
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
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    IntegrationMethod GetIntegrationMethod() const override
    {
        return mThisIntegrationMethod;   ////
    }

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation IDs
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
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
     * @brief
     */
    void GetValuesVector(Vector& rValues, int step = 0) const override; 
    
    /**
     * @brief Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(VectorType& rValues, int step = 0) const override;

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
      * @brief Calculation of the Geometric Stiffness Matrix. Kg = dB * S
      * @param StressVector The vector containing the stress components
      */
    void CalculateAndAddKg(
        MatrixType& rK,
        const MatrixType& DW_DX,
        const VectorType& StressVector,
        const double IntegrationWeight
    ) const;

    /**
      * @brief Calculation of the Material Stiffness Matrix. Km = B^T * D *B
      * @param rLeftHandSideMatrix The local LHS of the element
      * @param B The deformation matrix (Total Lagrangian Framework)
      * @param D The constitutive matrix
      * @param IntegrationWeight The integration weight of the corresponding Gauss point
      */
    void CalculateAndAddKm(
        MatrixType& rK,
        const Matrix& B,
        const Matrix& rConstitutiveMatrix,
        const double IntegrationWeight
    ) const; 

    /**
     * @brief
     */
    virtual void CalculateLinearMomentumResidualVector(
        VectorType& rRHSv,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const VectorType& rStressVector,
        const double weight
    );

    /**
     * @brief
    */
   virtual void CalculateGeometricalResidualVector(
        VectorType& rRHSF,
        KinematicVariables& rThisKinematicVariables,
        const double weight
   );

   /**
    * @brief This function calculates the K21 block of the LHS matrix
    * @details correspond to the derivative of the geometrical governing law for F with respect to the velocity DOFs
    */
   virtual void CalculateLeftHandSideK21Block(
        MatrixType& rK21,
        const KinematicVariables& rThisKinematicVariables,
        const double weight
    );

    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix The elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
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
    IntegrationMethod mThisIntegrationMethod;

    /**
     * @brief This function sets the used constitutive laws
     */
    void SetConstitutiveLaw(const ConstitutiveLaw::Pointer rThisConstitutiveLaw)
    {
        mThisConstitutiveLaw = rThisConstitutiveLaw;
    }

    /**
     * @brief Sets the used integration method
     * @details In SPH the inetgration points are the particles themselves,
     * this method is implement to maintain compatibility and display the results on the integration points 
    */
    void SetIntegrationMethod(const IntegrationMethod ThisIntegrationMethod)
    {
        mThisIntegrationMethod = ThisIntegrationMethod; ////
    }

    /**
     * @brief
     */
    virtual void AssembleLHS(
        MatrixType& rLHS,
        const MatrixType& rK11,
        const MatrixType& rK12,
        const MatrixType& rK21,
        const MatrixType& rK22
    );

    /**
     * @brief 
     */
    virtual void AssembleRHS(
        VectorType& rRHS,
        const VectorType& rRHSv,
        const VectorType& rRHSF
    );

    /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     */
    virtual void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues
        ) const;

    /**
     * @brief This functions updates the constitutive variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param ThisStressMeasure The stress measure considered
     */
    virtual void CalculateConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure = ConstitutiveLaw::StressMeasure_PK1
        ) const;
    
    /**
     * @brief Calculate the kinematics
     * @details This method calculates the kinematics of the element for a given integration point
     * @param rThisKinematicVariables Integration point kinematics container
     * @param PointNumber Integration point index
     * @param rIntegrationMethod Integration rule
     */
    virtual void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rProcessInfo
        );
    
    /**
     * @brief This function computes the deformation gradient of the particle
     * @details At the same time it stores the kernel and the kernel gradient values in the reference configuration
     * @param rDW_DX The matrix containing the kernel gradients in the reference configuration
     * @param rW The vector containing the kernels in the reference configuration
     */
    virtual void CalculateKernelsAndKernelGradients(
        MatrixType& rDW_DX,
        VectorType& rW,
        const ProcessInfo& rProcessInfo
    );

    /**
     * @brief
     */
    virtual void AssembleDeformationGradient(
        MatrixType& rF
    );

    /**
     * @brief It initializes the material 
     */
    void InitializeMaterial();

private:



};

}