//  ____  ____  _   _                   _ _           _   _             
// / ___||  _ \| | | | __ _ _ __  _ __ | (_) ___ __ _| |_(_) ___  _ __  
// \___ \| |_) | |_| |/ _` | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_ \ 
//  ___) |  __/|  _  | (_| | |_) | |_) | | | (_| (_| | |_| | (_) | | | |
// |____/|_|   |_| |_|\__,_| .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                         |_|   |_|                                    

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

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

template<class TKernelType, std::size_t TDim>
class KRATOS_API(SPH_APPLICATION) SmallDisplacementParticle : public Element
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

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementParticle);

    // Constructor void 
    SmallDisplacementParticle()
    {
    }

    // Constructor using an array of nodes 
    SmallDisplacementParticle(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    // Constructor using an array of nodes with properties 
    SmallDisplacementParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    { 
    }

    // Copy constructor
    SmallDisplacementParticle(SmallDisplacementParticle const& rOther)
        : BaseType(rOther),
        mThisConstitutiveLaw(rOther.mThisConstitutiveLaw)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SmallDisplacementParticle>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SmallDisplacementParticle>(NewId, pGeom, pProperties);
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
     * @brief Sets on rValues the nodal displacements
     */
    void GetValuesVector(VectorType& rValues, int step ) const override;

    /**
     * @brief Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(VectorType& rValues, int step = 0) const override;

    /**
     * @brief Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(VectorType& rValues, int step = 0) const override;
    
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
     * @brief This method returns if the element provides the strain
     */
    virtual bool UseElementProvidedStrain() const;

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
        const ProcessInfo& rProcessInfo,
        int Step = 0
    );

    /**
     * @brief This function calculates the external forces contribution
     * @param rBodyForce The Body force vector 
     */
    virtual void CalculateAndAddExternalForcesContribution(
        const VectorType& rW,
        const ProcessInfo& rProcessInfo,
        const VectorType& rBodyForce,
        VectorType& rRHS,
        const double weight
    ) const;

    /**
     * @brief This function is called to initialize the constitutive variables
     * @param ThisStressMeasure The stress measure to be used in the constitutive law
    */
    virtual void CalculateConstitutiveVariables(
        ConstitutiveVariables& rThisConstitutiveVariables,
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );

    /**
     * @brief This function is called to set the constitutive variables
     */
    virtual void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues
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

    /**
     * @brief This method gets a value directly in the CL avoiding code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained in the integration points
     * @tparam TType The type considered
     */

    template<class TType>
    void GetValueOnConstituitiveLaw(const Variable<TType>& rVariable, std::vector<TType>& rOutput){
        mThisConstitutiveLaw->GetValue(rVariable, rOutput[0]); 
    }

    /**
     * @brief This method computes directly in the CL
     * @details Avoids code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained in the integration points
     * @tparam TType The type considered
     */
    template<class TType>
    void CalculateOnConstitutiveLaw(
        const Variable<TType>& rVariable,
        std::vector<TType>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
        {

            const auto& r_geom = GetGeometry();
            const auto& r_props = GetProperties();
            const SizeType dimension = r_geom.WorkingSpaceDimension();
            const auto& r_neighbours = this->GetValue(NEIGHBOURS);

            rOutput.resize(1);
            const SizeType strain_size = mThisConstitutiveLaw->GetStrainSize();
            
            KinematicVariables this_kinematic_variables(strain_size, dimension, r_neighbours.size());
            ConstitutiveVariables this_constitutive_variables(strain_size);

            ConstitutiveLaw::Parameters cl_values(r_geom, r_props, rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& r_cl_options = cl_values.GetOptions();
            r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            cl_values.SetStrainVector(this_constitutive_variables.StrainVector);
            
            CalculateKinematicVariables(this_kinematic_variables, rCurrentProcessInfo);
            SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values);

            rOutput[0] = mThisConstitutiveLaw->CalculateValue(cl_values, rVariable, rOutput[0]);
        }
};

}