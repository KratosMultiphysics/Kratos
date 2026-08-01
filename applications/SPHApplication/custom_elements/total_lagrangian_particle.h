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

#include "custom_elements/small_displacement_particle.h"

namespace Kratos
{

using SizeType = std::size_t;

template<class TKernelType>
class KRATOS_API(SPH_APPLICATION) TotalLagrangianDisplacementParticle: public SmallDisplacementParticle<TKernelType>
{
public: 

    using BaseType = SmallDisplacementParticle<TKernelType>;

    using GeometryType = typename BaseType::GeometryType;
    using PropertiesType = typename BaseType::PropertiesType;
    using IndexType = typename BaseType::IndexType;
    using MatrixType = typename BaseType::MatrixType;
    using VectorType = typename BaseType::VectorType;
    using SizeType = typename BaseType::SizeType;
    using NodesArrayType = typename BaseType::NodesArrayType;
    
    using KinematicVariables = typename BaseType::KinematicVariables;
    using ConstitutiveVariables = typename BaseType::ConstitutiveVariables;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

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
        : BaseType(rOther)
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
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function tells the position of the particle in the list of neighbours
     */
    int GetNeighbourPosition(const std::vector<Element::Pointer>& rNeighbours) const
    {
        int i = 0; 
        
        while (i<rNeighbours.size() && this->Id() != rNeighbours[i]->Id()) i++;

        return i;
    }

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
     * @details At the same time it strores the kernel and the kernel gradient values in the reference configuration
     * @param rDW_DX The matrix containing the kernel gradients in the reference configuration
     * @param rW The vector containing the kernels in the reference configuration
     */
    virtual void CalculateDeformationGradient(
        MatrixType& rF,
        MatrixType& rDW_DX,
        VectorType& rW,
        const ProcessInfo& rProcessInfo
    );

    /**
     * @brief This function computes the deformation matrix B for 2D simulations
     */
    void Calculate2DB(
        MatrixType& rB,
        const MatrixType& rF,
        const MatrixType& rDW_DX
    );

    /**
     * @brief This function computes the deformation matrix B for 3D simulations 
     */
    void Calculate3DB(
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
      * @param StressVector The vector containing the stress components
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
      * @param B The deformation matrix (Total Lagrangian Framework)
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
      */
    void CalculateAndAddResidualVector(
        VectorType& rRightHandSideVector,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const VectorType& rBodyForce,
        const VectorType& rStressVector,
        const double IntegrationWeight
    ) const; 

    /**
     * @brief This function calculates the contribution of the penalization to the LHS and RHS
     * @details A pairwise elastic spring is introduced between neighboring particles to suppress spurious zero-energy 
     * modes and improve the stability of the discretization. 
     * @cite "A review of mesh-free Smoothed Particle Hydrodynamics for large strain solid dynamics: 
     * from displacement-based formulations to first-order conservation laws"  Chun Hean Lee et al. 2026
     */
    virtual void CalculateAndAddPenalization(
        MatrixType& rLHS,
        VectorType& rRHS,
        KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
    ); 

    /**
      * @brief This is called during the assembling process in order to calculate the elemental damping matrix
      * @param rDampingMatrix The elemental damping matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;
    
    /**
     * @brief This function calculates the contribution of the dissipation to the LHS and RHS
     * @details A viscous interaction is introduced between neighboring particles to damp high-frequency oscillations, reduce spurious zero-energy modes, 
     * and improve the numerical stability during shock-wave propagation and dynamic simulations.
     * @cite "A review of mesh-free Smoothed Particle Hydrodynamics for large strain solid dynamics: 
     * from displacement-based formulations to first-order conservation laws"  Chun Hean Lee et al. 2026
     */
    virtual void CalculateAndAddDissipation(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rProcessInfo
    );

    /**
     * @brief These functions calculates the values of variables in the integrations points.
     * In SPH case coincide with the neighbouring particles 
     * @details These functions expect a std::vector of values for the specified variable type
     * @param rVariable This parameter selects the output 
     * @param F_DEFORMATION_GRADIENT The function computes the deformation gradient in the neighbours  
     */

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

};

}