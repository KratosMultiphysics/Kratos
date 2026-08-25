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

#include "custom_elements/total_lagrangian_particle.h"

namespace Kratos
{
    
using SizeType = std::size_t;

template<class TKernelType, std::size_t TDim>
class KRATOS_API(SPH_APPLICATION) TotalLagrangianMixedStrainParticle : public TotalLagrangianDisplacementParticle<TKernelType, TDim>
{

public:

    using BaseType = TotalLagrangianDisplacementParticle<TKernelType, TDim>;

    using GeometryType = typename BaseType::GeometryType;
    using PropertiesType = typename BaseType::PropertiesType;
    using IndexType = typename BaseType::IndexType;
    using MatrixType = typename BaseType::MatrixType;
    using VectorType = typename BaseType::VectorType;
    using NodesArrayType = typename BaseType::NodesArrayType;
    using EquationIdVectorType = typename BaseType::EquationIdVectorType;
    using DofsVectorType = typename BaseType::DofsVectorType;
    
    using KinematicVariables = typename BaseType::KinematicVariables;
    using ConstitutiveVariables = typename BaseType::ConstitutiveVariables;

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
    }

    // Copy constructor
    TotalLagrangianMixedStrainParticle(TotalLagrangianMixedStrainParticle const& rOther)
        : BaseType(rOther)
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
     * @brief Sets on rValues the nodal velocities and the nodal deformation gradient components
     */
    void GetValuesVector(Vector& rValues, int step = 0) const override; 
    
    /**
     * @brief Sets on rValues the nodal acellerations and the nodal deformation gradient rate components
     */
    void GetFirstDerivativesVector(VectorType& rValues, int step = 0) const override;

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
      * @brief Calculation of the Geometric Stiffness Matrix component.
      * @param StressVector The vector containing the stress components
      */
    void CalculateAndAddKg(
        MatrixType& rK12,
        const MatrixType& rDW_DX,
        const VectorType& rStressVector,
        const double IntegrationWeight
    ) const;

    /**
      * @brief Calculation of the Material Stiffness Matrix component.
      * @param StressVector The vector containing the stress components
      */
    void CalculateAndAddKm(
        MatrixType& rK12,
        const KinematicVariables& rThisKinematicVariables,
        const MatrixType& rConstitutiveMatrix,
        const double IntegrationWeight
    ) const;

    /**
     * @brief This function calculates the linear momentum residual vector
     */
    virtual void CalculateLinearMomentumResidualVector(
        VectorType& rRHSv,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const VectorType& rStressVector,
        const double weight,
        const int Step = 0
    );

   /**
    * @brief This function calculates the K21 block of the LHS matrix
    * @details correspond to the derivative of the geometrical governing law for F with respect to the velocity DOFs
    */
   virtual void CalculateGeometricalTangentMatrix(
        MatrixType& rK21,
        const KinematicVariables& rThisKinematicVariables,
        const double weight
    );

    /**
     * @brief This function calculates the geometrical residual vector
     * @details correspond to the governing law for F
    */
   virtual void CalculateGeometricalResidualVector(
        VectorType& rRHSF,
        KinematicVariables& rThisKinematicVariables,
        const double weight,
        const int Step = 0
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
     * @brief This function calculates the upwind stabilization residual vector
     */
    virtual void CalculateAndAddUpwindStabilizationResidual(
        VectorType& rRHSv,
        KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rProcessInfo,
        int Step = 0
    );

    /**
     * @brief This function calculates the upwind stabilization tangent matrix
     */
    virtual void CalculateAndAddUpwindStabilizationTangent(
        MatrixType& rK11,
        KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rProcessInfo
    );

protected:

    /**
     * @brief This function assembles the LHS matrix of the element block by block 
     */
    virtual void AssembleLHS(
        MatrixType& rLHS,
        const MatrixType& rK11,
        const MatrixType& rK12,
        const MatrixType& rK21,
        const MatrixType& rK22
    );

    /**
     * @brief This function assembles the RHS matrix of the element block by block 
     */
    virtual void AssembleRHS(
        VectorType& rRHS,
        const VectorType& rRHSv,
        const VectorType& rRHSF
    );

    /**
     * @brief Calculate the kinematics
     * @details This method calculates the kinematics of the element for a given integration point
     * @param rThisKinematicVariables Integration point kinematics container
     */
    virtual void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rProcessInfo,
        const int Step = 0
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
     * @brief Assemble the deformation gradient of the particle using directly the dofs values of particle itself
     */
    virtual void AssembleDeformationGradient(
        MatrixType& rF,
        const int Step = 0
    );

private:


};

}