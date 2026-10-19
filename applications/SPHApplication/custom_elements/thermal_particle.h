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

namespace Kratos
{
using SizeType = std::size_t;

template<class TKernelType, std::size_t TDim>
class KRATOS_API(SPH_APPLICATION) ThermalParticle : public Element
{
public: 

    using BaseType = Element;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThermalParticle);

    // Constructor void 
    ThermalParticle()
    {
    }

    // Constructor using an array of nodes 
    ThermalParticle(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    // Constructor using an array of nodes with properties 
    ThermalParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    { 
    }

    // Copy constructor
    ThermalParticle(ThermalParticle const& rOther)
        : BaseType(rOther)
    {
    }

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ThermalParticle>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ThermalParticle>(NewId, pGeom, pProperties);
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override;

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
     * @brief Sets on rValues the support temperatures
     */
    void GetValuesVector(VectorType& rValues, int step = 0) const override;

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
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix The elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;
    
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rProcessInfo) override;
    
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

protected: 
        /**
         * @brief Internal variables used in the spatial calculations
         */
        struct SpatialVariables
        {
            VectorType W; 
            MatrixType DW_DX;

            /**
             * @brief Default constructor
             * @param DomainSize The size of the problem domain
             * @param NumberOfNeighbours The number of neighbours of the particle
             */
            SpatialVariables(SizeType DomainSize, SizeType NumberOfNeighbours)
            {
                W = ZeroVector(NumberOfNeighbours);
                DW_DX = ZeroMatrix(NumberOfNeighbours, DomainSize);
            }
        };

        /**
         * @brief
         */
        virtual void CalculateSpatialVariables(
            SpatialVariables& rThisSpatialVariables, 
            const ProcessInfo& rProcessInfo
        );

private: 
    

};

}
