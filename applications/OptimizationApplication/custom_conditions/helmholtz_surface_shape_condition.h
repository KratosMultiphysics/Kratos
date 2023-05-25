//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "includes/model_part.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "optimization_application_variables.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/**
 * @class HelmholtzSurfaceShapeCondition
 * @ingroup OptimizationApplication
 * @brief Helmholtz filtering condition for 3D geometries.
 * @details Implements Laplace–Beltrami operator on the surface boundaries of a bulk/solid domains. Also note that
 * this condition should be used with helmholtz_solid_shape_element
 * @author Reza Najian Asl
 */
class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzSurfaceShapeCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of HelmholtzSurfaceShapeCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(HelmholtzSurfaceShapeCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HelmholtzSurfaceShapeCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry);

    HelmholtzSurfaceShapeCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~HelmholtzSurfaceShapeCondition();

    ///@}
    ///@name Operations
    ///@{

   /**
     * @brief Creates a new condition
     * @param NewId The Id of the new created condition
     * @param pGeom The pointer to the geometry of the condition
     * @param pProperties The pointer to property
     * @return The pointer to the created condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition
     * @param NewId The Id of the new created condition
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief It creates a new condition pointer and clones the previous condition data
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief Sets on rResult the ID's of the condition degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on rConditionalDofList the degrees of freedom of the considered condition geometry
     * @param rConditionalDofList The vector containing the dof of the condition
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rConditionalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on rValues the nodal values
     * @param rValues The values of values
     * @param Step The step to be computed
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * @brief This function provides a more general interface to the condition.
     * @details It is designed so that rLHSvariables and rRHSvariables are passed to the condition thus telling what is the desired output
     * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
     * @param rRightHandSideVector container for the desired RHS output
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the conditional left hand side matrix only
     * @param rLeftHandSideMatrix the conditional left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the conditional right hand side vector only
      * @param rRightHandSideVector the conditional right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;


    /**
      * @brief This is called during the assembling process in order to calculate the conditional mass matrix
      * @param rMassMatrix The conditional mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during filtering to compute condition strain energy
      * @param rOutput The calculated variable
      * @param rCurrentProcessInfo The current process info instance
      */
    void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) const override;

    ///@}

protected:
    // Protected default constructor necessary for serialization
    HelmholtzSurfaceShapeCondition() : Condition()
    {
    }

private:
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}

    ///@name Private Operations
    ///@{
    /**
      * @brief This is called during the assembling process in order to calculate the condition stiffness matrix
      * @param rStiffnessMatrix The condition stiffness matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateStiffnessMatrix(
        MatrixType& rStiffnessMatrix,
        const ProcessInfo& rCurrentProcessInfo) const;

    /**
      * @brief This is called during the assembling process in order to calculate surface normal
      * @param rN The surface normal
      */
    void CalculateNormal(
        VectorType & rN) const;

    /**
      * @brief This is called during the assembling process in order to calculate shape function using parent solid element
      * @param rNMatrix The matrix of shape functions
      * @param rIntegrationMethod The given integrations method
      * @param rCurrentProcessInfo The current process info instance
      */
    void GetParentElementShapeFunctionsValues(
        MatrixType& rNMatrix,
        const IntegrationMethod& rIntegrationMethod,
        const ProcessInfo& rCurrentProcessInfo) const;

    /**
      * @brief This is called during the assembling process in order to calculate shape function gradients from parent solid element
      * @param rDN_DX The matrix of shape functions gradients
      * @param rIntegrationMethod The given integrations method
      * @param rCurrentProcessInfo The current process info instance
      */
    void GetParentElementShapeFunctionsGlobalGradients(
        MatrixType& rDN_DX,
        const IndexType PointNumber,
        const IntegrationMethod& rIntegrationMethod,
        const ProcessInfo& rCurrentProcessInfo) const;
    ///@}

}; // Class HelmholtzSurfaceShapeCondition

///@}

}  // namespace Kratos.



