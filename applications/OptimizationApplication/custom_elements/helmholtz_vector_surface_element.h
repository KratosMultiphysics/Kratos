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
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "geometries/pyramid_3d_5.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "optimization_application_variables.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/**
 * @class HelmholtzVectorSurfaceElement
 * @ingroup OptimizationApplication
 * @brief Helmholtz filtering element for 3D geometries.
 * @details Implements Helmholtz surface PDEs which involve Laplace-Beltrami operations for filtering a vector field on surface.
 * @author Reza Najian Asl
 */
class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzVectorSurfaceElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{
        typedef Node                                PointType;
        typedef Node::Pointer                       PointPtrType;
        typedef Geometry<PointType>                 GeometryType;
        typedef Pyramid3D5<PointType>               PyramidGeometryType;
        typedef Tetrahedra3D4<PointType>            TetrahedraGeometryType;
    /// Counted pointer of HelmholtzVectorSurfaceElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(HelmholtzVectorSurfaceElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HelmholtzVectorSurfaceElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);
    HelmholtzVectorSurfaceElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~HelmholtzVectorSurfaceElement();

    ///@}
    ///@name Operations
    ///@{

   /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

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
     * @brief Sets on rConditionalDofList the degrees of freedom of the considered element geometry
     * @param rConditionalDofList The vector containing the dof of the element
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
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
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
      * @brief This is called during filtering to compute element strain energy
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
    HelmholtzVectorSurfaceElement() : Element()
    {
    }

private:
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}

    ///@name Private Operations
    ///@{
    /**
      * @brief This is called during the assembling process in order to calculate stiffness matrix
      * @param rStiffnessMatrix The stiffness matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateStiffnessMatrix(
        MatrixType& rStiffnessMatrix,
        const ProcessInfo& rCurrentProcessInfo) const;

    /**
      * @brief This is called during the assembling process in order to calculate shape functions
      * from an pseudo bulk element
      * @param rNMatrix The shape functions matrix
      * @param rIntegrationMethod The given integration method
      * @param rCurrentProcessInfo The current process info instance
      */
    void GetPseudoBulkSurfaceShapeFunctionsValues(
        MatrixType& rNMatrix,
        const IntegrationMethod& rIntegrationMethod,
        const ProcessInfo& rCurrentProcessInfo) const;

    /**
      * @brief This is called during the assembling process in order to calculate shape functions
      * gradients from an pseudo bulk element at the given gp point index
      * @param rDN_DX The shape functions gradient matrix
      * @param rIntegrationMethod The given integration method.
      * @param PointNumber The index of the desired gp.
      * @param rCurrentProcessInfo The current process info instance.
      */
    void CalculatePseudoBulkSurfaceDN_DXMatrix(
        MatrixType& rDN_DX,
        const IntegrationMethod& rIntegrationMethod,
        const IndexType PointNumber,
        const ProcessInfo& rCurrentProcessInfo) const;

    /**
      * @brief This is called during the assembling process in order to element surface normal
      * @param rNormal The averaged surface normal
      */
    void CalculateAvgSurfUnitNormal(
        VectorType& rNormal) const;
    ///@}

}; // Class HelmholtzVectorSurfaceElement

///@}

}  // namespace Kratos.



