//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "custom_elements/mpm_updated_lagrangian.hpp"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * Large Displacement Lagrangian Element for 3D and 2D geometries. (base class)
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class MPMUpdatedLagrangianAxisym
    : public MPMUpdatedLagrangian
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of LargeDisplacementElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMUpdatedLagrangianAxisym );
    ///@}

public:


    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    MPMUpdatedLagrangianAxisym();

    /// Default constructors
    MPMUpdatedLagrangianAxisym(
        IndexType NewId,
        GeometryType::Pointer pGeometry
    );

    MPMUpdatedLagrangianAxisym(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
    );

    ///Copy constructor
    MPMUpdatedLagrangianAxisym(MPMUpdatedLagrangianAxisym const& rOther);

    /// Destructor.
    ~MPMUpdatedLagrangianAxisym() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MPMUpdatedLagrangianAxisym& operator=(MPMUpdatedLagrangianAxisym const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes
    ) const override;

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize(
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MPM Axisymmetric Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MPM Axisymmetric Element #" << Id();
    }

protected:

    void InitializeGeneralVariables(
        GeneralVariables & rVariables,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    void CalculateKinematics(
        GeneralVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight
    ) override;

    /**
     * Calculation of the Deformation Matrix  BL
     */
    virtual void CalculateDeformationMatrix(
        Matrix& rB,
        const Matrix& rDN_DX,
        const Matrix& rN
    ) override;

    /// Calculation of the Deformation Gradient F
    void CalculateDeformationGradient(
        const Matrix& rDN_DX,
        Matrix& rF,
        Matrix& rDeltaPosition
    ) override;

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class MPMUpdatedLagrangian

} // namespace Kratos.
