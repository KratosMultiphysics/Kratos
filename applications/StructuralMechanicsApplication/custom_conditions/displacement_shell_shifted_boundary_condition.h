// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Ricky Aristio   
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/condition.h"
#include "geometries/geometry.h"
#include "includes/variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DisplacementShellShiftedBoundaryCondition: public Condition
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of DisplacementShellShiftedBoundaryCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(DisplacementShellShiftedBoundaryCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with geometry
    DisplacementShellShiftedBoundaryCondition(
        IndexType NewId,
        Geometry<Node>::Pointer pGeometry);

    /// Constructor with geometry and properties
    DisplacementShellShiftedBoundaryCondition(
        IndexType NewId,
        Geometry<Node>::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Constructor with array of nodes and properties
    DisplacementShellShiftedBoundaryCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes);

    /// Copy constructor
    DisplacementShellShiftedBoundaryCondition(DisplacementShellShiftedBoundaryCondition const& rOther) = delete;

    /// Destructor
    ~DisplacementShellShiftedBoundaryCondition() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    DisplacementShellShiftedBoundaryCondition& operator=(DisplacementShellShiftedBoundaryCondition const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& ConditionalDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
protected:
    ///@name Protected Life Cycle
    ///@{

    // Internal default constructor for serialization
    DisplacementShellShiftedBoundaryCondition();

    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Computes the normal projection of the B transpose times C product
     * This function calculates the normal projection of the standard
     * B transpose (strain matrix tranpose) times C (constitutive tensor) product
     * @param rC Reference to the constituive tensor
     * @param rB Reference to the strain matrix
     * @param rUnitNormal Reference to the unit normal vector
     * @param rAuxMat Output result
     */
    void CalculateBtransCProjectionLinearisation(
        const Matrix& rC,
        const Matrix& rB,
        const array_1d<double, 3>& rUnitNormal,
        Matrix& rAuxMat);

    /**
     * @brief Auxiliary shape functions calculation
     * This function expands the shape functions vector in a matrix
     * according to the problem dimension
     * @param rN Vector containing the shape function values
     * @param rConstrainedDofs Per-DOF mask (ux,uy,uz,rx,ry,rz): 1.0 = row active, 0.0 = row left zero
     * @param rAuxMat Result shape functions matrix
     */
    void CalculateAuxShapeFunctionsMatrix(
        const Vector& rN,
        const array_1d<double, 6>& rConstrainedDofs,
        Matrix& rAuxMat);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
}; // Class DisplacementShellShiftedBoundaryCondition

///@}

///@} addtogroup block

}  // namespace Kratos.
