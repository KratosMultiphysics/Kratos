//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti, Ricky Aristio
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "iga_application_variables.h"
#include "base_3d_beam_element.h"

namespace Kratos{
///@name Kratos Classes
///@{
/// Short class definition.
/** 3D Bernoulli beam element.
*/
class KRATOS_API(IGA_APPLICATION) Bernoulli3DBeamElement
    : public Base3DBeamElement
{
protected:
    struct KinematicVariables
    {
        // Reference configuration
        array_1d<double, 3> R1;     // Reference tangent vector
        array_1d<double, 3> R2;     // Reference curvature vector  
        array_1d<double, 3> R3;     // Reference third derivative
        double A;                    // Reference length measure
        double B;                    // Reference curvature measure
        
        // Current configuration
        array_1d<double, 3> r1;     // Current tangent vector
        array_1d<double, 3> r2;     // Current curvature vector
        array_1d<double, 3> r3;     // Current third derivative
        double a;                    // Current length measure
        double b;                    // Current curvature measure
        
        // Cross-section directors
        array_1d<double, 3> N0, V0;  // Reference directors
        array_1d<double, 3> n, v;    // Current directors
        
        // Curvature components
        double B_n, B_v;             // Reference curvatures
        double b_n, b_v;             // Current curvatures
        double C_12, C_13;           // Reference twist
        double c_12, c_13;           // Current twist
        
        // Rotations
        double Phi, Phi_der;    // Reference rotation
        double phi, phi_der;    // Current rotation
        
        KinematicVariables(SizeType Dimension = 3)
        {
            R1 = ZeroVector(Dimension);
            R2 = ZeroVector(Dimension);
            R3 = ZeroVector(Dimension);
            r1 = ZeroVector(Dimension);
            r2 = ZeroVector(Dimension);
            r3 = ZeroVector(Dimension);
            N0 = ZeroVector(Dimension);
            V0 = ZeroVector(Dimension);
            n = ZeroVector(Dimension);
            v = ZeroVector(Dimension);
            A = B = a = b = 0.0;
            B_n = B_v = b_n = b_v = 0.0;
            C_12 = C_13 = c_12 = c_13 = 0.0;
            Phi = Phi_der = 0.0;
            phi = phi_der =  0.0;
        }
    };

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Bernoulli3DBeamElement);

    using BaseType = Base3DBeamElement;
    using SizeType = typename BaseType::SizeType;
    using IndexType = typename BaseType::IndexType;

    // GeometryType
    using GeometryType = Geometry<Node>;

    //@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    Bernoulli3DBeamElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Base3DBeamElement(NewId, pGeometry)
    {
    };

    /// Constructor using an array of nodes with properties
    Bernoulli3DBeamElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Base3DBeamElement(NewId, pGeometry, pProperties)
    {
    };

    /// Default constructor necessary for serialization
    Bernoulli3DBeamElement()
        : Base3DBeamElement()
    {
    };

    /// Destructor.
    ~Bernoulli3DBeamElement() override = default;

     /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<Bernoulli3DBeamElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<Bernoulli3DBeamElement>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Degrees of freedom
    ///@{

    /// Relates the degrees of freedom of the element geometry
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /// Sets the ID's of the element degrees of freedom
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    void ComputeGAxial(
        const IndexType IntegrationPointIndex,
        Matrix& rGAxial,
        KinematicVariables& rKinematicVariables) const;

    void ComputeGBending(
        const IndexType IntegrationPointIndex,
        Matrix& rGBending1,
        Matrix& rGBending2,
        KinematicVariables& rKinematicVariables) const;
    
    void ComputeGTorsion(
        const IndexType IntegrationPointIndex,
        Matrix& rGTorsion1,
        Matrix& rGTorsion2,
        KinematicVariables& rKinematicVariables) const;

    void ComputeBAxial(
        const IndexType IntegrationPointIndex,
        Matrix& rBAxial,
        KinematicVariables& rKinematicVariables) const;

    void ComputeBBending(
        const IndexType IntegrationPointIndex,
        Matrix& rBBending1,
        Matrix& rBBending2,
        KinematicVariables& rKinematicVariables) const;
    
    void ComputeBTorsion(
        const IndexType IntegrationPointIndex,
        Matrix& rBTorsion1,
        Matrix& rBTorsion2,
        KinematicVariables& rKinematicVariables) const;

private:

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = this->GetGeometry().size();
        const SizeType number_dofs_per_node = 4;
        const SizeType nb_dofs = nb_nodes * number_dofs_per_node;

        if (rRightHandSideVector.size() != nb_dofs)
            rRightHandSideVector.resize(nb_dofs);
        noalias(rRightHandSideVector) = ZeroVector(nb_dofs);

        MatrixType left_hand_side_matrix;

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = this->GetGeometry().size();
        const SizeType number_dofs_per_node = 4;
        const SizeType nb_dofs = nb_nodes * number_dofs_per_node;

        VectorType right_hand_side_vector;
         
        if (rLeftHandSideMatrix.size1() != nb_dofs && rLeftHandSideMatrix.size2() != nb_dofs)
            rLeftHandSideMatrix.resize(nb_dofs, nb_dofs);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(nb_dofs, nb_dofs);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = GetGeometry().size();
        const SizeType number_dofs_per_node = 4;
        const SizeType nb_dofs = nb_nodes * number_dofs_per_node;

        if (rRightHandSideVector.size() != nb_dofs)
            rRightHandSideVector.resize(nb_dofs);
        noalias(rRightHandSideVector) = ZeroVector(nb_dofs);

        if (rLeftHandSideMatrix.size1() != nb_dofs && rLeftHandSideMatrix.size2() != nb_dofs)
            rLeftHandSideMatrix.resize(nb_dofs, nb_dofs);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(nb_dofs, nb_dofs);

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    ) override
    {}

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
    }

    void load(Serializer& rSerializer) override
    {
    }
};
}