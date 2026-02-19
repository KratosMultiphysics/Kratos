// KRATOS  ___|  |                   |                   |                   
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |             
//             | |   |    |   | (    |   |   | |   (   | |             
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
#pragma once

#include "non_linear_timoshenko_beam_element_2D2N.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class NonLinearTimoshenkoBeamElement3D2N
 * @ingroup StructuralMechanicsApplication
 * @brief This is the 3D Timoshenko beam element of 2 nodes. It is formulated in a Total Lagrangian fashion.
 * This is not a corotational approach. The non-linearity is introduced by the use of the exact beam kinematics, which results in a non-linear strain-displacement relationship.
 * Reference: 
 * @author Alejandro Cornejo
 */
class NonLinearTimoshenkoBeamElement3D2N : public NonLinearTimoshenkoBeamElement2D2N
{
public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NonLinearTimoshenkoBeamElement3D2N);

    struct BeamElementData
    {
        double du; // displacement derivative in x
        double dv; // displacement derivative in y
        double dw; // displacement derivative in z

        double theta_x; // rotation about x
        double theta_y; // rotation about y
        double theta_z; // rotation about z

        double dtheta_x; // Derivative of rotation about x
        double dtheta_y; // Derivative of rotation about y
        double dtheta_z; // Derivative of rotation about z

        // Displacements shape functions
        Vector Nu;
        Vector Nv;
        Vector Nw;

        // Derivatives of displacements shape functions
        Vector dNu;
        Vector dNv;
        Vector dNw;

        // Rotations shape functions
        Vector Ntheta_x;
        Vector Ntheta_y;
        Vector Ntheta_z;

        // Derivatives of rotations shape functions
        Vector dNtheta_x;
        Vector dNtheta_y;
        Vector dNtheta_z;

        Matrix B; // Strain-displacement matrix

        double cx, cy, cz; // Cosines of the rotations
        double sx, sy, sz; // Sines of the rotations
        double L; // Reference length of the beam element

        void Initialize()
        {
            du = dv = dw = 0.0;
            theta_x = theta_y = theta_z = 0.0;
            dtheta_x = dtheta_y = dtheta_z = 0.0;

            Nu.resize(12);
            Nv.resize(12);
            Nw.resize(12);
            Nu.clear();
            Nv.clear();
            Nw.clear();

            dNu.resize(12);
            dNv.resize(12);
            dNw.resize(12);
            dNu.clear();
            dNv.clear();
            dNw.clear();

            Ntheta_x.resize(12);
            Ntheta_y.resize(12);
            Ntheta_z.resize(12);
            Ntheta_x.clear();
            Ntheta_y.clear();
            Ntheta_z.clear();

            dNtheta_x.resize(12);
            dNtheta_y.resize(12);
            dNtheta_z.resize(12);
            dNtheta_x.clear();
            dNtheta_y.clear();
            dNtheta_z.clear();

            B.resize(6, 12);
            B.clear();

            cx = cy = cz = 0.0;
            sx = sy = sz = 0.0;
            L = 0.0;
        }
    };

    ///@name Type Definitions
    ///@{
    using BaseType = NonLinearTimoshenkoBeamElement2D2N;

    NonLinearTimoshenkoBeamElement3D2N() {}

    NonLinearTimoshenkoBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {}

    NonLinearTimoshenkoBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {}

    NonLinearTimoshenkoBeamElement3D2N(NonLinearTimoshenkoBeamElement3D2N const& rOther)
        : BaseType(rOther)
    {}

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearTimoshenkoBeamElement3D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearTimoshenkoBeamElement3D2N>(NewId, pGeom, pProperties);
    }

    Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override;


    /**
     * @brief Calculate local system
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate left hand side
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate right hand side
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief General calculation method with flags
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLHS,
        const bool ComputeRHS
        );


    // void BuildShapeFunctionsValuesVectors(
    //     const double xi,
    //     const double du,
    //     const double dv,
    //     const double dw,
    //     const double theta_x,
    //     const double theta_y,
    //     const double theta_z,
    //     Vector& rdNu,
    //     Vector& rdNv,
    //     Vector& rdNw,
    // )

    // void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rProcessInfo) override;
};

} // namespace Kratos
