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

#include "timoshenko_beam_element_3D2N.h"
#include "utilities/quaternion.h"

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
class NonLinearTimoshenkoBeamElement3D2N : public LinearTimoshenkoBeamElement3D2N
{
public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NonLinearTimoshenkoBeamElement3D2N);

    ///@name Type Definitions
    ///@{
    using BaseType = LinearTimoshenkoBeamElement3D2N;

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

        // void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rProcessInfo) override;

private:
    /**
     * @brief This quaternion is located at the centre of the beam and
     * it will be updated at any time that the rotations are updated.
     * It is used to compute the rotation matrix that transforms the local internal forces to the global frame, and also to compute the local deformation of the beam.
     * It is initialized as the identity quaternion, and then it is updated at each iteration using the incremental rotations of the nodes.
     * The use of a quaternion allows to avoid singularities and to have a more robust update of the rotation matrix, especially for large rotations.
     * The use of a quaternion also allows to have a more efficient update of the rotation matrix.
     */
    Quaternion<double> mQuaternion = Quaternion<double>::Identity();


};

} // namespace Kratos
