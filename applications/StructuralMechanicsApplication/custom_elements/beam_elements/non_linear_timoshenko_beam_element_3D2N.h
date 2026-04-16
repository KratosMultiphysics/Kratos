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

#include "includes/element.h"

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
 * @brief This is the 3D Timoshenko beam element of 2 nodes. 
 * Reference:  I. Romero and F. Armero, "An objective finite element approximation o the kinetics of geometrically exact rods
 * and its use in the formulation of an energy-momentum conservative scheme in dynamics", IJNME, 2002.
 * DOI: 10.1002/nme.486 
 * @author Alejandro Cornejo
 */
class NonLinearTimoshenkoBeamElement3D2N : public Element
{
public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NonLinearTimoshenkoBeamElement3D2N);

    ///@name Type Definitions
    ///@{
    using BaseType = Element;

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

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief It initializes the material
     */
    void InitializeMaterial();

        // void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rProcessInfo) override;

    /**
     * @brief Indicates the amount of DoFs per node (u, v, w, theta_x, theta_y, theta_z)
     */
    IndexType GetDoFsPerNode() const
    {
        return 6; // 3 displacements and 3 rotations
    }

private:

    /* The rotation operators are built with [d1, d2, d3] as col vectors */
    std::vector<BoundedMatrix<double, 3, 3>> mRotationOperators; // the two rotation matrices, one per each IP.
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; // The vector containing the beam constitutive laws, one per each IP
    IntegrationMethod mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_LOBATTO_1; // By default the quadrature points are located at the nodes of the beam
};

} // namespace Kratos
