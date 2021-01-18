// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

#if !defined(KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

/* GEOMETRIES */
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/line_2d_2.h"

/* CONDITIONS */
// Mortar conditions
#include "custom_conditions/mesh_tying_mortar_condition.h"
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"
#include "custom_conditions/ALM_frictionless_components_mortar_contact_condition.h"
#include "custom_conditions/penalty_frictionless_mortar_contact_condition.h"
#include "custom_conditions/ALM_frictionless_mortar_contact_axisym_condition.h"
#include "custom_conditions/penalty_frictionless_mortar_contact_axisym_condition.h"
#include "custom_conditions/ALM_frictional_mortar_contact_condition.h"
#include "custom_conditions/penalty_frictional_mortar_contact_condition.h"
#include "custom_conditions/ALM_frictional_mortar_contact_axisym_condition.h"
#include "custom_conditions/penalty_frictional_mortar_contact_axisym_condition.h"
#include "custom_conditions/mpc_mortar_contact_condition.h"

/* CONSTRAINTS */
#include "custom_master_slave_constraints/contact_master_slave_constraint.h"

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
 * @class KratosContactStructuralMechanicsApplication
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This application features Elements, Conditions, Constitutive laws and Utilities for structural analysis problems with contact constraints
 * @details It implements different types of conditions (penalty, ALM with mortar integration) and utilities (search, active set, convergence...)
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) KratosContactStructuralMechanicsApplication
    : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the node type
    typedef Node<3> NodeType;

    /// Definition of the geometry type
    typedef Geometry<NodeType> GeometryType;

    /// Definition of the geometry type
    typedef GeometryType::Pointer GeometryPointerType;

    /// Definition of the points array type
    typedef GeometryType::PointsArrayType PointsArrayType;

    /// Line type definition
    typedef Line2D2<NodeType> LineType;

    /// Triangle type definition
    typedef Triangle3D3<NodeType> TriangleType;

    /// Quadrilateral type definition
    typedef Quadrilateral3D4<NodeType> QuadrilateralType;

    /// Pointer definition of KratosContactStructuralMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosContactStructuralMechanicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosContactStructuralMechanicsApplication();

    /// Destructor.
    ~KratosContactStructuralMechanicsApplication() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosContactStructuralMechanicsApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    /* CONDITIONS*/
    // Mesh tying mortar condition
    const MeshTyingMortarCondition<2, 3> mMeshTyingMortarCondition2D2NTriangle;                   // 2DLine/Triangle
    const MeshTyingMortarCondition<2, 4> mMeshTyingMortarCondition2D2NQuadrilateral;              // 2DLine/Quadrilateral
    const MeshTyingMortarCondition<3, 4> mMeshTyingMortarCondition3D3NTetrahedron;                // 3D Triangle/Tetrahedron
    const MeshTyingMortarCondition<3, 8> mMeshTyingMortarCondition3D4NHexahedron;                 // 3D Quadrilateral/Hexahedra
    const MeshTyingMortarCondition<3, 4, 8> mMeshTyingMortarCondition3D3NTetrahedron4NHexahedron; // 3D Triangle/Tetrahedron-Quadrilateral/Hexahedra
    const MeshTyingMortarCondition<3, 8, 4> mMeshTyingMortarCondition3D4NHexahedron3NTetrahedron; // 3D Quadrilateral/Hexahedra-Triangle/Tetrahedron

    // ALM Mortar contact conditions
    // Frictionless cases
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, false> mALMFrictionlessMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, true> mALMNVFrictionlessMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<2, false> mALMFrictionlessAxisymMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<2, true> mALMNVFrictionlessAxisymMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, false, 3> mALMFrictionlessMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, true,  3> mALMNVFrictionlessMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, false, 4> mALMFrictionlessMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, true,  4> mALMNVFrictionlessMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, false, 4> mALMFrictionlessMortarContactCondition3D3N4N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, true,  4> mALMNVFrictionlessMortarContactCondition3D3N4N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, false, 3> mALMFrictionlessMortarContactCondition3D4N3N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, true,  3> mALMNVFrictionlessMortarContactCondition3D4N3N;
    // Frictionless components cases
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<2, 2, false> mALMFrictionlessComponentsMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<2, 2, true> mALMNVFrictionlessComponentsMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 3, false, 3> mALMFrictionlessComponentsMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 3, true,  3> mALMNVFrictionlessComponentsMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 3, false, 4> mALMFrictionlessComponentsMortarContactCondition3D3N4N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 3, true,  4> mALMNVFrictionlessComponentsMortarContactCondition3D3N4N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 4, false, 4> mALMFrictionlessComponentsMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 4, true,  4> mALMNVFrictionlessComponentsMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 4, false, 3> mALMFrictionlessComponentsMortarContactCondition3D4N3N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 4, true,  3> mALMNVFrictionlessComponentsMortarContactCondition3D4N3N;
    // Frictional cases
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, false> mALMFrictionalMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, true> mALMNVFrictionalMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<2, false> mALMFrictionalAxisymMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<2, true> mALMNVFrictionalAxisymMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, false, 3> mALMFrictionalMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, true,  3> mALMNVFrictionalMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, false, 4> mALMFrictionalMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, true,  4> mALMNVFrictionalMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, false, 4> mALMFrictionalMortarContactCondition3D3N4N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, true,  4> mALMNVFrictionalMortarContactCondition3D3N4N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, false, 3> mALMFrictionalMortarContactCondition3D4N3N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, true,  3> mALMNVFrictionalMortarContactCondition3D4N3N;
    // Frictionless penalty cases
    const PenaltyMethodFrictionlessMortarContactCondition<2, 2, false> mPenaltyFrictionlessMortarContactCondition2D2N;
    const PenaltyMethodFrictionlessMortarContactCondition<2, 2, true> mPenaltyNVFrictionlessMortarContactCondition2D2N;
    const PenaltyMethodFrictionlessMortarContactAxisymCondition<2, false> mPenaltyFrictionlessAxisymMortarContactCondition2D2N;
    const PenaltyMethodFrictionlessMortarContactAxisymCondition<2, true> mPenaltyNVFrictionlessAxisymMortarContactCondition2D2N;
    const PenaltyMethodFrictionlessMortarContactCondition<3, 3, false, 3> mPenaltyFrictionlessMortarContactCondition3D3N;
    const PenaltyMethodFrictionlessMortarContactCondition<3, 3, true,  3> mPenaltyNVFrictionlessMortarContactCondition3D3N;
    const PenaltyMethodFrictionlessMortarContactCondition<3, 4, false, 4> mPenaltyFrictionlessMortarContactCondition3D4N;
    const PenaltyMethodFrictionlessMortarContactCondition<3, 4, true,  4> mPenaltyNVFrictionlessMortarContactCondition3D4N;
    const PenaltyMethodFrictionlessMortarContactCondition<3, 3, false, 4> mPenaltyFrictionlessMortarContactCondition3D3N4N;
    const PenaltyMethodFrictionlessMortarContactCondition<3, 3, true,  4> mPenaltyNVFrictionlessMortarContactCondition3D3N4N;
    const PenaltyMethodFrictionlessMortarContactCondition<3, 4, false, 3> mPenaltyFrictionlessMortarContactCondition3D4N3N;
    const PenaltyMethodFrictionlessMortarContactCondition<3, 4, true,  3> mPenaltyNVFrictionlessMortarContactCondition3D4N3N;
    // Frictional penalty cases
    const PenaltyMethodFrictionalMortarContactCondition<2, 2, false> mPenaltyFrictionalMortarContactCondition2D2N;
    const PenaltyMethodFrictionalMortarContactCondition<2, 2, true> mPenaltyNVFrictionalMortarContactCondition2D2N;
    const PenaltyMethodFrictionalMortarContactAxisymCondition<2, false> mPenaltyFrictionalAxisymMortarContactCondition2D2N;
    const PenaltyMethodFrictionalMortarContactAxisymCondition<2, true> mPenaltyNVFrictionalAxisymMortarContactCondition2D2N;
    const PenaltyMethodFrictionalMortarContactCondition<3, 3, false, 3> mPenaltyFrictionalMortarContactCondition3D3N;
    const PenaltyMethodFrictionalMortarContactCondition<3, 3, true,  3> mPenaltyNVFrictionalMortarContactCondition3D3N;
    const PenaltyMethodFrictionalMortarContactCondition<3, 4, false, 4> mPenaltyFrictionalMortarContactCondition3D4N;
    const PenaltyMethodFrictionalMortarContactCondition<3, 4, true,  4> mPenaltyNVFrictionalMortarContactCondition3D4N;
    const PenaltyMethodFrictionalMortarContactCondition<3, 3, false, 4> mPenaltyFrictionalMortarContactCondition3D3N4N;
    const PenaltyMethodFrictionalMortarContactCondition<3, 3, true,  4> mPenaltyNVFrictionalMortarContactCondition3D3N4N;
    const PenaltyMethodFrictionalMortarContactCondition<3, 4, false, 3> mPenaltyFrictionalMortarContactCondition3D4N3N;
    const PenaltyMethodFrictionalMortarContactCondition<3, 4, true,  3> mPenaltyNVFrictionalMortarContactCondition3D4N3N;

    // MPC Conditions
    const MPCMortarContactCondition<2, 2> mMPCMortarContactCondition2D2N;
    const MPCMortarContactCondition<3, 3, 3> mMPCMortarContactCondition3D3N;
    const MPCMortarContactCondition<3, 4, 4> mMPCMortarContactCondition3D4N;
    const MPCMortarContactCondition<3, 3, 4> mMPCMortarContactCondition3D3N4N;
    const MPCMortarContactCondition<3, 4, 3> mMPCMortarContactCondition3D4N3N;

    /// Constraints
    const ContactMasterSlaveConstraint mContactMasterSlaveConstraint;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosContactStructuralMechanicsApplication& operator=(KratosContactStructuralMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosContactStructuralMechanicsApplication(KratosContactStructuralMechanicsApplication const& rOther);


    ///@}

}; // Class KratosContactStructuralMechanicsApplication

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED  defined


