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
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "includes/variables.h"

/* CONDITIONS */
// Mortar conditions
#include "custom_conditions/mesh_tying_mortar_condition.h"
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"
#include "custom_conditions/ALM_frictionless_components_mortar_contact_condition.h"
#include "custom_conditions/ALM_frictionless_mortar_contact_axisym_condition.h"
#include "custom_conditions/ALM_frictional_mortar_contact_condition.h"
#include "custom_conditions/ALM_frictional_mortar_contact_axisym_condition.h"

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

/// Short class definition.
/**
 * This application features Elements, Conditions, Constitutive laws and Utilities
 * for structural analysis problems
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) KratosContactStructuralMechanicsApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


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
    const MeshTyingMortarCondition<2, 3, ScalarValue> mMeshTyingMortarCondition2D2NTriangleScalar;            // 2DLine/Triangle for scalar variables
    const MeshTyingMortarCondition<2, 4, ScalarValue> mMeshTyingMortarCondition2D2NQuadrilateralScalar;       // 2DLine/Quadrilateral for scalar variables
    const MeshTyingMortarCondition<2, 3, Vector2DValue> mMeshTyingMortarCondition2D2NTriangleComponents;      // 2DLine/Triangle for components variables
    const MeshTyingMortarCondition<2, 4, Vector2DValue> mMeshTyingMortarCondition2D2NQuadrilateralComponents; // 2DLine/Quadrilateral for scalar variables
    const MeshTyingMortarCondition<3, 4, ScalarValue> mMeshTyingMortarCondition3D3NTetrahedronScalar;         // 3D Triangle/Tetrahedron for scalar variables
    const MeshTyingMortarCondition<3, 8, ScalarValue> mMeshTyingMortarCondition3D4NHexahedronScalar;          // 3D Quadrilateral/Hexahedra for scalar variables
    const MeshTyingMortarCondition<3, 4, Vector3DValue> mMeshTyingMortarCondition3D3NTetrahedronComponents;   // 3D Triangle/Tetrahedron for components variables
    const MeshTyingMortarCondition<3, 8, Vector3DValue> mMeshTyingMortarCondition3D4NHexahedronComponents;    // 3D Quadrilateral/Hexahedra for components variables
    
    // ALM Mortar contact conditions
    // Frictionless cases
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, false> mALMFrictionlessMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2, true> mALMNVFrictionlessMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<2, false> mALMFrictionlessAxisymMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<2, true> mALMNVFrictionlessAxisymMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, false> mALMFrictionlessMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3, true> mALMNVFrictionlessMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, false> mALMFrictionlessMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4, true> mALMNVFrictionlessMortarContactCondition3D4N;
    // Frictionless components cases
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<2, 2, false> mALMFrictionlessComponentsMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<2, 2, true> mALMNVFrictionlessComponentsMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 3, false> mALMFrictionlessComponentsMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 3, true> mALMNVFrictionlessComponentsMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 4, false> mALMFrictionlessComponentsMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<3, 4, true> mALMNVFrictionlessComponentsMortarContactCondition3D4N;
    // Frictional cases
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, false> mALMFrictionalMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, true> mALMNVFrictionalMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<2, false> mALMFrictionalAxisymMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<2, true> mALMNVFrictionalAxisymMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, false> mALMFrictionalMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, true> mALMNVFrictionalMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, false> mALMFrictionalMortarContactCondition3D4N;
    const AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, true> mALMNVFrictionalMortarContactCondition3D4N;
    
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


