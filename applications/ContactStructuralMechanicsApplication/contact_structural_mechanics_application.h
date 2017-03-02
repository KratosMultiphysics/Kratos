// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
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

// OLD Mortar conditions 
#include "custom_conditions/mortar_contact_condition.h"

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

#if !defined(TENSOR_VALUE)
#define TENSOR_VALUE
    enum TensorValue {ScalarValue = 1, Vector2DValue = 2, Vector2DPScalarValue = 3, Vector3DValue = 3, Vector3DPScalarValue = 4 };
#endif
    
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
class KratosContactStructuralMechanicsApplication : public KratosApplication
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
    virtual ~KratosContactStructuralMechanicsApplication() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register();



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
    virtual std::string Info() const
    {
        return "KratosContactStructuralMechanicsApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        KRATOS_WATCH("in my application");
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



    //       static const ApplicationCondition  msApplicationCondition;

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
    const MeshTyingMortarCondition<3, 6, ScalarValue> mMeshTyingMortarCondition3D4NHexahedronScalar;          // 3D Quadrilateral/Hexahedra for scalar variables
    const MeshTyingMortarCondition<3, 4, Vector3DValue> mMeshTyingMortarCondition3D3NTetrahedronComponents;   // 3D Triangle/Tetrahedron for components variables
    const MeshTyingMortarCondition<3, 6, Vector3DValue> mMeshTyingMortarCondition3D4NHexahedronComponents;    // 3D Quadrilateral/Hexahedra for components variables
    
    // ALM Mortar contact conditions
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, 2> mALMFrictionlessMortarContactCondition2D2N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 3> mALMFrictionlessMortarContactCondition3D3N;
    const AugmentedLagrangianMethodFrictionlessMortarContactCondition<3, 4> mALMFrictionlessMortarContactCondition3D4N;
    
    // OLD Mortar contact conditions
    const MortarContactCondition<2, 2> mMortarContactCondition2D2N;
    const MortarContactCondition<2, 3> mMortarContactCondition2D3N;
    const MortarContactCondition<3, 3> mMortarContactCondition3D3N;
    const MortarContactCondition<3, 4> mMortarContactCondition3D4N;
    
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


