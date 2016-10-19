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

// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

namespace Kratos
{
KratosContactStructuralMechanicsApplication::KratosContactStructuralMechanicsApplication():
    /* CONDITIONS */
    // Contact mortar conditions
    mMortarContactCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mMortarContactCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mMortarContactCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mMortarContactCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6 ) ) ) ),
    mMortarContactCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
    mMortarContactCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8 ) ) ) ),
    mMortarContactCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9 ) ) ) )
{}

void KratosContactStructuralMechanicsApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    // CONDITIONS
    /* Mortar method */
    KRATOS_REGISTER_VARIABLE( CONTACT_CONTAINERS )                          // A vector of which contains the structure which defines the contact conditions
    KRATOS_REGISTER_VARIABLE( INTEGRATION_ORDER_CONTACT )                   // The integration order computed in the contact
    KRATOS_REGISTER_VARIABLE( ELEMENT_POINTER )                             // A pointer to the element belonging to this condition
    KRATOS_REGISTER_VARIABLE( MORTAR_CONTACT_OPERATOR )                     // Mortar Contact Operator
    KRATOS_REGISTER_VARIABLE( ACTIVE_CHECK_FACTOR )                         // The factor employed to serach an active/inactive node
    KRATOS_REGISTER_VARIABLE( NORMAL_AUGMENTATION_FACTOR )                  // The constant that is considered for the check of active or inactive (when 0 it doesn't accept traction)
    KRATOS_REGISTER_VARIABLE( TANGENT_AUGMENTATION_FACTOR )                 // The constant that is considered for the check if the node is slip/stick
    KRATOS_REGISTER_VARIABLE( WEIGHTED_GAP )                                // The integrated gap employed in mortar formulation
    KRATOS_REGISTER_VARIABLE( WEIGHTED_SLIP )                               // The integrated slip employed in mortar formulation
    KRATOS_REGISTER_VARIABLE( WEIGHTED_FRICTION )                           // The integrated friction employed in mortar formulation
    KRATOS_REGISTER_VARIABLE( DELTA_NORMAL )                                // Directional derivative of the normal
    KRATOS_REGISTER_VARIABLE( AUXILIAR_ACTIVE )                             // Auxiliar boolean to check if the node is active or not
    KRATOS_REGISTER_VARIABLE( AUXILIAR_SLIP )                               // Auxiliar boolean to check if the node is stick or not
    KRATOS_REGISTER_VARIABLE( IS_ACTIVE_SET )                               // A bool storing whether the node is in the active set or not               
    KRATOS_REGISTER_VARIABLE( GAP_GP )                                      // A double storing the gap of the GP
    KRATOS_REGISTER_VARIABLE( SLIP_GP )                                     // A double storing the slip of the GP

    // Register the conditions
    // Mortar contact condition
    KRATOS_REGISTER_CONDITION( "MortarContactCondition2D2N", mMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "MortarContactCondition2D3N", mMortarContactCondition2D3N );
    KRATOS_REGISTER_CONDITION( "MortarContactCondition3D3N", mMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "MortarContactCondition3D6N", mMortarContactCondition3D6N );
    KRATOS_REGISTER_CONDITION( "MortarContactCondition3D4N", mMortarContactCondition3D4N );
    KRATOS_REGISTER_CONDITION( "MortarContactCondition3D8N", mMortarContactCondition3D8N );
    KRATOS_REGISTER_CONDITION( "MortarContactCondition3D9N", mMortarContactCondition3D9N );
}

}  // namespace Kratos.
