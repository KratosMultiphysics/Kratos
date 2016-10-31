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

#if !defined(KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "solid_mechanics_application.h"
#include "structural_mechanics_application.h"
#include "contact_structural_mechanics_application.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
typedef array_1d<double,3> Vector3;

struct contact_container 
{
    Condition::Pointer                 condition;
  
    ~contact_container(){}
    
    void clear()
    {
        condition = nullptr;
    }
    
    void print()
    {
//        KRATOS_WATCH(condition);
       std::cout << " The condition: " << condition->Id() << " is MASTER: " << condition->Is(MASTER) << std::endl;
       std::cout << std::endl;
    }
    
    void save( Serializer& rSerializer ) const
    {
        rSerializer.save("condition",                               condition);
    }

    void load( Serializer& rSerializer )
    {
        rSerializer.load("condition",                              condition);
    }
};
// CONDITIONS
/* Mortar method */ 
KRATOS_DEFINE_VARIABLE( std::vector<contact_container>*, CONTACT_CONTAINERS ) // A vector of which contains the structure which defines the contact conditions
KRATOS_DEFINE_VARIABLE( Element::Pointer, ELEMENT_POINTER )                   // A pointer to the element belonging to this condition
KRATOS_DEFINE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                     // The integration order computed in the contact
KRATOS_DEFINE_VARIABLE( Matrix, MORTAR_CONTACT_OPERATOR )                     // Mortar Contact Operator
KRATOS_DEFINE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                         // The factor employed to serach an active/inactive node
KRATOS_DEFINE_VARIABLE( double, NORMAL_AUGMENTATION_FACTOR )                  // The constant that is considered for the check of active or inactive (when 0 it doesn't accept traction)
KRATOS_DEFINE_VARIABLE( double, TANGENT_AUGMENTATION_FACTOR )                 // The constant that is considered for the check if the node is slip/stick
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_GAP )                                // The integrated gap employed in mortar formulation
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_SLIP )                               // The integrated slip employed in mortar formulation
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_FRICTION )                           // The integrated friction employed in mortar formulation
KRATOS_DEFINE_VARIABLE( Matrix, DELTA_NORMAL )                                // Directional derivative of the normal
KRATOS_DEFINE_VARIABLE( bool, AUXILIAR_ACTIVE )                               // Auxiliar boolean to check if the node is active or not
KRATOS_DEFINE_VARIABLE( bool, AUXILIAR_SLIP )                                 // Auxiliar boolean to check if the node is stick or not
KRATOS_DEFINE_VARIABLE( bool, IS_ACTIVE_SET )                                 // A bool storing whether the node is in the active set or not
KRATOS_DEFINE_VARIABLE( double, GAP_GP )                                      // A double storing the gap of the GP
KRATOS_DEFINE_VARIABLE( double, SLIP_GP )                                     // A double storing the slip of the GP
}       

#endif	/* KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
