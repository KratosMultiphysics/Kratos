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

#if !defined(KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include <unordered_map>

namespace Kratos
{

#if !defined(SHARED_POINTER_HASHER)
#define SHARED_POINTER_HASHER
template<class TSharedPointer>
struct SharedPointerHasher
{
    size_t operator()(const TSharedPointer& pCond) const
    {
        return (size_t)pCond.get(); 
    }
};
#endif
#if !defined(SHARED_POINTER_COMPARATOR)
#define SHARED_POINTER_COMPARATOR
template<class TSharedPointer>
struct SharedPointerComparator
{
    bool operator()(const TSharedPointer& first, const TSharedPointer& second) const
    {
        return *first == *second;
    }
};
#endif

typedef array_1d<double,3> Vector3;
typedef std::unordered_map<Condition::Pointer, bool, SharedPointerHasher<Condition::Pointer>, SharedPointerComparator<Condition::Pointer> > ConditionHashMap;

struct contact_container // TODO: Remove this, deprecated
{
    Condition::Pointer condition;
    bool             active_pair;
  
    ~contact_container(){}
    
    void clear()
    {
        condition   = nullptr;
        active_pair = false;
    }
    
    void print()
    {
       std::cout << " The condition: " << condition->Id() << " is MASTER: " << condition->Is(MASTER) << "ACTIVE: " << active_pair << std::endl;
       std::cout << std::endl;
    }
    
    void save( Serializer& rSerializer ) const
    {
        rSerializer.save("condition",     condition);
        rSerializer.save("active_pair", active_pair);
    }

    void load( Serializer& rSerializer )
    {
        rSerializer.load("condition",     condition);
        rSerializer.load("active_pair", active_pair);
    }
};

struct ConditionMap : ConditionHashMap
{
    ~ConditionMap(){}
    
    typedef ConditionHashMap BaseType;
    
    void RemoveCondition(Condition::Pointer pCond)
    {
        BaseType::iterator Set = find(pCond);
        if(Set != end())
        {
            erase(Set);
        }
    }
    
    void AddNewCondition(Condition::Pointer pCond)
    {
        insert({pCond, true}); // True by default when adding a new one
    }
    
    void SetActive(Condition::Pointer pCond, const bool Active)
    {
        BaseType::iterator Set = find(pCond);
        if(Set != end())
        {
            Set->second = Active;
        }
    }
    
    bool IsActive(Condition::Pointer pCond)
    {
        BaseType::const_iterator Set = find(pCond);
        return (Set->second);
    }
    
    void print()
    {
        for ( auto it = begin(); it != end(); ++it )
        {
            std::cout << "The condition " << (it->first)->Id() << " is ACTIVE: " << it->second;
            
            KRATOS_WATCH((it->first)->GetGeometry());
        }
    }
    
    void save( Serializer& rSerializer ) const
    {
        // TODO: Fill if necessary
    }

    void load( Serializer& rSerializer )
    {
        // TODO: Fill if necessary
    }
};

// CONDITIONS
/* Mortar method */ 
KRATOS_DEFINE_VARIABLE( std::vector<contact_container>*, CONTACT_CONTAINERS )                                                   // A vector of which contains the structure which defines the contact conditions // TODO: Remove this, deprecated
KRATOS_DEFINE_VARIABLE( ConditionMap*, CONTACT_SETS )                                                                              // An unordened map of which contains the structure which defines the contact conditions
KRATOS_DEFINE_VARIABLE( Element::Pointer, ELEMENT_POINTER )                                                                     // A pointer to the element belonging to this condition
KRATOS_DEFINE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                                                                       // The integration order computed in the contact
KRATOS_DEFINE_VARIABLE( Matrix, MORTAR_CONTACT_OPERATOR )                                                                       // Mortar Contact Operator
KRATOS_DEFINE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                                                                           // The factor employed to serach an active/inactive node

/* The complementary values */
// NOTE: This will be eventually not necessary
KRATOS_DEFINE_VARIABLE( double, NORMAL_AUGMENTATION_FACTOR )                                                                    // The constant that is considered for the check of active or inactive (when 0 it doesn't accept traction)
KRATOS_DEFINE_VARIABLE( double, TANGENT_AUGMENTATION_FACTOR )                                                                   // The constant that is considered for the check if the node is slip/stick

/* Weighted values */
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_GAP )                                                                                  // The integrated gap employed in mortar formulation
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_SLIP )                                                                                 // The integrated slip employed in mortar formulation

/* Matrix to store the derivatives of the normal */
KRATOS_DEFINE_VARIABLE( Matrix, DELTA_NORMAL )                                                                                  // Directional derivative of the normal

/* Auxiliar booleans to store the change in active/inactive slip/stick */
KRATOS_DEFINE_VARIABLE( bool, AUXILIAR_ACTIVE )                                                                                 // Auxiliar boolean to check if the node is active or not
KRATOS_DEFINE_VARIABLE( bool, AUXILIAR_SLIP )                                                                                   // Auxiliar boolean to check if the node is stick or not

/* For ALM mortar condition */
KRATOS_DEFINE_VARIABLE( double, PENALTY_FACTOR )                                                                                // The penalty factor for the ALM
KRATOS_DEFINE_VARIABLE( double, SCALE_FACTOR )                                                                                  // The scale factor for the ALM

/* For mesh tying mortar condition */
KRATOS_DEFINE_VARIABLE( std::string, TYING_VARIABLE )                                                                           // The variable name for the mesh tying  

}       

#endif	/* KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
