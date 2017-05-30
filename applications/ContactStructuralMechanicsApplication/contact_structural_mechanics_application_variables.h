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
#include <unordered_set>

namespace Kratos
{

#if !defined(SHARED_POINTER_HASHER)
#define SHARED_POINTER_HASHER
template<class TSharedPointer>
struct SharedPointerHasher
{
    size_t operator()(const TSharedPointer& pCond) const
    {
        return reinterpret_cast<size_t>(pCond.get());
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
		return first.get() == second.get();
    }
};
#endif

typedef array_1d<double,3> Vector3;

struct ConditionSet : std::unordered_set<Condition::Pointer, SharedPointerHasher<Condition::Pointer>, SharedPointerComparator<Condition::Pointer> >
{
    ~ConditionSet(){}
    
    typedef std::unordered_set<Condition::Pointer, SharedPointerHasher<Condition::Pointer>, SharedPointerComparator<Condition::Pointer> > BaseType;
    
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
        insert(pCond);
    }
    
    void print()
    {
        for (auto it = begin(); it != end(); ++it )
        {            
            KRATOS_WATCH((*it)->GetGeometry());
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
KRATOS_DEFINE_VARIABLE( boost::shared_ptr<ConditionSet>, CONTACT_SETS )         // An unordened map of which contains the structure which defines the contact conditions
KRATOS_DEFINE_VARIABLE( Element::Pointer, ELEMENT_POINTER )                     // A pointer to the element belonging to this condition
KRATOS_DEFINE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                       // The integration order computed in the contact
KRATOS_DEFINE_VARIABLE( Matrix, MORTAR_CONTACT_OPERATOR )                       // Mortar Contact Operator
KRATOS_DEFINE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                           // The factor employed to serach an active/inactive node

/* Weighted values */
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_GAP )                                  // The integrated gap employed in mortar formulation
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_SLIP )                                 // The integrated slip employed in mortar formulation
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_FRICTION )                             // The integrated friction coefficient
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_SCALAR_RESIDUAL )                      // The integrated scalar residual  
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_VECTOR_RESIDUAL )           // The integrated vector residual         

/* Matrix to store the derivatives of the normal */
KRATOS_DEFINE_VARIABLE( Matrix, DELTA_NORMAL )                                  // Directional derivative of the normal

/* For ALM mortar condition */
KRATOS_DEFINE_VARIABLE( double, AUGMENTED_NORMAL_CONTACT_PRESSURE )            // The resultant augmented pressure in the normal direction
KRATOS_DEFINE_VARIABLE( double, AUGMENTED_TANGENT_CONTACT_PRESSURE )           // The resultant augmented pressure in the tangent direction
KRATOS_DEFINE_VARIABLE( double, PENALTY_PARAMETER )                            // The penalty factor for the ALM
KRATOS_DEFINE_VARIABLE( double, SCALE_FACTOR )                                 // The scale factor for the ALM
KRATOS_DEFINE_VARIABLE( double, TANGENT_FACTOR )                               // The proportion between the tangent and normal penalty
KRATOS_DEFINE_VARIABLE( bool, CONSIDER_NORMAL_VARIATION )                      // A value used to check if consider normal variation or not
KRATOS_DEFINE_VARIABLE( bool, CONSIDER_PAIR_VARIATION )                        // A value used to check if consider variation or not in the active inactive pairs

/* For mesh tying mortar condition */
KRATOS_DEFINE_VARIABLE( std::string, TYING_VARIABLE )                          // The variable name for the mesh tying  

}       

#endif	/* KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
