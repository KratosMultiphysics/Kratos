// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(CONDITION_MAP_DEFINED )
#define  CONDITION_MAP_DEFINED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/condition.h"

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
    
/** @brief Custom Point container to be used by the mapper
 */
class ConditionMap : public std::unordered_map<Condition::Pointer, bool, SharedPointerHasher<Condition::Pointer>, SharedPointerComparator<Condition::Pointer> >
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of ConditionMap
    KRATOS_CLASS_POINTER_DEFINITION( ConditionMap );

    typedef std::unordered_map<Condition::Pointer, bool, SharedPointerHasher<Condition::Pointer>, SharedPointerComparator<Condition::Pointer> > BaseType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    ConditionMap(){}

    /// Destructor
    virtual ~ConditionMap(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * It removes one particular condition from the map
     * @param pCond: The condition to remove
     */     
    void RemoveCondition(Condition::Pointer pCond)
    {
        BaseType::iterator set = find(pCond);
        if(set != end())
        {
            erase(set);
        }
    }
    
    /**
     * It adds one new condition
     * @param pCond: The condition to set
     */
    void AddNewCondition(Condition::Pointer pCond)
    {
        insert({pCond, true}); // True by default when adding a new one
    }
    
    /**
     * It adds one new condition, as active
     * @param pCond: The condition to set
     */
    void AddNewActiveCondition(Condition::Pointer pCond)
    {
        insert({pCond, true});
    }
    
    /**
     * It adds one new condition, as inactive
     * @param pCond: The condition to set
     */
    void AddNewInactiveCondition(Condition::Pointer pCond)
    {
        insert({pCond, false});
    }
    
    /**
     * It sets one particular condition as active or not
     * @param pCond: The condition to set
     * @param Active: The flag, true if active, false otherwise
     */
    void SetActive(Condition::Pointer pCond, const bool Active = true)
    {
        BaseType::iterator set = find(pCond);
        if(set != end())
        {
            set->second = Active;
        }
    }
    
    /**
     * It checks if one particular condition is active
     * @param pCond: The condition to check
     * @return True if it is active, false otherwise
     */
    bool IsActive(Condition::Pointer pCond) const 
    {
        BaseType::const_iterator set = find(pCond);
        return (set->second);
    }
    
    /**
     * It checks if at least one pair is active
     * @return True if at least one pair is active, false otherwise
     */
    bool AtLeastOnePairActive()
    {
        for ( auto it = begin(); it != end(); ++it )
        {
            if (it->second == true)
            {
                return true;
            }
        }
        
        return false;
    }
    
    /**
     * Print the map information
     */
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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class ConditionMap 

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // CONDITION_MAP_DEFINED  defined
