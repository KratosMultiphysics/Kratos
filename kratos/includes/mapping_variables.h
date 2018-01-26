//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MAPPING_VARIABLES_H_INCLUDED)
#define KRATOS_MAPPING_VARIABLES_H_INCLUDED

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "includes/kratos_components.h"

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

namespace Kratos
{
///@name Kratos Globals
///@{
   
///@}
///@name Type Definitions
///@{

    typedef std::size_t IndexType;
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{
    
    /** @brief Custom Point container to be used by the mapper
    */
    class IndexSet : public std::unordered_set<IndexType>
    {
    public:

        ///@name Type Definitions
        ///@{
        /// Counted pointer of IndexSet
        KRATOS_CLASS_POINTER_DEFINITION( IndexSet );

        typedef std::unordered_set<IndexType> BaseType;
        
        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructors
        IndexSet(){}

        /// Destructor
        virtual ~IndexSet(){}

        ///@}
        ///@name Operators
        ///@{
        
        ///@}
        ///@name Operations
        ///@{
        
        /**
        * It checks if an ID exists in the map
        * @param IndexOrigin The condition ID to remove
        * @return If the ID already exists or not
        */     
        bool Has(const IndexType IndexOrigin)
        {
            BaseType::iterator set = find(IndexOrigin);
            if(set != end())
            {
                return true;
            }
            
            return false;
        }
        
        /**
        * It adds a new ID to the map
        * @param IndexOrigin The condition ID to remove
        */     
        void AddId(const IndexType IndexOrigin)
        {
            insert({IndexOrigin});
        }
        
        /**
        * It removes one particular pair from the map
        * @param IndexOrigin The condition ID to remove
        */     
        void RemoveId(const IndexType IndexOrigin)
        {
            BaseType::iterator set = find(IndexOrigin);
            if(set != end())
            {
                erase(set);
            }
        }

        /**
         * Turn back information as a string.
         */
        std::string Info() const //override
        {
            std::stringstream buffer;
            for ( auto it = begin(); it != end(); ++it )
                buffer << "The condition " << *it << std::endl;
            return buffer.str();
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
    }; // Class IndexSet 
    
    KRATOS_DEFINE_VARIABLE( IndexSet::Pointer, INDEX_SET )         // An unordened map of which contains the indexes with the paired 
    KRATOS_DEFINE_VARIABLE( double, TANGENT_FACTOR )               // The factor between the tangent and normal behaviour

} // namespace Kratos

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_MAPPING_VARIABLES_H_INCLUDED defined
