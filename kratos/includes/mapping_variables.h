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
#include <unordered_map>

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
    
    /**
     * @class IndexDatabase
     * @ingroup KratosCore
     * @brief Base class to derive common methods
     * @author Vicente Mataix Ferrandiz
    */
    class IndexDatabase
    {
    public:
        ///@name Type Definitions
        ///@{
        /// Counted pointer of IndexDatabase
        KRATOS_CLASS_POINTER_DEFINITION( IndexDatabase );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructors
        IndexDatabase(){}

        /// Destructor
        virtual ~IndexDatabase(){}

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

       /**
        * @brief It checks if an ID exists in the map
        * @param IndexOrigin The condition ID to remove
        * @return If the ID already exists or not
        */
        virtual bool Has(const IndexType IndexOrigin) {return false;}

       /**
        * @brief It adds a new ID to the map
        * @param IndexOrigin The condition ID to remove
        * @param IndexNewEntity The new created entity ID
        */
        virtual void AddId(
            const IndexType IndexOrigin,
            const IndexType IndexNewEntity = 0
            ) {}

       /**
        * @brief It removes one particular pair from the map
        * @param IndexOrigin The condition ID to remove
        */
        virtual void RemoveId(const IndexType IndexOrigin) {}

        /**
        * @brief It returns the new created entity ID
        * @param IndexOrigin The condition ID to remove
        * @return The new created entity ID
        */
        virtual IndexType GetNewEntityId(const IndexType IndexOrigin) {return 0;}

        /**
        * @brief It sets the new created entity ID
        * @param IndexOrigin The condition ID to remove
        * @param IndexNewEntity The new created entity ID
        */
        virtual void SetNewEntityId(
            const IndexType IndexOrigin,
            const IndexType IndexNewEntity = 0
            ) {}

        /**
         * @brief Turn back information as a string.
         */
        virtual std::string Info() const {return "IndexDatabase";}

        virtual void save( Serializer& rSerializer ) const {}

        virtual void load( Serializer& rSerializer ) {}
    }; // Class IndexDatabase

    /**
     * @class IndexSet
     * @ingroup KratosCore
     * @brief Custom unordered set container to be used by the mapper
     * @details Contains a set of IDs of paired conditions
     * @author Vicente Mataix Ferrandiz
    */
    class IndexSet
        : public std::unordered_set<IndexType>, public IndexDatabase
    {
    public:
        ///@name Type Definitions
        ///@{
        /// Counted pointer of IndexSet
        KRATOS_CLASS_POINTER_DEFINITION( IndexSet );

        typedef std::unordered_set<IndexType> BaseType;
        typedef iterator IteratorType;

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
        * @brief Returns the id corresponding a iterator
        * @param ThisIterator The iterator of the class
        * @return The ID
        */
        IndexType GetId(IteratorType ThisIterator)
        {
            return *ThisIterator;
        }

        /**
        * @brief Returns the new entity id corresponding a iterator
        * @param ThisIterator The iterator of the class
        * @return The ID of the new generated entity
        */
        IndexType GetOtherId(IteratorType ThisIterator)
        {
            return 0;
        }

        /**
        * @brief It checks if an ID exists in the map
        * @param IndexOrigin The condition ID to remove
        * @return If the ID already exists or not
        */
        bool Has(const IndexType IndexOrigin) override
        {
            BaseType::iterator set = find(IndexOrigin);
            if(set != end())
                return true;

            return false;
        }

        /**
        * @brief It adds a new ID to the map
        * @param IndexOrigin The condition ID to remove
        * @param IndexNewEntity The new created entity ID
        */
        void AddId(
            const IndexType IndexOrigin,
            const IndexType IndexNewEntity = 0
            ) override
        {
            insert({IndexOrigin});
        }

        /**
        * @brief It removes one particular pair from the map
        * @param IndexOrigin The condition ID to remove
        */
        void RemoveId(const IndexType IndexOrigin) override
        {
            BaseType::iterator set = find(IndexOrigin);
            if(set != end())
                erase(set);
        }

        /**
         * @brief Turn back information as a string.
         */
        std::string Info() const override
        {
            std::stringstream buffer;
            for ( auto it = begin(); it != end(); ++it )
                buffer << "The condition " << *it << std::endl;
            return buffer.str();
        }

        void save( Serializer& rSerializer ) const override
        {
//             KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        }

        void load( Serializer& rSerializer ) override
        {
//             KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        }
    }; // Class IndexSet

    /**
     * @class IndexMap
     * @ingroup KratosCore
     * @brief Custom unordered map container to be used by the mapper
     * @details Contains a map of IDs of paired conditions, and the new created condition
     * @author Vicente Mataix Ferrandiz
    */
    class IndexMap
        : public std::unordered_map<IndexType, IndexType>, public IndexDatabase
    {
    public:
        ///@name Type Definitions
        ///@{
        /// Counted pointer of IndexMap
        KRATOS_CLASS_POINTER_DEFINITION( IndexMap );

        typedef std::unordered_map<IndexType, IndexType> BaseType;
        typedef iterator IteratorType;
        
        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructors
        IndexMap(){}

        /// Destructor
        virtual ~IndexMap(){}

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

        /**
        * @brief Returns the id corresponding a iterator
        * @param ThisIterator The iterator of the class
        * @return The ID
        */
        IndexType GetId(IteratorType ThisIterator)
        {
            return ThisIterator->first;
        }

        /**
        * @brief Returns the new entity id corresponding a iterator
        * @param ThisIterator The iterator of the class
        * @return The ID of the new generated entity
        */
        IndexType GetOtherId(IteratorType ThisIterator)
        {
            return ThisIterator->second;
        }

        /**
        * @brief It checks if an ID exists in the map
        * @param IndexOrigin The condition ID to remove
        * @return If the ID already exists or not
        */     
        bool Has(const IndexType IndexOrigin) override
        {
            BaseType::iterator map = find(IndexOrigin);
            if(map != end())
                return true;
            
            return false;
        }

        /**
        * @brief It adds a new ID to the map
        * @param IndexOrigin The condition ID to remove
        * @param IndexNewEntity The new created entity ID
        */     
        void AddId(
            const IndexType IndexOrigin,
            const IndexType IndexNewEntity = 0
            ) override
        {
            insert({IndexOrigin, IndexNewEntity});
        }
        
        /**
        * @brief It removes one particular pair from the map
        * @param IndexOrigin The condition ID to remove
        */     
        void RemoveId(const IndexType IndexOrigin) override
        {
            BaseType::iterator map = find(IndexOrigin);
            if(map != end())
                erase(map);
        }

        /**
        * @brief It returns the new created entity ID
        * @param IndexOrigin The condition ID to remove
        * @return The new created entity ID
        */
        IndexType GetNewEntityId(const IndexType IndexOrigin) override
        {
            BaseType::iterator map = find(IndexOrigin);
            if(map != end())
                return map->second;

            return 0;
        }

        /**
        * @brief It sets the new created entity ID
        * @param IndexOrigin The condition ID to remove
        * @param IndexNewEntity The new created entity ID
        */
        void SetNewEntityId(
            const IndexType IndexOrigin,
            const IndexType IndexNewEntity
            ) override
        {
            BaseType::iterator map = find(IndexOrigin);
            if(map != end())
                map->second = IndexNewEntity;
        }

        /**
         * @brief Turn back information as a string.
         */
        std::string Info() const override
        {
            std::stringstream buffer;
            for ( auto it = begin(); it != end(); ++it )
                buffer << "The condition " << it->first << " related with the new condition " << it->second << std::endl;
            return buffer.str();
        }
        
        void save( Serializer& rSerializer ) const override
        {
//             KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        }

        void load( Serializer& rSerializer ) override
        {
//             KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        }
    }; // Class IndexMap
    
    KRATOS_DEFINE_VARIABLE( IndexSet::Pointer, INDEX_SET )         // An unordened set of which contains the indexes with the paired
    KRATOS_DEFINE_VARIABLE( IndexMap::Pointer, INDEX_MAP )         // An unordened map of which contains the indexes with the paired
    KRATOS_DEFINE_VARIABLE( double, TANGENT_FACTOR )               // The factor between the tangent and normal behaviour

} // namespace Kratos

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_MAPPING_VARIABLES_H_INCLUDED defined
