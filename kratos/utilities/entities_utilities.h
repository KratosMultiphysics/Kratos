//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/string_utilities.h"

namespace Kratos
{
/**
 * @namespace EntitiesUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for the computation of entities functions in a efficient way
 * @author Vicente Mataix Ferrandiz
 * @author Philipp Bucher
 */
namespace EntitiesUtilities
{
    // Enum class representing different types of definitions
    enum class DefinitionType {
        Single,     // Single definition
        Multiple,   // Multiple definitions
        Templated   // Templated definition
    };

    /**
     * @brief Template struct for entity identifier.
     * @brief This struct is used to identify and retrieve entity types based on their names and definitions.
     * @tparam TEntity The entity type.
     */
    template<class TEntity>
    class KRATOS_API(KRATOS_CORE) EntitityIdentifier
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Geometry type definition
        using GeometryType = typename TEntity::GeometryType;

        /// Pointer definition of ReplaceElementsAndConditionsProcess
        KRATOS_CLASS_POINTER_DEFINITION(EntitityIdentifier);

        ///@}
        ///@name Life Cycle
        ///@{
            
        /** 
         * @brief Default constructor 
         */
        EntitityIdentifier() :
            mName("")
        {
        }

        /**
         * @brief Constructor
         * @param rName The name of the entity.
         */
        EntitityIdentifier(const std::string& rName);

        /**
         * @brief Copy constructor
         */
        EntitityIdentifier(const EntitityIdentifier& rOther)
            : mpPrototypeEntity(rOther.mpPrototypeEntity),
              mDefinitionType(rOther.mDefinitionType),
              mName(rOther.mName),
              mTypes(rOther.mTypes)
        {
        }

        ///@}
        ///@name Operators
        ///@{

        /**
         * @brief Assignment operator.
         */
        EntitityIdentifier& operator=(const EntitityIdentifier& rOther)
        {
            mpPrototypeEntity = rOther.mpPrototypeEntity;
            mDefinitionType = rOther.mDefinitionType;
            mName = rOther.mName;
            mTypes = rOther.mTypes;

            return *this;
        }

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Checks if the object is initialized.
         * @return true if the object is initialized, false otherwise.
         */
        bool IsInitialized();

        /**
         * @brief Get the prototype entity.
         * @param pGeometry The pointer to the geometry.
         * @return const TEntity& The prototype entity.
         */
        const TEntity& GetPrototypeEntity(typename GeometryType::Pointer pGeometry);

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const
        {
            return "EntitityIdentifier";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "EntitityIdentifier";
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const;

        ///@}
    private:
        ///@name Member Variables
        ///@{

        typename TEntity::Pointer mpPrototypeEntity = nullptr;                    /// The prototype entity
        DefinitionType mDefinitionType = DefinitionType::Single;                  /// The type of definition
        std::string mName;                                                        /// The name of the entity
        std::unordered_map<GeometryData::KratosGeometryType, const TEntity*> mTypes; /// The settings of the entities

        ///@}
    }; ///  Class EntitityIdentifier

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    template<class TEntity>
    inline std::istream& operator >> (std::istream& rIStream,
                                    EntitityIdentifier<TEntity>& rThis);

    /// output stream function
    template<class TEntity>
    inline std::ostream& operator << (std::ostream& rOStream,
                                    const EntitityIdentifier<TEntity>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    /**
     * @brief This method initializes all the active entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method calls InitializeSolution for all the entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeSolutionStepAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method calls FinalizeSolutionStep for all the entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) FinalizeSolutionStepAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method calls InitializeNonLinearIteration for all the entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeNonLinearIterationAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method calls FinalizeNonLinearIteration for all the entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) FinalizeNonLinearIterationAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method returns the appropriate TEntitytype container (elements, conditions, and nodes) from model part
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    KRATOS_API(KRATOS_CORE) PointerVectorSet<TEntityType, IndexedObject>& GetEntities(ModelPart& rModelPart);

    /**
     * @brief This method initializes all the active entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void InitializeEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // Array of entities
        auto& r_entities_array = GetEntities<TEntityType>(rModelPart);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            r_entities_array,
            [&r_current_process_info](TEntityType& rEntity) {
                // If the entity is active
                if (rEntity.IsActive()) {
                    rEntity.Initialize(r_current_process_info);
                }
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This method calls InitializeSolutionStep for all the entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void InitializeSolutionStepEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            GetEntities<TEntityType>(rModelPart),
            [&r_current_process_info](TEntityType& rEntity){
                rEntity.InitializeSolutionStep(r_current_process_info);
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This method calls FinalizeSolutionStep for all the entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void FinalizeSolutionStepEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            GetEntities<TEntityType>(rModelPart),
            [&r_current_process_info](TEntityType& rEntity){
                rEntity.FinalizeSolutionStep(r_current_process_info);
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This method calls InitializeNonLinearIteration for all the entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void InitializeNonLinearIterationEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            GetEntities<TEntityType>(rModelPart),
            [&r_current_process_info](TEntityType& rEntity){
                rEntity.InitializeNonLinearIteration(r_current_process_info);
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This method calls FinalizeNonLinearIteration for all the entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void FinalizeNonLinearIterationEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            GetEntities<TEntityType>(rModelPart),
            [&r_current_process_info](TEntityType& rEntity){
                rEntity.FinalizeNonLinearIteration(r_current_process_info);
            }
        );

        KRATOS_CATCH("")
    }

    ///@}

}; // namespace EntitiesUtilities
}  // namespace Kratos
