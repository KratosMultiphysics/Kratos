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

        /// Length of array definition
        constexpr static std::size_t LengthArray = static_cast<std::size_t>(GeometryData::KratosGeometryType::NumberOfGeometryTypes);

        /// Pointer definition of ReplaceElementsAndConditionsProcess
        KRATOS_CLASS_POINTER_DEFINITION(EntitityIdentifier);

        ///@}
        ///@name Life Cycle
        ///@{

        /**
         * @brief Default constructor
         */
        EntitityIdentifier() = default;

        /**
         * @brief Constructor
         * @param rName The name of the entity.
         */
        EntitityIdentifier(const std::string& rName);

        /**
         * @brief Copy constructor
         */
        EntitityIdentifier(const EntitityIdentifier& rOther)
            : mTypes(rOther.mTypes),
              mIsInitialized(rOther.mIsInitialized)
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
            mTypes = rOther.mTypes;
            mIsInitialized = rOther.mIsInitialized;

            return *this;
        }

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Checks if the object is initialized.
         * @return true if the object is initialized, false otherwise.
         */
        bool IsInitialized() const;

        /**
         * @brief Get the prototype entity.
         * @param rGeometry The reference to the geometry.
         * @return true there is a prototype for the provided entity.
         * @return false there is no prototype for the provided entity.
         */
        bool HasPrototypeEntity(const GeometryType& rGeometry) const;

        /**
         * @brief Get the prototype entity.
         * @param pGeometry The pointer to the geometry.
         * @return const TEntity& The prototype entity.
         */
        const TEntity& GetPrototypeEntity(typename GeometryType::Pointer pGeometry) const;

        /**
         * @brief Get the prototype entity.
         * @param rGeometry The reference to the geometry.
         * @return const TEntity& The prototype entity.
         */
        const TEntity& GetPrototypeEntity(const GeometryType& rGeometry) const;

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
        ///@name Private Member Variables
        ///@{

        std::array<const TEntity*, LengthArray> mTypes;          /// The settings of the entities
        bool mIsInitialized = false;                             /// If the object is initialized

        ///@}
        ///@name Private Operations
        ///@{

        /**
         * @brief Returns the name of the entity type for the given entity.
         * @return the name of the entity type as a string
         * @throws std::logic_error if the entity type is not supported
         */
        std::string GetEntityTypeName() const;

        /**
         * @brief Generate single type prototype entity.
         * @param rName The name to generate single type from
         * @throws std::runtime_error if entity name is not found in KratosComponents
         */
        void GenerateSingleType(const std::string& rName);

        /**
         * @brief Generate multiple types entities map.
         * @param rName The name to generate multiple types from
         * @throws std::runtime_error if entity name is not found in KratosComponents
         */
        void GenerateMultipleTypes(const std::string& rName);

        /**
         * @brief Generate templated types entities map.
         * @param rName The name to generate templated types from
         * @throws std::runtime_error if entity name is not found in KratosComponents
         */
        void GenerateTemplatedTypes(const std::string& rName);

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
