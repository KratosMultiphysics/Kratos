//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define_registry.h"
#include "includes/model_part.h"
#include "modeler/modeler.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class ConnectivityPreserveModeler
 * @ingroup KratosCore
 * @brief This class is used to generate a copy of a ModelPart, sharing the same nodes as the original.
 * @details The elements and conditions of the new ModelPart will have the same connectivity (and id) as in the original ModelPart,
 * but their type is determined by the reference element and condition provided during construction.
 * Note that both ModelParts will share the same nodes, as well as ProcessInfo and tables.
 * SubModelParts and, in MPI, communicator data will be replicated in the new ModelPart.
 * @author Riccardo Rossi
 */
class KRATOS_API(KRATOS_CORE) ConnectivityPreserveModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ConnectivityPreserveModeler);

    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.KratosMultiphysics", Modeler, ConnectivityPreserveModeler)

    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.All", Modeler, ConnectivityPreserveModeler)

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor
    ConnectivityPreserveModeler() = default;

    /// Factory Constructor
    ConnectivityPreserveModeler(Model& rModeler, Parameters Settings);

    /// Copy constructor.
    ConnectivityPreserveModeler(ConnectivityPreserveModeler const& rOther) = delete;

    /// Destructor.
    ~ConnectivityPreserveModeler() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ConnectivityPreserveModeler & operator=(ConnectivityPreserveModeler const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Factory initialization
    Modeler::Pointer Create(Model& rModel, const Parameters Settings) const override;

    /**
     * @brief Generate a copy of rOriginModelPart in rDestinationModelPart, using the given element and condition types.
     * @details This function fills rDestinationModelPart using data obtained from rOriginModelPart. The elements
     *  and conditions of rDestinationModelPart part use the same connectivity (and id) as in
     *  rOriginModelPart but their type is determined by rReferenceElement and rReferenceBoundaryCondition.
     *  Note that both ModelParts will share the same nodes, as well as ProcessInfo and tables.
     *  SubModelParts and, in MPI, communicator data will be replicated in DestinationModelPart.
     * @param rOriginModelPart The source ModelPart.
     * @param rDestinationModelPart The ModelPart to be filled by this function.
     * @param rReferenceElement The Element type for rDestinationModelPart.
     * @param rReferenceBoundaryCondition The Condition type for rDestinationModelPart.
     */
    void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart,
        const Element& rReferenceElement,
        const Condition& rReferenceBoundaryCondition
        ) override;

    /**
     * @brief Generate a copy of rOriginModelPart in rDestinationModelPart, using the given element type.
     * @details This function fills rDestinationModelPart using data obtained from
     *  rOriginModelPart. The elements of rDestinationModelPart part use the
     *  same connectivity (and id) as in rOriginModelPart but their type is
     *  determined by rReferenceElement. In this variant, conditions are not
     *  copied.
     *  Note that both ModelParts will share the same nodes, as well as
     *  ProcessInfo and tables. SubModelParts and, in MPI, communicator data
     *  will be replicated in DestinationModelPart.
     * @param rOriginModelPart The source ModelPart.
     * @param rDestinationModelPart The ModelPart to be filled by this function.
     * @param rReferenceElement The Element type for rDestinationModelPart.
     */
    virtual void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart,
        const Element& rReferenceElement
        );

    /**
     * @brief Generate a copy of rOriginModelPart in rDestinationModelPart, using the given condition type.
     * @details This function fills rDestinationModelPart using data obtained from
     *  rOriginModelPart. The conditions of rDestinationModelPart part use the
     *  same connectivity (and id) as in rOriginModelPart but their type is
     *  determined by rReferenceCondition. In this variant, elements are not
     *  copied.
     *  Note that both ModelParts will share the same nodes, as well as
     *  ProcessInfo and tables. SubModelParts and, in MPI, communicator data
     *  will be replicated in DestinationModelPart.
     * @param rOriginModelPart The source ModelPart.
     * @param rDestinationModelPart The ModelPart to be filled by this function.
     * @param rReferenceCondition The Condition type for rDestinationModelPart.
     */
    virtual void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart,
        const Condition& rReferenceCondition
        );

    /**
     * @brief Generate a copy of rOriginModelPart in rDestinationModelPart.
     * @details This function fills rDestinationModelPart using data obtained from
     *  rOriginModelPart. The geometries of the rDestinationModelPart use
     *  the same connectivity (and id) as in rOriginModelPart and their type
     *  is determined according to the corresponding geometry of the entity.
     *  Note that both ModelParts will share the same nodes, geometries, as
     *  well as ProcessInfo and tables. SubModelParts and, in MPI, communicator
     *  data will be replicated in DestinationModelPart.
     * @param rOriginModelPart The source ModelPart.
     * @param rDestinationModelPart The ModelPart to be filled by this function
     */
    virtual void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart
        );

    /**
     * @brief Generate a copy of rOriginModelPart in rDestinationModelPart.
     * @details This function fills rDestinationModelPart using data obtained from
     *  rOriginModelPart. It is equivalent to one of the GenerateModelPart
     *  functions, depending on whether an element and/or a condition
     *  have been defined in the Parameters during construction.
     */
    void SetupModelPart() override;

    /// Defines the expected structure for the Parameters of this class.
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    ///@}
private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Checks if the variable lists between the origin and destination ModelParts match.
     * @details This function checks if the list of solution step variables in the destination ModelPart
     * matches with those in the origin ModelPart, and logs warnings if there are any mismatches.
     * @param rOriginModelPart The origin ModelPart.
     * @param rDestinationModelPart The destination ModelPart.
     */
    void CheckVariableLists(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const;

    /**
     * @brief Resets the destination ModelPart by setting the TO_ERASE flag for nodes, elements, and conditions.
     * @details This function marks all nodes, elements, and conditions in the destination ModelPart with the
     * TO_ERASE flag and removes them from all levels.
     * @param rDestinationModelPart The destination ModelPart to reset.
     */
    void ResetModelPart(ModelPart& rDestinationModelPart) const;

    /**
     * @brief Copies common data (such as solution step variables, process info, nodes, and geometries) 
     * from the origin ModelPart to the destination ModelPart.
     * @details This function handles copying essential data from the origin to the destination ModelPart. It 
     * performs checks to ensure the destination ModelPart is not a SubModelPart and raises errors 
     * if any attempts are made to modify attributes that would break the parent ModelPart.
     * @param rOriginModelPart The origin ModelPart to copy data from.
     * @param rDestinationModelPart The destination ModelPart to copy data to.
     */
    void CopyCommonData(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart
        ) const;

    /**
     * @brief Duplicates elements from the origin ModelPart to the destination ModelPart using a reference element.
     * @details This function duplicates elements from the origin to the destination ModelPart, reusing the 
     * geometry of the reference element for memory efficiency.
     * @param rOriginModelPart The origin ModelPart containing elements to duplicate.
     * @param rDestinationModelPart The destination ModelPart where the elements will be copied to.
     * @param rReferenceElement The reference element used for creating new elements.
     */
    void DuplicateElements(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Element& rReferenceElement
        ) const;

    /**
     * @brief Duplicates conditions from the origin ModelPart to the destination ModelPart using a reference condition.
     * @details This function duplicates conditions from the origin to the destination ModelPart, reusing the 
     * geometry of the reference condition for memory efficiency.
     * @param rOriginModelPart The origin ModelPart containing conditions to duplicate.
     * @param rDestinationModelPart The destination ModelPart where the conditions will be copied to.
     * @param rReferenceBoundaryCondition The reference condition used for creating new conditions.
     */
    void DuplicateConditions(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Condition& rReferenceBoundaryCondition
        ) const;

    /**
     * @brief Duplicates communicator data from the origin ModelPart to the destination ModelPart.
     * @details This function duplicates communicator data, including information about the mesh, colors, 
     * and node lists for parallel computing scenarios. The elements and conditions will be added 
     * later using the new elements.
     * @param rOriginModelPart The origin ModelPart from which communicator data will be copied.
     * @param rDestinationModelPart The destination ModelPart to which communicator data will be copied.
     */
    void DuplicateCommunicatorData(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart
        ) const;

    /**
     * @brief Duplicates sub-model parts recursively from the origin ModelPart to the destination ModelPart.
     * @details This function copies all sub-model parts, including nodes, elements, conditions, geometries, constraints,
     * and communicator data, ensuring that all child sub-model parts are duplicated as well.
     * @param rOriginModelPart The origin ModelPart containing sub-model parts to duplicate.
     * @param rDestinationModelPart The destination ModelPart to which the sub-model parts will be copied.
     */
    void DuplicateSubModelParts(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart
        ) const;

    ///@}
    ///@name Private members
    ///@{

    /// The pointer to the Model object.
    Model* mpModel = nullptr;

    ///@}
};

///@}

} // namespace Kratos.