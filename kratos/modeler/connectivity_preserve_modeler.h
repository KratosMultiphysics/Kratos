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

#if !defined( KRATOS_CONNECTIVITY_PRESERVE_MODELER_H )
#define  KRATOS_CONNECTIVITY_PRESERVE_MODELER_H



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "modeler/modeler.h"


namespace Kratos
{
///@}
///@name Kratos Classes
///@{

/// A tool to generate a copy of a ModelPart, sharing the same nodes as the original.
class KRATOS_API(KRATOS_CORE) ConnectivityPreserveModeler : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ConnectivityPreserveModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor
    ConnectivityPreserveModeler() = default;

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

    /// Generate a copy of rOriginModelPart in rDestinationModelPart, using the given element and condtion types.
    /** This function fills rDestinationModelPart using data obtained from rOriginModelPart. The elements
     *  and conditions of rDestinationModelPart part use the same connectivity (and id) as in
     *  rOriginModelPart but their type is determined by rReferenceElement and rReferenceBoundaryCondition.
     *  Note that both ModelParts will share the same nodes, as well as ProcessInfo and tables.
     *  SubModelParts and, in MPI, communicator data will be replicated in DestinationModelPart.
     *  @param rOriginModelPart The source ModelPart.
     *  @param rDestinationModelPart The ModelPart to be filled by this function.
     *  @param rReferenceElement The Element type for rDestinationModelPart.
     *  @param rReferenceBoundaryCondition The Condition type for rDestinationModelPart.
     */
    void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart,
        const Element& rReferenceElement,
        const Condition& rReferenceBoundaryCondition
    ) override;

    /// Generate a copy of rOriginModelPart in rDestinationModelPart, using the given element type.
    /** This function fills rDestinationModelPart using data obtained from
     *  rOriginModelPart. The elements of rDestinationModelPart part use the
     *  same connectivity (and id) as in rOriginModelPart but their type is
     *  determined by rReferenceElement. In this variant, conditions are not
     *  copied.
     *  Note that both ModelParts will share the same nodes, as well as
     *  ProcessInfo and tables. SubModelParts and, in MPI, communicator data
     *  will be replicated in DestinationModelPart.
     *  @param rOriginModelPart The source ModelPart.
     *  @param rDestinationModelPart The ModelPart to be filled by this function.
     *  @param rReferenceElement The Element type for rDestinationModelPart.
     */
    virtual void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart,
        const Element& rReferenceElement
    );

    /// Generate a copy of rOriginModelPart in rDestinationModelPart, using the given condition type.
    /** This function fills rDestinationModelPart using data obtained from
     *  rOriginModelPart. The conditions of rDestinationModelPart part use the
     *  same connectivity (and id) as in rOriginModelPart but their type is
     *  determined by rReferenceCondition. In this variant, elements are not
     *  copied.
     *  Note that both ModelParts will share the same nodes, as well as
     *  ProcessInfo and tables. SubModelParts and, in MPI, communicator data
     *  will be replicated in DestinationModelPart.
     *  @param rOriginModelPart The source ModelPart.
     *  @param rDestinationModelPart The ModelPart to be filled by this function.
     *  @param rReferenceCondition The Condition type for rDestinationModelPart.
     */
    virtual void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart,
        const Condition& rReferenceCondition
    );

    ///@}

private:
    ///@name Private member Variables
    ///@{

    IndexType mElementOffset;
    IndexType mConditionOffset;

    ///@}

    ///@name Private Operations
    ///@{

    void CheckVariableLists(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const;

    void ResetModelPart(ModelPart& rDestinationModelPart) const;

    void CopyCommonData(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart
    ) const;

    void DuplicateElements(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Element& rReferenceElement
    );

    void DuplicateConditions(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Condition& rReferenceBoundaryCondition
    );

    void DuplicateCommunicatorData(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart
    ) const;

    void DuplicateSubModelParts(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart
    ) const;

    template <typename TContainer>
    IndexType GetMaxId(const TContainer& rContainer, const Communicator& rCommunicator) const
    {
        KRATOS_TRY

        const int number_of_items = rContainer.size();

        IndexType max_item_id{0}, local_id{0};
        // using customized max reduction since windows doesnt support max reduction in omp
#pragma omp parallel firstprivate(local_id)
        {
#pragma omp for
            for (int i_item = 0; i_item < number_of_items; ++i_item)
            {
                const auto& r_item = *(rContainer.begin() + i_item);
                local_id = std::max(local_id, r_item.Id());
            }
#pragma omp critical
            {
                max_item_id = std::max(max_item_id, local_id);
            }
        }
        return rCommunicator.GetDataCommunicator().MaxAll(max_item_id);

        KRATOS_CATCH("");
    }

    ///@}
};

///@}

} // namespace Kratos.

#endif //KRATOS_GENERATE_MODEL_PART_MODELER_INCLUDED  defined
