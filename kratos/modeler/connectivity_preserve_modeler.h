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
#include "includes/node.h"
#include "modeler/modeler.h"


namespace Kratos
{

/// A tool to generate a copy of a ModelPart, sharing the same nodes as the original.
class KRATOS_API(KRATOS_CORE) ConnectivityPreserveModeler : public Modeler
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(ConnectivityPreserveModeler);

    /// Constructor
    ConnectivityPreserveModeler();

    /// Destructor.
    virtual ~ConnectivityPreserveModeler() override;

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
    virtual void GenerateModelPart(
        ModelPart& OriginModelPart,
        ModelPart& DestinationModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition
    ) override;

private:

    void CheckVariableLists(ModelPart &rOriginModelPart, ModelPart &rDestinationModelPart);

    void ResetModelPart(ModelPart &rDestinationModelPart);

    void CopyCommonData(
        ModelPart &rOriginModelPart,
        ModelPart &rDestinationModelPart
    );

    void DuplicateElements(
        ModelPart &rOriginModelPart,
        ModelPart &rDestinationModelPart,
        Element const &rReferenceElement
    );

    void DuplicateConditions(
        ModelPart &rOriginModelPart,
        ModelPart &rDestinationModelPart,
        Condition const &rReferenceBoundaryCondition
    );

    void DuplicateCommunicatorData(
        ModelPart &rOriginModelPart,
        ModelPart &rDestinationModelPart
    );

    void DuplicateSubModelParts(
        ModelPart &rOriginModelPart,
        ModelPart &rDestinationModelPart
    );

    ConnectivityPreserveModeler & operator=(ConnectivityPreserveModeler const& rOther);

    /// Copy constructor.
    ConnectivityPreserveModeler(ConnectivityPreserveModeler const& rOther);


};

} // namespace Kratos.

#endif //KRATOS_GENERATE_MODEL_PART_MODELER_INCLUDED  defined
