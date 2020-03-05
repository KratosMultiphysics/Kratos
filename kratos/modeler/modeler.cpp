//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//


// Project includes
#include "modeler.h"


namespace Kratos
{
    ///@name Operations
    ///@{

    void Modeler::GenerateModelPart(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for mesh generation." << std::endl;
    }

    void Modeler::GenerateMesh(
        ModelPart& ThisModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for mesh generation." << std::endl;
    }

    void Modeler::GenerateNodes(ModelPart& ThisModelPart)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for node generation." << std::endl;
    }

    ///@}
    ///@name Generate Elements and Conditions
    ///@{



    ///@}
}