//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Sonja Schneider
//

// From S. Schneider Implementation of a Chimera Technique TUM Master Thesis 2015

#include "interpolation_utility.h"
#include "includes/kratos_components.h"

namespace Kratos {

template<std::size_t TDim>
void InterpolationUtility<TDim>::ExtractBoundaryFaces(ModelPart &rModelPart, ModelPart &rOutModelPart)
{
    rOutModelPart.GetNodalSolutionStepVariablesList() = rModelPart.GetNodalSolutionStepVariablesList();
    rOutModelPart.SetBufferSize(rModelPart.GetBufferSize());
    rOutModelPart.SetNodes(rModelPart.pNodes());
    rOutModelPart.SetProperties(rModelPart.pProperties());

    const Condition& rRefCondition = KratosComponents<Condition>::Get("ChimeraFluidCouplingCondition3D");
    int Id = 1;

    for (ModelPart::ElementIterator itElem = rModelPart.ElementsBegin(); itElem != rModelPart.ElementsEnd(); itElem++)
    {
        Condition::NodesArrayType Nodes = Condition::NodesArrayType();
        Geometry< Node<3> >& rGeom = itElem->GetGeometry();
        for( std::size_t i = 0; i != rGeom.PointsNumber(); i++)
        {
            if (rGeom[i].FastGetSolutionStepValue(FLAG_VARIABLE) == 1.0)
            {
                Nodes.push_back(rGeom(i));
            }
        }

        if (Nodes.size()==3)
        {
            rOutModelPart.Conditions().push_back( rRefCondition.Create(Id++,Nodes,itElem->pGetProperties()) );
        }
    }
}

template class InterpolationUtility<2>;
template class InterpolationUtility<3>;

}
