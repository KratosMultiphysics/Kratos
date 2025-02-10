//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "assembly_utilities.h"

namespace Kratos
{

void AssemblyUtilities::SetUpSparseGraph(
    const ModelPart& rModelPart,
    SparseGraph<typename ModelPart::IndexType>& rSparseGraph)
{
    if (!rSparseGraph.IsEmpty()) {
        KRATOS_WARNING("AssemblyUtilities") << "Provided sparse graph is not empty and will be cleared." << std::endl;
        rSparseGraph.Clear();
    }

    Element::EquationIdVectorType eq_ids;
    for (auto& r_elem : rModelPart.Elements()) {
        r_elem.EquationIdVector(eq_ids, rModelPart.GetProcessInfo());
        rSparseGraph.AddEntries(eq_ids);
    }
    for (auto& r_cond : rModelPart.Conditions()) {
        r_cond.EquationIdVector(eq_ids, rModelPart.GetProcessInfo());
        rSparseGraph.AddEntries(eq_ids);
    }
}

}  // namespace Kratos.
