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

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/parallel_utilities.h"
#include "processes/find_conditions_neighbours_process.h"

namespace Kratos
{

FindConditionsNeighboursProcess::FindConditionsNeighboursProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : mrModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))
{
    // Checking MPI
    KRATOS_ERROR_IF(mrModelPart.IsDistributed()) << "ModelPart cannot be distributed!. Current implementation is serial only" << std::endl;

    // Now validate against defaults -- this also ensures no type mismatch
    Parameters default_parameters = GetDefaultParameters();
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    // Setting the rest
    mAverageConditions = ThisParameters["average_conditions"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

FindConditionsNeighboursProcess::FindConditionsNeighboursProcess(
    ModelPart& rModelPart,
    const int Dim,
    const unsigned int AverageConditions
    ) : mrModelPart(rModelPart),
        mAverageConditions(AverageConditions),
        mDim(Dim)
{
    // Checking MPI
    KRATOS_ERROR_IF(mrModelPart.IsDistributed()) << "ModelPart cannot be distributed!. Current implementation is serial only" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void FindConditionsNeighboursProcess::Execute()
{
    // Compute dimension
    ComputeDimension();

    // Entities arrays
    auto& r_nodes_array = mrModelPart.Nodes();
    auto& r_conditions_array = mrModelPart.Conditions();

    // First of all the neighbour nodes and conditions array are initialized to the guessed size and empties the old entries
    block_for_each(r_nodes_array, [this](Node& rNode){
        GlobalPointersVector<Condition> neighbour_conditions;
        neighbour_conditions.reserve(mAverageConditions);
        rNode.SetValue(NEIGHBOUR_CONDITIONS, neighbour_conditions);
    });
    block_for_each(r_conditions_array, [this](Condition& rCond){
        GlobalPointersVector<Condition> neighbour_conditions;
        neighbour_conditions.resize(mDim);
        rCond.SetValue(NEIGHBOUR_CONDITIONS, neighbour_conditions);
    });

    // Add the neighbour conditions to all the nodes in the mesh
    for(auto it_cond = r_conditions_array.begin(); it_cond!=r_conditions_array.end(); it_cond++) {
        auto& r_geom = it_cond->GetGeometry();
        for(unsigned int i = 0; i < r_geom.size(); i++) {
            (r_geom[i].GetValue(NEIGHBOUR_CONDITIONS)).push_back( Condition::WeakPointer(*(it_cond.base())));
        }
    }

    // Shrink to fit
    block_for_each(r_nodes_array, [](Node& rNode) {
        rNode.GetValue(NEIGHBOUR_CONDITIONS).shrink_to_fit();
    });

    // Adding the neighbouring conditions to the condition loop over faces
    if (mDim == 3) {
        for(auto it_cond = r_conditions_array.begin(); it_cond!=r_conditions_array.end(); it_cond++) {
            // Face nodes
            auto& r_geom = it_cond->GetGeometry();
            if (r_geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) { // 3D Triangle
                // Vector of the 3 faces around the given face
                auto& r_neighb_faces = it_cond->GetValue(NEIGHBOUR_CONDITIONS);
                // r_neighb_faces is the vector containing pointers to the three faces around it_cond
                // r_neighb_faces[0] = neighbour face over edge 1-2 of element it_cond
                // r_neighb_faces[1] = neighbour face over edge 2-0 of element it_cond
                // r_neighb_faces[2] = neighbour face over edge 0-1 of element it_cond
                r_neighb_faces(0) = CheckForNeighbourFaces(r_geom[1].Id(), r_geom[2].Id(), r_geom[1].GetValue(NEIGHBOUR_CONDITIONS), it_cond->Id());
                r_neighb_faces(1) = CheckForNeighbourFaces(r_geom[2].Id(), r_geom[0].Id(), r_geom[2].GetValue(NEIGHBOUR_CONDITIONS), it_cond->Id());
                r_neighb_faces(2) = CheckForNeighbourFaces(r_geom[0].Id(), r_geom[1].Id(), r_geom[0].GetValue(NEIGHBOUR_CONDITIONS), it_cond->Id());
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FindConditionsNeighboursProcess::Clear()
{
    block_for_each(mrModelPart.Nodes(), [](Node& rNode) {
        if (rNode.Has(NEIGHBOUR_CONDITIONS)) {
            auto& r_neighbour_conditions = rNode.GetValue(NEIGHBOUR_CONDITIONS);
            r_neighbour_conditions.erase(r_neighbour_conditions.begin(),r_neighbour_conditions.end());
        }
    });
    block_for_each(mrModelPart.Conditions(), [](Condition& rCond) {
        if (rCond.Has(NEIGHBOUR_CONDITIONS)) {
            auto& r_neighbour_conditions = rCond.GetValue(NEIGHBOUR_CONDITIONS);
            r_neighbour_conditions.erase(r_neighbour_conditions.begin(),r_neighbour_conditions.end());
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

void FindConditionsNeighboursProcess::ClearNeighbours()
{
    // Call Clear method
    KRATOS_WARNING("FindConditionsNeighboursProcess") << "ClearNeighbours is a legacy method, please call Clear instead" << std::endl;
    Clear();
}

/***********************************************************************************/
/***********************************************************************************/

Condition::WeakPointer FindConditionsNeighboursProcess::CheckForNeighbourFaces(
    const unsigned int Id1,
    const unsigned int Id2,
    GlobalPointersVector<Condition>& rNeighbourFace,
    const unsigned int Face
    )
{
    // Look for the faces around node Id1
    for( auto it_face =rNeighbourFace.begin(); it_face != rNeighbourFace.end(); it_face++) {
        // Look for the nodes of the neighbour faces
        auto& r_neigh_face_geometry = it_face->GetGeometry();
        for( unsigned int i_node = 0 ; i_node < r_neigh_face_geometry.size(); i_node++) {
            if (r_neigh_face_geometry[i_node].Id() == Id2) {
                if(it_face->Id() != Face) {
                    return *(it_face.base());
                }
            }
        }
    }
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

void FindConditionsNeighboursProcess::ComputeDimension()
{
    // Retrieve from geometry if not defined
    if (mrModelPart.NumberOfConditions() > 0) {
        if (mDim < 0) {
            mDim = mrModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
        }
        if (mDim != 2 && mDim != 3) {
            KRATOS_ERROR << "FindConditionsNeighboursProcess: invalid dimension " << mDim << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters FindConditionsNeighboursProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "help"                         : "This process finds the neighboring conditions for each node and condition in a given model part.",
        "model_part_name"              : "please_provide_model_part_name",
        "average_conditions"           : 10
    })" );
    return default_parameters;
}

}  // namespace Kratos.