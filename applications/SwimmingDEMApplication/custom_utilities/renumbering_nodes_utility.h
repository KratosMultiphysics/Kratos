//Author: Miguel Angel Celigueta. maceli@cimne.upc.edu

#if !defined(KRATOS_RENUMBERING_NODES_UTILITY)
#define KRATOS_RENUMBERING_NODES_UTILITY

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/model_part.h"

#include "pybind11/stl.h"


namespace Kratos
{
class RenumberingNodesUtility
{

public:

KRATOS_CLASS_POINTER_DEFINITION(RenumberingNodesUtility);

RenumberingNodesUtility(ModelPart& mp1) {
    mListOfModelParts.push_back(&mp1);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}

RenumberingNodesUtility(ModelPart& mp1, ModelPart& mp2): RenumberingNodesUtility(mp1){
    mListOfModelParts.push_back(&mp2);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}

RenumberingNodesUtility(ModelPart& mp1, ModelPart& mp2, ModelPart& mp3): RenumberingNodesUtility(mp1, mp2) {
    mListOfModelParts.push_back(&mp3);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}

RenumberingNodesUtility(ModelPart& mp1, ModelPart& mp2, ModelPart& mp3, ModelPart& mp4): RenumberingNodesUtility(mp1, mp2, mp3) {
    mListOfModelParts.push_back(&mp4);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}

RenumberingNodesUtility(ModelPart& mp1, ModelPart& mp2, ModelPart& mp3, ModelPart& mp4, ModelPart& mp5): RenumberingNodesUtility(mp1, mp2, mp3, mp4) {
    mListOfModelParts.push_back(&mp5);
    std::map<int,int> aux_map;
    mListOfMapsOfIdsNewToOld.push_back(aux_map);
}


virtual ~RenumberingNodesUtility(){}

void Renumber() {
    int id = 1;
    for (int i=0; i<(int)mListOfModelParts.size(); i++){
        ModelPart& mp = *mListOfModelParts[i];
        std::map<int,int>& new_to_old = mListOfMapsOfIdsNewToOld[i];

        for (int j = 0; j < (int)mp.Nodes().size(); j++){
            auto it = mp.NodesBegin() + j;
            new_to_old[id] = it->Id();
            it->SetId(id);
            id++;
        }
    }
}

void UndoRenumber() {
    for (int i=0; i<(int)mListOfModelParts.size(); i++){
        ModelPart& mp = *mListOfModelParts[i];
        std::map<int,int>& new_to_old = mListOfMapsOfIdsNewToOld[i];

        for (int j = 0; j < (int)mp.Nodes().size(); j++){
            auto it = mp.NodesBegin() + j;
            it->SetId(new_to_old[it->Id()]);
        }
    }
}



private:

std::vector<ModelPart*> mListOfModelParts;
std::vector<std::map<int,int> > mListOfMapsOfIdsNewToOld;

}; // Class RenumberingNodesUtility

} // namespace Kratos.

#endif // KRATOS_RENUMBERING_NODES_UTILITY  defined

