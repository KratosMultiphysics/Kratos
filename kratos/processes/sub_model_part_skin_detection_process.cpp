//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// Project includes
#include "sub_model_part_skin_detection_process.h"

namespace Kratos
{

template<unsigned int TDim>
SubModelPartSkinDetectionProcess<TDim>::SubModelPartSkinDetectionProcess(
    ModelPart& rModelPart, Parameters Settings)
    : SkinDetectionProcess<TDim>(rModelPart, Settings)
{}

template<> class SubModelPartSkinDetectionProcess<2>;
template<> class SubModelPartSkinDetectionProcess<3>;

}  // namespace Kratos.


