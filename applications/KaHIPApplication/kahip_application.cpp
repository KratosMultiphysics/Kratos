//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __ \
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "kahip_application.h"

namespace Kratos {

KratosKaHIPApplication::KratosKaHIPApplication()
    : KratosApplication("KaHIPApplication"),
      mKaHIPPartitioningModeler() {}

void KratosKaHIPApplication::Register()
{
    KRATOS_INFO("") << "     __ __      __  __________  ___                ___            __  _           \n"
                    << "    / //_/___ _/ / / /  _/ __ \\/   |  ____  ____  / (_)________ _/ /_(_)___  ____ \n"
                    << "   / ,< / __ `/ /_/ // // /_/ / /| | / __ \\/ __ \\/ / / ___/ __ `/ __/ / __ \\/ __ \\\n"
                    << "  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /\n"
                    << " /_/ |_\\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\\___/\\__,_/\\__/_/\\____/_/ /_/ \n"
                    << "                                  /_/   /_/                                       \n"
                    << "Initializing KratosKaHIPApplication..." << std::endl;

    KRATOS_REGISTER_MODELER("KaHIPPartitioningModeler", mKaHIPPartitioningModeler);
}

}  // namespace Kratos
