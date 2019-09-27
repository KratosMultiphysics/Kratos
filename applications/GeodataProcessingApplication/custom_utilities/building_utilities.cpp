//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "building_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void BuildingUtilities::CheckOverlap(){

    }


    void BuildingUtilities::SplitBuilding(){

    }


    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const BuildingUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
