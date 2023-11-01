// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"


namespace Kratos {

///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(MED_APPLICATION) MedTestingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MedTestingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(MedTestingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MedTestingUtilities() = delete;

    /// Copy constructor.
    MedTestingUtilities(MedTestingUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MedTestingUtilities& operator=(MedTestingUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static void CheckModelPartsAreEqual(
        const ModelPart& rModelPart1,
        const ModelPart& rModelPart2,
        const bool CheckSubModelParts=true);

    static void AddGeometriesFromElements(
        ModelPart& rModelPart);

    static double ComputeLength(const ModelPart& rModelPart);

    static double ComputeArea(const ModelPart& rModelPart);

    static double ComputeVolume(const ModelPart& rModelPart);

    static double ComputeDomainSize(const ModelPart& rModelPart);

    ///@}

}; // Class MedTestingUtilities

///@}

///@} addtogroup block

} // namespace Kratos
