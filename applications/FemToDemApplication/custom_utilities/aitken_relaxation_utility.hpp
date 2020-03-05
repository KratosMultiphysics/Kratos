//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes

// External includes
#include "utilities/math_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}

///@name Type Definitions
///@{
///@}

///@name  Enum's
///@{
///@}

///@name  Functions
///@{
///@}

///@name Kratos Classes
///@{

/** @brief Aitken relaxation technique for FSI PFEM-FEM-DEM coupling
 */
class AitkenRelaxationUtility
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(AitkenRelaxationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

   /**
     * Constructor.
     * AitkenRelaxationUtility
     */
    AitkenRelaxationUtility(const double OmegaOld)
    {
        mOmega_old = OmegaOld;
    }


























} // namespace Kratos