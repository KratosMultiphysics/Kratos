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

#if !defined(KRATOS_FEMDEM_COUPLING_UTILITIES)
#define KRATOS_FEMDEM_COUPLING_UTILITIES

// System includes

// External includes

// Project includes
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

/// The size type definition
typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class FEMDEMCouplingUtilities
 * @ingroup FemToDemApplication
 * @brief This class includes several utilities necessaries for the coupling between the FEM and the DEM
 * @details The methods are static, so it can be called without constructing the class
 * @tparam TDim The dimension of the problem
 * @author Alejandro Cornejo
 */
class FEMDEMCouplingUtilities
{
  public:
    ///@name Type definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(FEMDEMCouplingUtilities);
    /// The index type definition
    typedef std::size_t IndexType;

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Default constructor.
    FEMDEMCouplingUtilities(){}

    void SaveStructuralSolution(ModelPart& rStructureModelPart);

    void InterpolateStructuralSolution(
        ModelPart &rStructureModelPart,
        const double FemDeltaTime,
        const double FemTime,
        const double DemDeltaTime,
        const double DemTime);

    void RestoreStructuralSolution(ModelPart &rStructureModelPart);

    void AddExplicitImpulses(ModelPart &rStructureModelPart, const double DEMTimeStep);

    void ComputeAndTranferAveragedContactTotalForces(ModelPart &rStructureModelPart, const double FEMtimeStep);

    void ResetContactImpulses(ModelPart &rStructureModelPart);
    


}; // class FEMDEMCouplingUtilities
} // namespace Kratos
#endif /* KRATOS_FEMDEM_COUPLING_UTILITIES defined */