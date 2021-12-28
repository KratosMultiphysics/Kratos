//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Rahul Kikkeri Nagaraja
//

#if !defined(KRATOS_RANS_CHIMERA_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_CHIMERA_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "utilities/binbased_fast_point_locator.h"

// Application includes
#include "rans_formulation_process.h"
#include "custom_processes/apply_rans_chimera_process_monolithic.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

template <int TDim>
class KRATOS_API(RANS_APPLICATION) RansChimeraWallDistanceCalculationProcess
: public RansFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{
    // change typedef to using later. 
    KRATOS_CLASS_POINTER_DEFINITION(RansChimeraWallDistanceCalculationProcess);
    using PointLocatorType = BinBasedFastPointLocator<TDim>;   // Change it to <TDim> and use templates 
    using PointLocatorPointerType = typename PointLocatorType::Pointer;
    using PointLocatorsMapType = std::map<std::string, PointLocatorPointerType>;
    using ChimeraProcessType = ApplyRANSChimeraProcessMonolithic<TDim>;
    using NodeType = ModelPart::NodeType;

    ///@}
    /// Pointer definition of RansChimeraWallDistanceCalculationProcess

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansChimeraWallDistanceCalculationProcess(
        Model& rModel,
        Parameters rParameters,
        ChimeraProcessType& rChimeraProcess);

    /// Destructor.
    ~RansChimeraWallDistanceCalculationProcess() override = default;

    /// Assignment operator.
    RansChimeraWallDistanceCalculationProcess& operator=(RansChimeraWallDistanceCalculationProcess const& rOther) = delete;

    /// Copy constructor.
    RansChimeraWallDistanceCalculationProcess(RansChimeraWallDistanceCalculationProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    // void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;
    std::string mMainModelPartName;
    std::string mWallModelPartName;

    int mMaxLevels;
    int mEchoLevel;
    std::string mDistanceVariableName;
    std::string mNodalAreaVariableName;
    bool mRecalculateAtEachTimeStep;
    bool mIsFormulated;
    double mMaxDistance;
    ChimeraProcessType& mrChimeraProcess;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Calculates wall distances
     *
     * This method first calculates wall distances on nodes adjacent to walls (they do not fall on the walls)
     * based on condition NORMAL. Afterwards, remaining wall adjacent nodes wall distances are calculated based
     * on distributed condition NORMAL on nodes which falls on wall. Then nodes which falls on walls are given negative
     * small wall distance so ParallelLevelSetDistance calculator can calculate wall distances properly. Afterwards,
     * wall nodes' distances are reverted to 0.0.
     *
     * Wall conditions are need to have mWallFlagVariableName Flag with value mWallFlagVariableValue
     * Wall nodes are also need to have mWallFlagVariableName Flag with value mWallFlagVariableValue
     *
     * Nodal NORMAL values are not required, nor them are altered. Only condition NORMAL values are used.
     *
     * Works in MPI as well.
     *
     */
    void CalculateWallDistances();

    void CalculateAnalyticalWallDistances();


    /**
     * @brief Searches for a given node using given locator.
     * @param rBinLocator The bin based locator formulated on the background.
     * This is used to locate nodes on rBoundaryModelPart.
     * @param pNodeToFind The patch node which is to be found.
     * @param[out] prHostElement The element where the node is found.
     * @param[out] rWeights the values of the shape functions at the node inside
     * the elements.
     */
    bool SearchNode(PointLocatorType& rBinLocator,
                    NodeType& rNodeToFind,
                    Element::Pointer& prHostElement,
                    Vector& rWeights);

    ///@}

}; // Class RansChimeraWallDistanceCalculationProcess

///@}
///@name Input and output
///@{

/// output stream function
template<int TDim>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansChimeraWallDistanceCalculationProcess<TDim>& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED defined
