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
//

#if !defined(KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"

// Application includes
#include "rans_formulation_process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(RANS_APPLICATION) RansWallDistanceCalculationProcess
: public RansFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansWallDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansWallDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansWallDistanceCalculationProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansWallDistanceCalculationProcess() override = default;

    /// Assignment operator.
    RansWallDistanceCalculationProcess& operator=(RansWallDistanceCalculationProcess const& rOther) = delete;

    /// Copy constructor.
    RansWallDistanceCalculationProcess(RansWallDistanceCalculationProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

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
    double mMaxDistance;

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

    ///@}

}; // Class RansWallDistanceCalculationProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansWallDistanceCalculationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED defined
