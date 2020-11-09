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

protected:
    ///@name Protected Operations
    ///@{

    virtual Process::Pointer GetWallDistanceCalculationProcess(
        ModelPart& rModelPart,
        Parameters LinearSolverParameters,
        const int MaxIterations);

    ///@}

private:
    ///@name Member Variables
    ///@{
    Model& mrModel;
    Parameters mLinearSolverParameters;
    std::string mModelPartName;
    int mMaxIterations;
    int mEchoLevel;
    std::string mWallFlagVariableName;
    bool mWallFlagVariableValue;
    bool mRecalculateAtEachTimeStep;
    bool mCorrectDistancesUsingNeighbors;
    Process::Pointer mpVariationalDistanceCalculationProcess;

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateWallDistances();

    /**
     * @brief Corrects wall distances
     *
     * VariationalDistanceCalculationProcess gives incorrect negative wall distances
     * when mesh is not refined enough. This method iterates through NEIGHBOUR nodes for a node with negative
     * DISTANCE and corrects it by calculated positive DISTANCE average from the NEIGHBOURs
     *
     * Works in MPI as well.
     *
     */
    void CorrectWallDistances();

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
