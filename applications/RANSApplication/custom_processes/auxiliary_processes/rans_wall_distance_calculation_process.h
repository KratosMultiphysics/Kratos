//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

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

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class KRATOS_API(RANS_APPLICATION) RansWallDistanceCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansWallDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansWallDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansWallDistanceCalculationProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~RansWallDistanceCalculationProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Model& mrModel;
    Parameters mrParameters;
    std::string mModelPartName;

    typename TLinearSolver::Pointer mpLinearSolver;

    int mMaxIterations;
    int mEchoLevel;
    std::string mWallFlagVariableName;
    bool mWallFlagVariableValue;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateWallDistances();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RansWallDistanceCalculationProcess& operator=(RansWallDistanceCalculationProcess const& rOther);

    /// Copy constructor.
    RansWallDistanceCalculationProcess(RansWallDistanceCalculationProcess const& rOther);

    ///@}

}; // Class RansWallDistanceCalculationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED defined
