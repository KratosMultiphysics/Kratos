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

#if !defined(KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_UPDATE_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_UPDATE_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"

// Application includes
#include "rans_formulation_process.h"


namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(RANS_APPLICATION) RansNutYPlusWallFunctionUpdateProcess
: public RansFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansNutYPlusWallFunctionUpdateProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutYPlusWallFunctionUpdateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutYPlusWallFunctionUpdateProcess(
        Model& rModel,
        Parameters rParameters);

    RansNutYPlusWallFunctionUpdateProcess(
        Model& rModel,
        const std::string& rModelPartName,
        const double VonKarman,
        const double MinValue,
        const int EchoLevel);

    /// Destructor.
    ~RansNutYPlusWallFunctionUpdateProcess() override = default;

    /// Assignment operator.
    RansNutYPlusWallFunctionUpdateProcess& operator=(RansNutYPlusWallFunctionUpdateProcess const& rOther) = delete;

    /// Copy constructor.
    RansNutYPlusWallFunctionUpdateProcess(RansNutYPlusWallFunctionUpdateProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteAfterCouplingSolveStep() override;

    const Parameters GetDefaultParameters() const override;

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

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;
    std::string mModelPartName;
    double mVonKarman;
    double mMinValue;
    int mEchoLevel;
    bool mIsInitialized = false;

    ///@}

}; // Class RansNutYPlusWallFunctionUpdateProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansNutYPlusWallFunctionUpdateProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_Y_PLUS_WALL_FUNCTION_UPDATE_PROCESS_H_INCLUDED defined
