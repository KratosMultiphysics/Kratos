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

#if !defined(KRATOS_RANS_NUT_NODAL_UPDATE_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_NODAL_UPDATE_PROCESS_H_INCLUDED

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

/**
 * @brief Updates nodal VISCOSITY
 *
 * This process is used to update nodal historical VISCOSITY based
 * on nodal historical TURBULENT_VISCOSITY and model part's first element's properties's
 * DYNAMIC_VISCOSITY and DENSITY. This is a workaround to support VMS and FractionalStep elements
 *
 * QSVMS formulation does not need this process.
 */

class KRATOS_API(RANS_APPLICATION) RansNutNodalUpdateProcess
: public RansFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansNutNodalUpdateProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutNodalUpdateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutNodalUpdateProcess(
        Model& rModel,
        Parameters rParameters);

    RansNutNodalUpdateProcess(
        Model& rModel,
        const std::string& rModelPartName,
        const int EchoLevel);

    /// Destructor.
    ~RansNutNodalUpdateProcess() override = default;

    /// Assignment operator.
    RansNutNodalUpdateProcess& operator=(RansNutNodalUpdateProcess const& rOther) = delete;

    /// Copy constructor.
    RansNutNodalUpdateProcess(RansNutNodalUpdateProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

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
    int mEchoLevel;
    bool mIsInitialized = false;

    ///@}

}; // Class RansNutNodalUpdateProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansNutNodalUpdateProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_UPDATE_PROCESS_H_INCLUDED defined
