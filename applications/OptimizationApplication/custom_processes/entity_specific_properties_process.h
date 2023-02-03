//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "processes/process.h"

// Application includes

namespace Kratos {
///@addtogroup OptimizationApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) EntitySpecificPropertiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EntitySpecificPropertiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(EntitySpecificPropertiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    EntitySpecificPropertiesProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~EntitySpecificPropertiesProcess() override = default;

    /// Assignment operator.
    EntitySpecificPropertiesProcess& operator=(EntitySpecificPropertiesProcess const& rOther) = delete;

    /// Copy constructor.
    EntitySpecificPropertiesProcess(EntitySpecificPropertiesProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitializeSolutionStep() override;

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

    bool mIsPropertiesReplaced = false;
    std::vector<const Variable<double>*> mPropertiesSpecificVariablePointersList;


    int mEchoLevel;

    ///@}
    ///@name Private operations
    ///@{

    template<class ContainerType>
    void CreateEntitySpecificProperties(ContainerType& rContainer);

    template<class ContainerType>
    void UpdateEntitySpecificProperties(ContainerType& rContainer);

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const EntitySpecificPropertiesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.
