//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"

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

/**
 * @class CheckSameModelPartUsingSkinDistanceProcess
 * @ingroup KratosCore
 * @brief Checks that the model part is the same using the skin distance
 * @details Using skin distance, this process checks that the model part is the same
 * @author Vicente Mataix Ferrandiz
 */
template<std::size_t TDim>
class KRATOS_API(KRATOS_CORE) CheckSameModelPartUsingSkinDistanceProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CheckSameModelPartUsingSkinDistanceProcess
    KRATOS_CLASS_POINTER_DEFINITION(CheckSameModelPartUsingSkinDistanceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit CheckSameModelPartUsingSkinDistanceProcess(
        Model& rModel,
        Parameters ThisParameters = Parameters(R"({})")
        )
        : mrSkinModelPart1(rModel.GetModelPart(ThisParameters["skin_model_part_1_name"].GetString())),
          mrSkinModelPart2(rModel.GetModelPart(ThisParameters["skin_model_part_2_name"].GetString())),
          mThisParameters(ThisParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = GetDefaultParameters();
        mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~CheckSameModelPartUsingSkinDistanceProcess() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    /**
     * @brief This method creates an pointer of the process
     * @details We consider as input a Model and a set of Parameters for the sake of generality
     * @warning Must be overrided in each process implementation
     * @param rModel The model to be consider
     * @param ThisParameters The configuration parameters
     */
    Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters
        ) override
    {
        return Kratos::make_shared<CheckSameModelPartUsingSkinDistanceProcess<TDim>>(rModel, ThisParameters);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

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
    std::string Info() const override
    {
        return "CheckSameModelPartUsingSkinDistanceProcess" + std::to_string(TDim) + "D";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CheckSameModelPartUsingSkinDistanceProcess" << TDim << "D";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrSkinModelPart1; /// The modelpart containing the skin of the first model part
    ModelPart& mrSkinModelPart2; /// The modelpart containing the skin of the second model part
    Parameters mThisParameters;  /// The parameters containing the settings

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
    CheckSameModelPartUsingSkinDistanceProcess& operator=(CheckSameModelPartUsingSkinDistanceProcess const& rOther);

    /// Copy constructor.
    //CheckSameModelPartUsingSkinDistanceProcess(CheckSameModelPartUsingSkinDistanceProcess const& rOther);

    ///@}
}; // Class CheckSameModelPartUsingSkinDistanceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
template<std::size_t TDim>
inline std::istream& operator >> (std::istream& rIStream,
                                  CheckSameModelPartUsingSkinDistanceProcess<TDim>& rThis);

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CheckSameModelPartUsingSkinDistanceProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
