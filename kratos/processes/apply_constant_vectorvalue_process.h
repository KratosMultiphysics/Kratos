//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class ApplyConstantVectorValueProcess
 * @brief The base class for all processes in Kratos.
 * @details This function applies a constant Value (and fixity) to all of the nodes in a given mesh
 * @todo Still segfaults if the mesh to which it is applied is not existing
 * @ingroup KratosCore
 * @author Riccardo Rossi
 */
class KRATOS_API(KRATOS_CORE) ApplyConstantVectorValueProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the local flags
    KRATOS_DEFINE_LOCAL_FLAG(X_COMPONENT_FIXED);
    KRATOS_DEFINE_LOCAL_FLAG(Y_COMPONENT_FIXED);
    KRATOS_DEFINE_LOCAL_FLAG(Z_COMPONENT_FIXED);

    /// Pointer definition of ApplyConstantVectorValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantVectorValueProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor for ApplyConstantVectorValueProcess.
     * @details This constructor initializes an instance of the ApplyConstantVectorValueProcess class with the given model part
     * and parameters.
     * @param rModelPart Reference to the model part to which the process is applied.
     * @param parameters Parameters defining the behavior of the process.
     */
    ApplyConstantVectorValueProcess(ModelPart& rModelPart,
                                    Parameters ThisParameters
                                   );

    /**
    * @brief Constructor for ApplyConstantVectorValueProcess.
    * @details This constructor initializes an instance of the ApplyConstantVectorValueProcess class with the given model part,
    * variable, modulus, direction, mesh ID, and options.
    * @param rModelPart Reference to the model part to which the process is applied.
    * @param rVariable Reference to the variable to which the constant vector value is applied.
    * @param Modulus The modulus of the constant vector value.
    * @param Direction The direction of the constant vector value.
    * @param MeshId The ID of the mesh to which the process is applied.
    * @param Options Flags specifying additional options for the process.
    */
    ApplyConstantVectorValueProcess(ModelPart& rModelPart,
                              const Variable< array_1d<double, 3 > >& rVariable,
                              const double Modulus,
                              const Vector& rDirection,
                              std::size_t MeshId,
                              const Flags Options
                              );

    /// Destructor.
    ~ApplyConstantVectorValueProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitialize() override;

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
        return "ApplyConstantVectorValueProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantVectorValueProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected Member Variables
    ///@{

    /// The computed model part
    ModelPart& mrModelPart;

    /// The variable to be fixed
    std::string mVariableName;

    /// The modulus of the vector
    double mModulus = 1.0;

    /// The direction of the vector
    Vector mDirection;

    /// The mesh id (LEGACY, to be removed in the future)
    std::size_t mMeshId = 0;

    ///@}
private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Apply a constant value to a specified variable for nodes in a model part.
     * @details This method applies a constant value to a specified variable for all nodes in a given model part.
     * Optionally, the value can be fixed for the nodes.
     * @param rVariable The variable to which the constant value is applied.
     * @param ToBeFixed Flag indicating whether the variable should be fixed for the nodes.
     * @param Value The constant value to be applied.
     */
    void InternalApplyValue(
        const Variable<double>& rVariable,
        const bool ToBeFixed,
        const double Value
        );

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ApplyConstantVectorValueProcess& operator=(ApplyConstantVectorValueProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantVectorValueProcess(ApplyConstantVectorValueProcess const& rOther);

    ///@}

}; // Class ApplyConstantVectorValueProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantVectorValueProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantVectorValueProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.