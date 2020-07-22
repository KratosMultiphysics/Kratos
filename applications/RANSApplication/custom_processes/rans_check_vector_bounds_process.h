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

#if !defined(KRATOS_RANS_CHECK_VECTOR_BOUNDS_PROCESS_H_INCLUDED)
#define KRATOS_RANS_CHECK_VECTOR_BOUNDS_PROCESS_H_INCLUDED

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

///@name Kratos Classes
///@{

/**
 * @brief Checks lower and upper bounds of a vector
 *
 * This process checks lower and upper bounds of a vector variable in nodes in a given modelpart
 *
 */

class KRATOS_API(RANS_APPLICATION) RansCheckVectorBoundsProcess : public Process
{
    enum VectorComponent { Magnitude, X, Y, Z };

public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    /// Pointer definition of RansCheckVectorBoundsProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansCheckVectorBoundsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansCheckVectorBoundsProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansCheckVectorBoundsProcess() override = default;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void Execute() override;

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
    Parameters mrParameters;
    std::string mModelPartName;
    std::string mVariableName;

    VectorComponent mVectorComponent;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RansCheckVectorBoundsProcess& operator=(RansCheckVectorBoundsProcess const& rOther);

    /// Copy constructor.
    RansCheckVectorBoundsProcess(RansCheckVectorBoundsProcess const& rOther);

    ///@}

}; // Class RansCheckVectorBoundsProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansCheckVectorBoundsProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_CHECK_VECTOR_BOUNDS_PROCESS_H_INCLUDED defined
