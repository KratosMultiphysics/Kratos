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

#if !defined(KRATOS_RANS_COMPUTE_REACTIONS_PROCESS_H_INCLUDED)
#define KRATOS_RANS_COMPUTE_REACTIONS_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"

// Application includes
#include "rans_point_execution_formulation_process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Clips given scalar variable to a range
 *
 * This process clips a given scalar variable to a range in all nodes in the model part.
 *
 */

class KRATOS_API(RANS_APPLICATION) RansComputeReactionsProcess
: public RansPointExecutionFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = RansPointExecutionFormulationProcess;

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using ConditionType = ModelPart::ConditionType;

    /// Pointer definition of RansComputeReactionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansComputeReactionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansComputeReactionsProcess(
        Model& rModel,
        Parameters rParameters);

    RansComputeReactionsProcess(
        Model& rModel,
        const std::string& rModelPartName,
        const std::vector<std::string>& rExecutionPoints,
        const int EchoLevel);

    /// Destructor.
    ~RansComputeReactionsProcess() override = default;

    /// Assignment operator.
    RansComputeReactionsProcess& operator=(RansComputeReactionsProcess const& rOther) = delete;

    /// Copy constructor.
    RansComputeReactionsProcess(RansComputeReactionsProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

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
    ///@name Private Operations
    ///@{

    void Initialize();

    template<unsigned int TDim>
    void CalculateStrainRate(
        Vector& rStrainRate,
        const GeometryType& rElementGeometry,
        const Matrix& rdNdX) const;

    template<unsigned int TDim>
    void CalculateViscousStressTensorReactionContribution(
        array_1d<double, 3>& rReaction,
        const Vector& rViscousStress,
        const array_1d<double, 3>& rNormal) const;

    template<unsigned int TDim>
    void CalculateReactions(
        ModelPart& rModelPart) const;

    void ExecuteOperation() override;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansComputeReactionsProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_COMPUTE_REACTIONS_PROCESS_H_INCLUDED defined
