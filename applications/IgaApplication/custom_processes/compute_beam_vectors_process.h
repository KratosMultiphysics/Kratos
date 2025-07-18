//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_COMPUTE_BEAM_VECTORS_PROCESS_H_INCLUDED )
#define  KRATOS_COMPUTE_BEAM_VECTORS_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "geometries/nurbs_curve_geometry.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "iga_application_variables.h"
#include <string>


namespace Kratos
{

///@name Kratos Classes
///@{

/* @class ComputeBeamVectorsProcess
 * @ingroup IgaApplication
 * @brief This process computes the tangential (T0) and normal (N0) vectors for isogeometric beam elements during initialization. */
class KRATOS_API(IGA_APPLICATION) ComputeBeamVectorsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeBeamVectorsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeBeamVectorsProcess);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ComputeBeamVectorsProcess(Model& rModel, Parameters ThisParameters);

    /// Destructor.
    ~ComputeBeamVectorsProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    /// Called once before the solution loop to compute T0 and N0 vectors
    void ExecuteInitialize() override;

    /// Compute T0 and N0 vectors for a given NURBS curve geometry
    void ComputeT0AndN0(
        const NurbsCurveGeometry<3, PointerVector<Node>>& rCurve,
        array_1d<double, 3>& rT0,
        array_1d<double, 3>& rN0);


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ComputeBeamVectorsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeBeamVectorsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@name Member Variables
    ///@{

    /// Model part and different settings
    Model& mrModel;             /// The main model part
    Parameters mThisParameters; /// The parameters (can be used for general pourposes)
    ModelPart& mrThisModelPart; /// The model part to compute
    ///@}
}; // Class ComputeBeamVectorsProcess
}  // namespace Kratos.

#endif // KRATOS_COMPUTE_BEAM_VECTORS_PROCESS_H_INCLUDED  defined