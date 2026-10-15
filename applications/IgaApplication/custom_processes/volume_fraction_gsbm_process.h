//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) VolumeFractionGsbmProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(VolumeFractionGsbmProcess);

    using SizeType = std::size_t;
    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    VolumeFractionGsbmProcess(
        Model& rModel,
        Parameters ThisParameters)
        : Process()
        , mpModel(&rModel)
        , mParameters(ThisParameters)
    {
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    /// Destructor.
    ~VolumeFractionGsbmProcess() override = default;

    // Get the default parameters
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Operations
    ///@{

    /// Assigns the VOLUME_FRACTION of their knot span to the elements and to
    /// both sides of the facets created by IgaModelerGsbm. Assumes that the
    /// parameter space of the background patch is the physical space (uv = xy)
    void ExecuteInitialize() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "VolumeFractionGsbmProcess";
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;

    ///@}
};

///@}

} // namespace Kratos
