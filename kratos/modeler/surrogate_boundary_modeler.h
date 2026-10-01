//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//


#pragma once

// Project includes
#include "modeler/modeler.h"

namespace Kratos {

/**
 * @class SurrogateBoundaryModeler
 * @brief Creates a structured FEM mesh and its oriented surrogate boundary.
 * @details The initial implementation constructs a 2D Quad4 background mesh
 *          around a closed Line2 skin. Its topology kernel is dimension
 *          independent and is intended to be shared by the future 3D and IGA
 *          adapters. Projection coordinates are stored per face Gauss point in
 *          SURROGATE_BOUNDARY_PROJECTION.
 */
class KRATOS_API(KRATOS_CORE) SurrogateBoundaryModeler : public Modeler
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SurrogateBoundaryModeler);

    SurrogateBoundaryModeler() = default;

    SurrogateBoundaryModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters());

    ~SurrogateBoundaryModeler() override = default;

    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<SurrogateBoundaryModeler>(rModel, ModelParameters);
    }

    void SetupModelPart() override;

    const Parameters GetDefaultParameters() const override;

    std::string Info() const override
    {
        return "SurrogateBoundaryModeler";
    }

private:
    Model* mpModel = nullptr;

    void SetupModelPart2D();
};

} // namespace Kratos
