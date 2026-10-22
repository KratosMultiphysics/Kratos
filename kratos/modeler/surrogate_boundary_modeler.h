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

// System includes
#include <vector>

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
 *          SURROGATE_BOUNDARY_PROJECTION. Derived modelers can additionally
 *          request the two face-endpoint projections without duplicating the
 *          closest-point search infrastructure.
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

    /**
     * @brief Updates a surrogate representation for a moved skin.
     * @details This operation is available when "is_moving" is true. The
     *          Cartesian background created by SetupModelPart is retained,
     *          while the active volume view and surrogate faces are rebuilt.
     *          Keeping the background nodes and elements alive preserves their
     *          nodal solution histories across active-set changes.
     */
    void PrepareGeometryModel() override;

    const Parameters GetDefaultParameters() const override;

    std::string Info() const override
    {
        return "SurrogateBoundaryModeler";
    }

protected:
    bool IsMovingModelInitialized() const
    {
        return mMovingModelIsInitialized;
    }

private:
    Model* mpModel = nullptr;
    bool mMovingModelIsInitialized = false;
    std::vector<Node::Pointer> mGridNodes;
    std::vector<Element::Pointer> mBackgroundElements;

    void SetupModelPart2D(const bool IsUpdate);

    /** Update closest-point data while retaining the current active topology. */
};

} // namespace Kratos
