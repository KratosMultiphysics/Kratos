//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "output_gismo_process.h"

namespace Kratos
{

OutputGismoProcess::OutputGismoProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    std::string model_part_name = ThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = rModel.GetModelPart(model_part_name);

    const SizeType thb_surface_id = (r_model_part.ElementsBegin()->GetGeometry().GetGeometryParent(0)).Id();
    auto p_thb_surface = dynamic_pointer_cast<THBSurfaceGeometry<3, PointerVector<Node>>>(r_model_part.GetParentModelPart().pGetGeometry(thb_surface_id)->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));

    gismo::gsTHBSpline<2> thb_geom = p_thb_surface->ToGismo();
    m_thb_geom = thb_geom;
}

void OutputGismoProcess::ExecuteFinalize()
{
    std::string model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);

    const SizeType thb_surface_id = (r_model_part.ElementsBegin()->GetGeometry().GetGeometryParent(0)).Id();
    auto p_thb_surface = dynamic_pointer_cast<THBSurfaceGeometry<3, PointerVector<Node>>>(r_model_part.GetParentModelPart().pGetGeometry(thb_surface_id)->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));

    gismo::gsTHBSpline<2> thb_geom = p_thb_surface->ToGismo();

    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    mp.addPatch(m_thb_geom);
    mp.addAutoBoundaries();

    mp_def.addPatch(thb_geom);

    gsMultiPatch<> deformation = mp_def;

    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    gsField<> solField(mp, deformation);
    gsWriteParaview<>( solField, "solution", 1000, true);
}

const Parameters OutputGismoProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"            : ""
    })" );
    return default_parameters;
}

} // namespace Kratos
