//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// Project includes
#include "cad_tessellation_modeler.h"

namespace Kratos
{
    ///@name Stages
    ///@{

    void CadTessellationModeler::SetupModelPart()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name"))
            << "Missing \"cad_model_part\" section" << std::endl;
        ModelPart& cad_model_part =
            mpModel->GetModelPart(mParameters["cad_model_part_name"].GetString());
        KRATOS_ERROR_IF_NOT(mParameters.Has("analysis_model_part_name"))
            << "Missing \"analysis_model_part_name\" section" << std::endl;
        ModelPart& analysis_model_part =
            mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

        const auto& r_geometries = cad_model_part.Geometries();

        for (auto it = r_geometries.begin(); it != r_geometries.end(); ++it)
        {
            if (it->GetGeometryType() == GeometryData::Kratos_Brep_Curve)
            {

                // BrepCurveOnSurface *p_aux = &it;
                // p_aux->pGetCurveOnSurface();
                // auto tessellation = NurbsCurveTessellation<2, PointerVector<Point>>::ComputeTessellation(
                //     p_aux->pGetCurveOnSurface(),
                //     p_aux->PolynomialDegree(0),
                //     p_aux->DomainInterval(),
                //     p_aux->KnotSpanIntervals(),
                //     1e-2);
            }
        }
    }
} // namespace Kratos
