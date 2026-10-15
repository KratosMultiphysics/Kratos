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

// System includes
#include <array>
#include <string>
#include <vector>

// External includes

// Project includes
#include "custom_processes/volume_fraction_gsbm_process.h"
#include "custom_utilities/volume_fraction_utilities.h"
#include "iga_application_variables.h"

namespace Kratos
{

void VolumeFractionGsbmProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const std::array<std::string, 4> required_parameters{
        "model_part_name", "geometry_name",
        "elements_model_part_name", "facets_model_part_name"
    };
    for (const std::string& r_key : required_parameters) {
        KRATOS_ERROR_IF(mParameters[r_key].GetString().empty())
            << "VolumeFractionGsbmProcess: \"" << r_key << "\" must be provided." << std::endl;
    }

    // Get the background patch
    const ModelPart& r_model_part = mpModel->GetModelPart(
        mParameters["model_part_name"].GetString()
    );
    const auto& r_background = r_model_part.GetGeometry(
        mParameters["geometry_name"].GetString()
    );
    KRATOS_ERROR_IF(r_background.LocalSpaceDimension() != 2)
        << "VolumeFractionGsbmProcess: only 2D background patches are supported." << std::endl;

    // Get the knot span boundaries, built such that uv = xy
    std::array<std::vector<double>, 2> spans;
    for (IndexType dir = 0; dir < 2; ++dir) {
        r_background.SpansLocalSpace(spans[dir], dir);
    }

    // Get the skins of the domain
    Model empty_model;
    const ModelPart& r_empty_skin = empty_model.CreateModelPart("EmptySkin");
    const std::string outer_skin_name =
        mParameters["skin_model_part_outer_initial_name"].GetString();
    const std::string inner_skin_name =
        mParameters["skin_model_part_inner_initial_name"].GetString();
    const ModelPart& r_outer_skin = outer_skin_name.empty()
        ? r_empty_skin
        : mpModel->GetModelPart(outer_skin_name);
    const ModelPart& r_inner_skin = inner_skin_name.empty()
        ? r_empty_skin
        : mpModel->GetModelPart(inner_skin_name);

    // Get the volume fraction of every knot span
    const Matrix volume_fractions = VolumeFractionUtilities::ComputeKnotSpanVolumeFractions(
        r_outer_skin, r_inner_skin, spans[0], spans[1]
    );

    // Lambda for volume fraction assignment
    const auto assign_volume_fraction = [&volume_fractions](GeometryType& rGeometry) {
        KRATOS_ERROR_IF_NOT(rGeometry.Has(KNOT_SPAN_INDICES))
            << "VolumeFractionGsbmProcess: an entity has no KNOT_SPAN_INDICES. Elements and "
            << "facets must be the ones created by IgaModelerGsbm." << std::endl;
        const auto& r_indices = rGeometry.GetValue(KNOT_SPAN_INDICES);
        rGeometry.SetValue(VOLUME_FRACTION, volume_fractions(r_indices[0], r_indices[1]));
    };

    // Get the elements and assign volume fractions
    ModelPart& r_elements = mpModel->GetModelPart(
        mParameters["elements_model_part_name"].GetString()
    );
    for (auto& r_element : r_elements.Elements()) {
        assign_volume_fraction(r_element.GetGeometry());
    }

    // Get the facets and assign volume fractions
    ModelPart& r_facets = mpModel->GetModelPart(
        mParameters["facets_model_part_name"].GetString()
    );
    for (auto& r_facet : r_facets.Conditions()) {
        for (IndexType side = 0; side < 2; ++side) {
            // Part 0: - side (lower span); part 1: + side (upper span)
            assign_volume_fraction(r_facet.GetGeometry().GetGeometryPart(side));
        }
    }

    // Number of cut knot spans with 0 < alpha < 1
    SizeType number_of_cut_spans = 0;
    for (IndexType i = 0; i < volume_fractions.size1(); ++i) {
        for (IndexType j = 0; j < volume_fractions.size2(); ++j) {
            const double alpha = volume_fractions(i, j);
            number_of_cut_spans += (alpha > 0.0 && alpha < 1.0);
        }
    }

    KRATOS_INFO_IF("VolumeFractionGsbmProcess", mParameters["echo_level"].GetInt() > 0)
        << "Assigned volume fractions to " << r_elements.NumberOfElements()
        << " elements and " << r_facets.NumberOfConditions() << " facets ("
        << number_of_cut_spans << " cut knot spans)." << std::endl;

    KRATOS_CATCH("")
}

const Parameters VolumeFractionGsbmProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "echo_level"                         : 0,
        "model_part_name"                    : "",
        "geometry_name"                      : "",
        "elements_model_part_name"           : "",
        "facets_model_part_name"             : "",
        "skin_model_part_outer_initial_name" : "",
        "skin_model_part_inner_initial_name" : ""
    })");
}

} // namespace Kratos
