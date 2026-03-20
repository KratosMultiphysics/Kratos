//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . <Kratos> utility for beam section properties
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Siemer
//

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/kratos_components.h"
#include "includes/global_variables.h"
#include "utilities/beam_section_properties_utility.h"

namespace Kratos
{
namespace
{

using SizeType = std::size_t;
using SectionProperties = BeamSectionPropertiesUtility::SectionProperties;

double GetDimension(
    const std::vector<double>& rDimensions,
    const SizeType Index,
    const std::string& rSectionType)
{
    KRATOS_ERROR_IF(Index >= rDimensions.size())
        << "Missing dimension " << Index + 1 << " for CROSS_SECTION_TYPE \""
        << rSectionType << "\"." << std::endl;

    const double value = rDimensions[Index];
    KRATOS_ERROR_IF(value <= 0.0)
        << "Dimension " << Index + 1 << " must be larger than zero for CROSS_SECTION_TYPE \""
        << rSectionType << "\"." << std::endl;

    return value;
}

std::string NormalizeSectionType(const std::string& rSectionType)
{
    std::string normalized_type(rSectionType);
    std::transform(
        normalized_type.begin(),
        normalized_type.end(),
        normalized_type.begin(),
        [](unsigned char c){return static_cast<char>(std::toupper(c));});
    return normalized_type;
}

double ClampShearFactor(const double Value)
{
    return std::max(0.0, std::min(1.0, Value));
}

std::array<double, 2> ApproximateShearFactors(
    const double Area,
    const double ShearResistingAreaY,
    const double ShearResistingAreaZ)
{
    // The NX Nastran element library reference used here does not provide K1/K2
    // for all PBARL section families, so engineering approximations are used below.
    return {
        ClampShearFactor(ShearResistingAreaY / (1.2 * Area)),
        ClampShearFactor(ShearResistingAreaZ / (1.2 * Area))
    };
}

SectionProperties CalculateRod(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double radius = GetDimension(rDimensions, 0, rSectionType);

    SectionProperties section_properties;
    section_properties.Area = Globals::Pi * std::pow(radius, 2);
    section_properties.I22 = Globals::Pi * std::pow(radius, 4) / 4.0;
    section_properties.I33 = section_properties.I22;
    section_properties.TorsionalInertia = section_properties.I22 + section_properties.I33;
    section_properties.ShearFactorY = 0.9;
    section_properties.ShearFactorZ = 0.9;

    return section_properties;
}

SectionProperties CalculateTube(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double outer_radius = GetDimension(rDimensions, 0, rSectionType);
    const double inner_radius = GetDimension(rDimensions, 1, rSectionType);

    KRATOS_ERROR_IF(inner_radius >= outer_radius)
        << "Dimension 2 must be smaller than dimension 1 for CROSS_SECTION_TYPE \""
        << rSectionType << "\"." << std::endl;

    SectionProperties section_properties;
    section_properties.Area = Globals::Pi * (std::pow(outer_radius, 2) - std::pow(inner_radius, 2));
    section_properties.I22 = Globals::Pi * (std::pow(outer_radius, 4) - std::pow(inner_radius, 4)) / 4.0;
    section_properties.I33 = section_properties.I22;
    section_properties.TorsionalInertia = section_properties.I22 + section_properties.I33;
    section_properties.ShearFactorY = 0.5;
    section_properties.ShearFactorZ = 0.5;

    return section_properties;
}

SectionProperties CalculateBar(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double width = GetDimension(rDimensions, 0, rSectionType);
    const double height = GetDimension(rDimensions, 1, rSectionType);

    SectionProperties section_properties;
    section_properties.Area = width * height;
    section_properties.I22 = width * std::pow(height, 3) / 12.0;
    section_properties.I33 = height * std::pow(width, 3) / 12.0;
    section_properties.TorsionalInertia =
        width * std::pow(height, 3) *
        (1.0 / 3.0 - 0.21 * height / width * (1.0 - std::pow(height, 4) / (12.0 * std::pow(width, 4))));
    section_properties.ShearFactorY = 5.0 / 6.0;
    section_properties.ShearFactorZ = 5.0 / 6.0;

    return section_properties;
}

SectionProperties CalculateBox(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double width = GetDimension(rDimensions, 0, rSectionType);
    const double height = GetDimension(rDimensions, 1, rSectionType);
    const double thickness_y = GetDimension(rDimensions, 2, rSectionType);
    const double thickness_z = GetDimension(rDimensions, 3, rSectionType);

    KRATOS_ERROR_IF(2.0 * thickness_y >= height || 2.0 * thickness_z >= width)
        << "Invalid wall thicknesses for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    SectionProperties section_properties;
    section_properties.Area = width * height - (width - 2.0 * thickness_z) * (height - 2.0 * thickness_y);
    section_properties.I22 =
        width * std::pow(height, 3) / 12.0
        - (width - 2.0 * thickness_z) * std::pow(height - 2.0 * thickness_y, 3) / 12.0;
    section_properties.I33 =
        height * std::pow(width, 3) / 12.0
        - (height - 2.0 * thickness_y) * std::pow(width - 2.0 * thickness_z, 3) / 12.0;
    section_properties.TorsionalInertia =
        2.0 * thickness_z * thickness_y * std::pow(width - thickness_z, 2) * std::pow(height - thickness_y, 2)
        / (width * thickness_z + height * thickness_y - std::pow(thickness_z, 2) - std::pow(thickness_y, 2));
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        2.0 * thickness_z * (height - thickness_y),
        2.0 * thickness_y * (width - thickness_z));
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateI(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double l1 = GetDimension(rDimensions, 1, rSectionType);
    const double l2 = GetDimension(rDimensions, 3, rSectionType);
    const double l3 = GetDimension(rDimensions, 2, rSectionType);
    const double h1 = GetDimension(rDimensions, 4, rSectionType);
    const double h3 = GetDimension(rDimensions, 5, rSectionType);
    const double total_height = GetDimension(rDimensions, 0, rSectionType);
    const double h2 = total_height - h1 - h3;

    KRATOS_ERROR_IF(h2 <= 0.0)
        << "Invalid flange/web dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = h1 / 2.0;
    const double y2 = h2 / 2.0 + h1;
    const double y3 = total_height - h3 / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2 + l3 * h3;

    const double y_bar = (l1 * h1 * y1 + l2 * h2 * y2 + l3 * h3 * y3) / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2)
        + l3 * std::pow(h3, 3) / 12.0 + l3 * h3 * std::pow(y_bar - y3, 2);
    section_properties.I33 = h1 * std::pow(l1, 3) / 12.0 + h2 * std::pow(l2, 3) / 12.0 + h3 * std::pow(l3, 3) / 12.0;
    section_properties.TorsionalInertia =
        (std::pow(h3, 3) * l3 + std::pow(h1, 3) * l1 + std::pow(l2, 3) * (total_height - (h1 + h3) / 2.0)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l2 * h2,
        l1 * h1 + l3 * h3);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateChan(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double l1 = GetDimension(rDimensions, 0, rSectionType);
    const double total_height = GetDimension(rDimensions, 1, rSectionType);
    const double l2 = GetDimension(rDimensions, 2, rSectionType);
    const double wall_thickness = GetDimension(rDimensions, 3, rSectionType);
    const double l3 = l1;
    const double h1 = wall_thickness;
    const double h2 = total_height - 2.0 * wall_thickness;
    const double h3 = wall_thickness;

    KRATOS_ERROR_IF(h2 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = wall_thickness / 2.0;
    const double y2 = total_height / 2.0;
    const double y3 = total_height - wall_thickness / 2.0;
    const double z1 = l1 / 2.0;
    const double z2 = l2 / 2.0;
    const double z3 = l3 / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2 + l3 * h3;

    const double y_bar = (l1 * h1 * y1 + l2 * h2 * y2 + l3 * h3 * y3) / section_properties.Area;
    const double z_bar = (l1 * h1 * z1 + l2 * h2 * z2 + l3 * h3 * z3) / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2)
        + l3 * std::pow(h3, 3) / 12.0 + l3 * h3 * std::pow(y_bar - y3, 2);
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2)
        + h3 * std::pow(l3, 3) / 12.0 + h3 * l3 * std::pow(z_bar - z3, 2);
    section_properties.TorsionalInertia =
        (2.0 * (l1 - l2 / 2.0) * std::pow(wall_thickness, 3) + (total_height - wall_thickness) * std::pow(l2, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l2 * h2,
        l1 * h1 + l3 * h3);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateT(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double flange_width = GetDimension(rDimensions, 0, rSectionType);
    const double total_height = GetDimension(rDimensions, 1, rSectionType);
    const double flange_thickness = GetDimension(rDimensions, 2, rSectionType);
    const double web_thickness = GetDimension(rDimensions, 3, rSectionType);
    const double h1 = total_height - flange_thickness;
    const double h2 = flange_thickness;
    const double y1 = h1 / 2.0;
    const double y2 = total_height - flange_thickness / 2.0;

    KRATOS_ERROR_IF(h1 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    SectionProperties section_properties;
    section_properties.Area = web_thickness * h1 + flange_width * h2;
    const double y_bar = (web_thickness * h1 * y1 + flange_width * h2 * y2) / section_properties.Area;

    section_properties.I22 =
        web_thickness * std::pow(h1, 3) / 12.0 + web_thickness * h1 * std::pow(y_bar - y1, 2)
        + flange_width * std::pow(h2, 3) / 12.0 + flange_width * h2 * std::pow(y_bar - y2, 2);
    section_properties.I33 = h1 * std::pow(web_thickness, 3) / 12.0 + h2 * std::pow(flange_width, 3) / 12.0;
    section_properties.TorsionalInertia =
        (std::pow(flange_thickness, 3) * flange_width + std::pow(web_thickness, 3) * (total_height - flange_thickness / 2.0)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        web_thickness * h1,
        flange_width * h2);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateCross(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double vertical_height = GetDimension(rDimensions, 0, rSectionType);
    const double horizontal_width = GetDimension(rDimensions, 1, rSectionType);
    const double vertical_thickness = GetDimension(rDimensions, 2, rSectionType);
    const double horizontal_thickness = GetDimension(rDimensions, 3, rSectionType);

    SectionProperties section_properties;
    section_properties.Area = vertical_thickness * horizontal_width + horizontal_thickness * vertical_height;
    section_properties.I22 =
        horizontal_width * std::pow(vertical_thickness, 3) / 12.0
        + vertical_height * std::pow(horizontal_thickness, 3) / 12.0;
    section_properties.I33 =
        vertical_thickness * std::pow(horizontal_width, 3) / 12.0
        + 2.0 * (horizontal_thickness * std::pow(vertical_height / 2.0, 3) / 12.0
        + horizontal_thickness * vertical_height / 2.0 * std::pow(horizontal_width / 2.0 + vertical_height / 4.0, 2));
    section_properties.TorsionalInertia =
        (vertical_height * std::pow(horizontal_thickness, 3) + vertical_thickness * horizontal_width) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        horizontal_thickness * vertical_height,
        vertical_thickness * horizontal_width);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateH(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double vertical_height = GetDimension(rDimensions, 0, rSectionType);
    const double flange_width = GetDimension(rDimensions, 1, rSectionType);
    const double flange_thickness = GetDimension(rDimensions, 2, rSectionType);
    const double web_thickness = GetDimension(rDimensions, 3, rSectionType);

    SectionProperties section_properties;
    section_properties.Area = flange_width * flange_thickness + web_thickness * vertical_height;
    section_properties.I22 =
        vertical_height * std::pow(web_thickness, 3) / 12.0
        + flange_width * std::pow(flange_thickness, 3) / 12.0;
    section_properties.I33 =
        web_thickness * std::pow(vertical_height, 3) / 12.0
        + 2.0 * (flange_thickness * std::pow(flange_width / 2.0, 3) / 12.0
        + flange_thickness * flange_width / 2.0 * std::pow(vertical_height / 2.0 + flange_width / 4.0, 2));
    section_properties.TorsionalInertia =
        (2.0 * flange_thickness * std::pow(flange_width / 2.0, 3) + vertical_height * std::pow(web_thickness, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        web_thickness * vertical_height,
        flange_width * flange_thickness);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateT1(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double h2 = GetDimension(rDimensions, 0, rSectionType);
    const double l1 = GetDimension(rDimensions, 1, rSectionType);
    const double l2 = GetDimension(rDimensions, 2, rSectionType);
    const double h1 = GetDimension(rDimensions, 3, rSectionType);
    const double z1 = l1 / 2.0;
    const double z2 = l1 + l2 / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2;
    const double z_bar = (l1 * h1 * z1 + l2 * h2 * z2) / section_properties.Area;

    section_properties.I22 = l1 * std::pow(h1, 3) / 12.0 + l2 * std::pow(h2, 3) / 12.0;
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2);
    section_properties.TorsionalInertia = (h1 * std::pow(l1, 3) + l1 * std::pow(h1, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l2 * h2,
        l1 * h1);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateI1(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double l2 = GetDimension(rDimensions, 0, rSectionType);
    const double l1 = GetDimension(rDimensions, 1, rSectionType);
    const double middle_gap = GetDimension(rDimensions, 2, rSectionType);
    const double h1 = GetDimension(rDimensions, 3, rSectionType);
    const double h2 = (h1 - middle_gap) / 2.0;
    const double y1 = (h1 - middle_gap) / 4.0;
    const double z1 = l2 / 4.0 + l1 / 2.0;

    KRATOS_ERROR_IF(h2 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + 2.0 * l2 * h2;
    const double y_bar = h1 / 2.0;

    section_properties.I22 = l1 * std::pow(h1, 3) / 12.0 + 2.0 * (l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y1, 2));
    section_properties.I33 = h1 * std::pow(l1, 3) / 12.0 + 4.0 * (h2 * std::pow(l2 / 2.0, 3) / 12.0 + h2 * (l2 / 2.0) * std::pow(z1, 2));
    section_properties.TorsionalInertia = (2.0 * l1 * std::pow(h1, 3) + h1 * std::pow(l1, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l1 * h1,
        2.0 * l2 * h2);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateChan1(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double l2 = GetDimension(rDimensions, 0, rSectionType);
    const double l1 = GetDimension(rDimensions, 1, rSectionType);
    const double middle_gap = GetDimension(rDimensions, 2, rSectionType);
    const double total_height = GetDimension(rDimensions, 3, rSectionType);
    const double l3 = l2;
    const double h1 = total_height;
    const double h2 = total_height - middle_gap / 2.0;
    const double h3 = (total_height - middle_gap) / 2.0;

    KRATOS_ERROR_IF(h3 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = total_height / 2.0;
    const double y2 = total_height - (total_height - middle_gap) / 4.0;
    const double y3 = (total_height - middle_gap) / 4.0;
    const double z1 = l1 / 2.0;
    const double z2 = l1 + l2 / 2.0;
    const double z3 = l1 + l2 / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2 + l3 * h3;
    const double y_bar = (l1 * h1 * y1 + l2 * h2 * y2 + l3 * h3 * y3) / section_properties.Area;
    const double z_bar = (l1 * h1 * z1 + l2 * h2 * z2 + l3 * h3 * z3) / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2)
        + l3 * std::pow(h3, 3) / 12.0 + l3 * h3 * std::pow(y_bar - y3, 2);
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2)
        + h3 * std::pow(l3, 3) / 12.0 + h3 * l3 * std::pow(z_bar - z3, 2);
    section_properties.TorsionalInertia = (2.0 * l2 * std::pow(h2, 3) + h1 * std::pow(l1, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l1 * h1,
        l2 * h2 + l3 * h3);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateZ(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double l1 = GetDimension(rDimensions, 0, rSectionType);
    const double l2 = GetDimension(rDimensions, 1, rSectionType);
    const double middle_gap = GetDimension(rDimensions, 2, rSectionType);
    const double total_height = GetDimension(rDimensions, 3, rSectionType);
    const double l3 = l1;
    const double h1 = (total_height - middle_gap) / 2.0;
    const double h2 = total_height;
    const double h3 = (total_height - middle_gap) / 2.0;

    KRATOS_ERROR_IF(h1 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = total_height - (total_height - middle_gap) / 4.0;
    const double y2 = total_height / 2.0;
    const double y3 = (total_height - middle_gap) / 4.0;
    const double z1 = -l1 / 2.0;
    const double z2 = l2 / 2.0;
    const double z3 = l2 + l1 / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2 + l3 * h3;
    const double y_bar = (l1 * h1 * y1 + l2 * h2 * y2 + l3 * h3 * y3) / section_properties.Area;
    const double z_bar = (l1 * h1 * z1 + l2 * h2 * z2 + l3 * h3 * z3) / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2)
        + l3 * std::pow(h3, 3) / 12.0 + l3 * h3 * std::pow(y_bar - y3, 2);
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2)
        + h3 * std::pow(l3, 3) / 12.0 + h3 * l3 * std::pow(z_bar - z3, 2);
    section_properties.TorsionalInertia = (2.0 * l1 * std::pow(h1, 3) + h2 * std::pow(l2, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l2 * h2,
        l1 * h1 + l3 * h3);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateChan2(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double l2 = GetDimension(rDimensions, 0, rSectionType);
    const double h1 = GetDimension(rDimensions, 1, rSectionType);
    const double total_height = GetDimension(rDimensions, 2, rSectionType);
    const double l1 = GetDimension(rDimensions, 3, rSectionType);
    const double l3 = l2;
    const double h2 = total_height - h1;
    const double h3 = total_height - h1;

    KRATOS_ERROR_IF(h2 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = h1 / 2.0;
    const double y2 = h1 + (total_height - h1) / 2.0;
    const double y3 = h1 + (total_height - h1) / 2.0;
    const double z1 = l1 / 2.0;
    const double z2 = l2 / 2.0;
    const double z3 = l1 - l2 / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2 + l3 * h3;
    const double y_bar = (l1 * h1 * y1 + l2 * h2 * y2 + l3 * h3 * y3) / section_properties.Area;
    const double z_bar = (l1 * h1 * z1 + l2 * h2 * z2 + l3 * h3 * z3) / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2)
        + l3 * std::pow(h3, 3) / 12.0 + l3 * h3 * std::pow(y_bar - y3, 2);
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2)
        + h3 * std::pow(l3, 3) / 12.0 + h3 * l3 * std::pow(z_bar - z3, 2);
    section_properties.TorsionalInertia = (2.0 * h1 * std::pow(l1, 3) + l1 * std::pow(h1, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l2 * h2 + l3 * h3,
        l1 * h1);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateT2(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double l1 = GetDimension(rDimensions, 0, rSectionType);
    const double total_height = GetDimension(rDimensions, 1, rSectionType);
    const double h1 = GetDimension(rDimensions, 2, rSectionType);
    const double l2 = GetDimension(rDimensions, 3, rSectionType);
    const double h2 = total_height - h1;

    KRATOS_ERROR_IF(h2 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = h1 / 2.0;
    const double y2 = h1 + h2 / 2.0;
    const double z1 = l1 / 2.0;
    const double z2 = l1 / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2;
    const double y_bar = (l1 * h1 * y1 + l2 * h2 * y2) / section_properties.Area;
    const double z_bar = (l1 * h1 * z1 + l2 * h2 * z2) / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2);
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2);
    section_properties.TorsionalInertia = (l1 * std::pow(h1, 3) + h2 * std::pow(l2, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l2 * h2,
        l1 * h1);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateBox1(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double l1 = GetDimension(rDimensions, 0, rSectionType);
    const double total_height = GetDimension(rDimensions, 1, rSectionType);
    const double h4 = GetDimension(rDimensions, 2, rSectionType);
    const double h1 = GetDimension(rDimensions, 3, rSectionType);
    const double l3 = GetDimension(rDimensions, 4, rSectionType);
    const double l2 = GetDimension(rDimensions, 5, rSectionType);
    const double l4 = l1;
    const double h2 = total_height - h4 - h1;
    const double h3 = h2;

    KRATOS_ERROR_IF(h2 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = h1 / 2.0;
    const double y2 = h1 + h2 / 2.0;
    const double y3 = h1 + h3 / 2.0;
    const double y4 = total_height - h4 / 2.0;
    const double z1 = l1 / 2.0;
    const double z2 = l2 / 2.0;
    const double z3 = l1 - l3 / 2.0;
    const double z4 = l1 / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2 + l3 * h3 + l4 * h4;
    const double y_bar = (l1 * h1 * y1 + l2 * h2 * y2 + l3 * h3 * y3 + l4 * h4 * y4) / section_properties.Area;
    const double z_bar = (l1 * h1 * z1 + l2 * h2 * z2 + l3 * h3 * z3 + l4 * h4 * z4) / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2)
        + l3 * std::pow(h3, 3) / 12.0 + l3 * h3 * std::pow(y_bar - y3, 2)
        + l4 * std::pow(h4, 3) / 12.0 + l4 * h4 * std::pow(y_bar - y4, 2);
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2)
        + h3 * std::pow(l3, 3) / 12.0 + h3 * l3 * std::pow(z_bar - z3, 2)
        + h4 * std::pow(l4, 3) / 12.0 + h4 * l4 * std::pow(z_bar - z4, 2);

    const double tl = h4;
    const double tr = h1;
    const double tt = l3;
    const double tb = l2;
    const double area_midline = (l1 - (tl + tr) / 2.0) * (total_height - (tt + tb) / 2.0);
    const double sum_l_over_t =
        (total_height - (tt + tb) / 2.0) / tl
        + (total_height - (tt + tb) / 2.0) / tr
        + (l1 - (tl + tr) / 2.0) / tt
        + (l1 - (tl + tr) / 2.0) / tb;
    section_properties.TorsionalInertia = 4.0 * std::pow(area_midline, 2) / sum_l_over_t;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l2 * h2 + l3 * h3,
        l1 * h1 + l4 * h4);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateHexa(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double side = GetDimension(rDimensions, 0, rSectionType);
    const double width = GetDimension(rDimensions, 1, rSectionType);
    const double thickness = GetDimension(rDimensions, 2, rSectionType);

    SectionProperties section_properties;
    section_properties.Area = (width - 2.0 * side) * thickness + side * thickness;
    section_properties.I22 =
        (width - 2.0 * side) * std::pow(thickness, 3) / 12.0
        + 4.0 * (thickness * side / 4.0 / 6.0 * std::pow(thickness / 2.0, 2));
    section_properties.I33 =
        thickness * std::pow(width - 2.0 * side, 3) / 12.0
        + 2.0 * (thickness * std::pow(side, 3) / 36.0 + thickness * side / 2.0 * std::pow(width / 2.0 - 2.0 * side / 3.0, 2));
    section_properties.TorsionalInertia = 0.1154 * std::pow(width, 4);
    section_properties.ShearFactorY = 5.0 / 6.0;
    section_properties.ShearFactorZ = 5.0 / 6.0;

    return section_properties;
}

SectionProperties CalculateHat(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double total_height = GetDimension(rDimensions, 0, rSectionType);
    const double wall_thickness = GetDimension(rDimensions, 1, rSectionType);
    const double middle_width = GetDimension(rDimensions, 2, rSectionType);
    const double lip_width = GetDimension(rDimensions, 3, rSectionType);

    const double l1 = lip_width;
    const double l2 = wall_thickness;
    const double l3 = middle_width - 2.0 * wall_thickness;
    const double l4 = wall_thickness;
    const double l5 = lip_width;
    const double h1 = wall_thickness;
    const double h2 = total_height;
    const double h3 = wall_thickness;
    const double h4 = total_height;
    const double h5 = wall_thickness;

    KRATOS_ERROR_IF(l3 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = wall_thickness / 2.0;
    const double y2 = total_height / 2.0;
    const double y3 = total_height - wall_thickness / 2.0;
    const double y4 = total_height / 2.0;
    const double y5 = wall_thickness / 2.0;
    const double z1 = lip_width / 2.0;
    const double z2 = lip_width + wall_thickness / 2.0;
    const double z3 = lip_width + middle_width / 2.0;
    const double z4 = lip_width + middle_width - wall_thickness / 2.0;
    const double z5 = lip_width + middle_width + lip_width / 2.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2 + l3 * h3 + l4 * h4 + l5 * h5;
    const double y_bar = (l1 * h1 * y1 + l2 * h2 * y2 + l3 * h3 * y3 + l4 * h4 * y4 + l5 * h5 * y5) / section_properties.Area;
    const double z_bar = (l1 * h1 * z1 + l2 * h2 * z2 + l3 * h3 * z3 + l4 * h4 * z4 + l5 * h5 * z5) / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2)
        + l3 * std::pow(h3, 3) / 12.0 + l3 * h3 * std::pow(y_bar - y3, 2)
        + l4 * std::pow(h4, 3) / 12.0 + l4 * h4 * std::pow(y_bar - y4, 2)
        + l5 * std::pow(h5, 3) / 12.0 + l5 * h5 * std::pow(y_bar - y5, 2);
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2)
        + h3 * std::pow(l3, 3) / 12.0 + h3 * l3 * std::pow(z_bar - z3, 2)
        + h4 * std::pow(l4, 3) / 12.0 + h4 * l4 * std::pow(z_bar - z4, 2)
        + h5 * std::pow(l5, 3) / 12.0 + h5 * l5 * std::pow(z_bar - z5, 2);
    section_properties.TorsionalInertia =
        (2.0 * lip_width * std::pow(wall_thickness, 3)
        + 2.0 * (total_height - 2.0 * wall_thickness) * std::pow(wall_thickness, 3)
        + middle_width * std::pow(wall_thickness, 3)) / 3.0;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l2 * h2 + l4 * h4,
        l1 * h1 + l3 * h3 + l5 * h5);
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

SectionProperties CalculateHat1(const std::vector<double>& rDimensions, const std::string& rSectionType)
{
    const double width = GetDimension(rDimensions, 0, rSectionType);
    const double total_height = GetDimension(rDimensions, 1, rSectionType);
    const double top_width = GetDimension(rDimensions, 2, rSectionType);
    const double wall_thickness = GetDimension(rDimensions, 3, rSectionType);
    const double bottom_thickness = GetDimension(rDimensions, 4, rSectionType);

    const double l1 = width;
    const double l2 = (width - top_width) / 2.0;
    const double l3 = wall_thickness;
    const double l4 = top_width;
    const double l5 = wall_thickness;
    const double l6 = (width - top_width) / 2.0;
    const double h1 = bottom_thickness;
    const double h2 = wall_thickness;
    const double h3 = total_height - bottom_thickness - wall_thickness;
    const double h4 = top_width;
    const double h5 = total_height - bottom_thickness - wall_thickness;
    const double h6 = wall_thickness;

    KRATOS_ERROR_IF(l2 <= 0.0 || h3 <= 0.0)
        << "Invalid dimensions for CROSS_SECTION_TYPE \"" << rSectionType << "\"." << std::endl;

    const double y1 = bottom_thickness / 2.0;
    const double y2 = bottom_thickness + wall_thickness / 2.0;
    const double y3 = bottom_thickness + h3 / 2.0;
    const double y4 = total_height - wall_thickness / 2.0;
    const double y5 = bottom_thickness + h3 / 2.0;
    const double y6 = bottom_thickness + wall_thickness / 2.0;
    const double z1 = width / 2.0;
    const double z2 = (width - top_width) / 4.0;
    const double z3 = (width - top_width + wall_thickness) / 2.0;
    const double z4 = width / 2.0;
    const double z5 = (width + top_width - wall_thickness) / 2.0;
    const double z6 = width - (width - top_width) / 4.0;

    SectionProperties section_properties;
    section_properties.Area = l1 * h1 + l2 * h2 + l3 * h3 + l4 * h4 + l5 * h5 + l6 * h6;
    const double y_bar =
        (l1 * h1 * y1 + l2 * h2 * y2 + l3 * h3 * y3 + l4 * h4 * y4 + l5 * h5 * y5 + l6 * h6 * y6)
        / section_properties.Area;
    const double z_bar =
        (l1 * h1 * z1 + l2 * h2 * z2 + l3 * h3 * z3 + l4 * h4 * z4 + l5 * h5 * z5 + l6 * h6 * z6)
        / section_properties.Area;

    section_properties.I22 =
        l1 * std::pow(h1, 3) / 12.0 + l1 * h1 * std::pow(y_bar - y1, 2)
        + l2 * std::pow(h2, 3) / 12.0 + l2 * h2 * std::pow(y_bar - y2, 2)
        + l3 * std::pow(h3, 3) / 12.0 + l3 * h3 * std::pow(y_bar - y3, 2)
        + l4 * std::pow(h4, 3) / 12.0 + l4 * h4 * std::pow(y_bar - y4, 2)
        + l5 * std::pow(h5, 3) / 12.0 + l5 * h5 * std::pow(y_bar - y5, 2)
        + l6 * std::pow(h6, 3) / 12.0 + l6 * h6 * std::pow(y_bar - y6, 2);
    section_properties.I33 =
        h1 * std::pow(l1, 3) / 12.0 + h1 * l1 * std::pow(z_bar - z1, 2)
        + h2 * std::pow(l2, 3) / 12.0 + h2 * l2 * std::pow(z_bar - z2, 2)
        + h3 * std::pow(l3, 3) / 12.0 + h3 * l3 * std::pow(z_bar - z3, 2)
        + h4 * std::pow(l4, 3) / 12.0 + h4 * l4 * std::pow(z_bar - z4, 2)
        + h5 * std::pow(l5, 3) / 12.0 + h5 * l5 * std::pow(z_bar - z5, 2)
        + h6 * std::pow(l6, 3) / 12.0 + h6 * l6 * std::pow(z_bar - z6, 2);

    const double a =
        4.0 * bottom_thickness * wall_thickness * std::pow(top_width - wall_thickness, 2)
        * std::pow(total_height - wall_thickness / 2.0 - bottom_thickness / 2.0, 2);
    const double b =
        (top_width - wall_thickness) * (wall_thickness + bottom_thickness)
        + bottom_thickness * (2.0 * total_height - wall_thickness - bottom_thickness);
    const double c = (width - top_width) * std::pow(wall_thickness + bottom_thickness, 3) / 3.0;
    section_properties.TorsionalInertia = a / b + c;
    const auto shear_factors = ApproximateShearFactors(
        section_properties.Area,
        l3 * h3 + l5 * h5,
        section_properties.Area - (l3 * h3 + l5 * h5));
    section_properties.ShearFactorY = shear_factors[0];
    section_properties.ShearFactorZ = shear_factors[1];

    return section_properties;
}

} // namespace

BeamSectionPropertiesUtility::SectionProperties BeamSectionPropertiesUtility::CalculateProperties(
    const std::string& rSectionType,
    const std::vector<double>& rDimensions)
{
    const std::string section_type = NormalizeSectionType(rSectionType);

    if (section_type == "ROD") {
        return CalculateRod(rDimensions, section_type);
    } else if (section_type == "TUBE") {
        return CalculateTube(rDimensions, section_type);
    } else if (section_type == "BAR") {
        return CalculateBar(rDimensions, section_type);
    } else if (section_type == "BOX") {
        return CalculateBox(rDimensions, section_type);
    } else if (section_type == "I") {
        return CalculateI(rDimensions, section_type);
    } else if (section_type == "CHAN") {
        return CalculateChan(rDimensions, section_type);
    } else if (section_type == "T") {
        return CalculateT(rDimensions, section_type);
    } else if (section_type == "CROSS") {
        return CalculateCross(rDimensions, section_type);
    } else if (section_type == "H") {
        return CalculateH(rDimensions, section_type);
    } else if (section_type == "T1") {
        return CalculateT1(rDimensions, section_type);
    } else if (section_type == "I1") {
        return CalculateI1(rDimensions, section_type);
    } else if (section_type == "CHAN1") {
        return CalculateChan1(rDimensions, section_type);
    } else if (section_type == "Z") {
        return CalculateZ(rDimensions, section_type);
    } else if (section_type == "CHAN2") {
        return CalculateChan2(rDimensions, section_type);
    } else if (section_type == "T2") {
        return CalculateT2(rDimensions, section_type);
    } else if (section_type == "BOX1") {
        return CalculateBox1(rDimensions, section_type);
    } else if (section_type == "HEXA") {
        return CalculateHexa(rDimensions, section_type);
    } else if (section_type == "HAT") {
        return CalculateHat(rDimensions, section_type);
    } else if (section_type == "HAT1") {
        return CalculateHat1(rDimensions, section_type);
    }

    KRATOS_ERROR << "Unsupported CROSS_SECTION_TYPE \"" << section_type << "\"." << std::endl;
}

} // namespace Kratos
