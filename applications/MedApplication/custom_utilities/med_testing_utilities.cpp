// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "med_testing_utilities.h"


namespace Kratos {

namespace { // helpers namespace

using NodeType = ModelPart::NodeType;
using GeometryType = ModelPart::GeometryType;
using GeometryContainerType = ModelPart::GeometryContainerType;

template<class T>
bool contains(
    const std::vector<T>& rVec,
    const T& rValue)
{
    return std::find(rVec.begin(), rVec.end(), rValue) != rVec.end();
}

// check if the geometries have the same nodes
// for this to be the case, two conditions need to be fulflled:
// 1: geometries have same number of nodes
// 2: nodes of both geometries are in same order and in same coords
bool have_same_nodes(
    const GeometryType& rGeom1,
    const GeometryType& rGeom2)
{
    if (rGeom1.PointsNumber() != rGeom2.PointsNumber()) return false;

    for (std::size_t i=0; i<rGeom1.PointsNumber(); ++i) {
        if (rGeom1[i].Distance(rGeom2[i]) > 1e-12) return false;
    }

    return true;
}

void CheckEntitiesAreEqual(
    const NodeType& rNode1,
    const NodeType& rNode2)
{
    KRATOS_TRY

    // KRATOS_CHECK_EQUAL(rNode1.Id(), rNode2.Id()); // the MEDApp deliberately does not care about IDs

    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.X(),  rNode2.X());
    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.X0(), rNode2.X0());

    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.Y(),  rNode2.Y());
    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.Y0(), rNode2.Y0());

    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.Z(),  rNode2.Z());
    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.Z0(), rNode2.Z0());

    KRATOS_CATCH("")
}

void CheckGeometriesAreEqual(
    const GeometryType& rGeom1,
    const GeometryType& rGeom2)
{
    KRATOS_TRY

    // KRATOS_CHECK_EQUAL(rGeom1.Id(), rGeom2.Id()); // the MEDApp deliberately does not care about IDs

    KRATOS_CHECK_EQUAL(rGeom1.PointsNumber(), rGeom2.PointsNumber());

    // make sure to not accidentially compare base-geometries
    KRATOS_CHECK_NOT_EQUAL(rGeom1.GetGeometryType(), GeometryData::KratosGeometryType::Kratos_generic_type);
    KRATOS_CHECK(GeometryType::IsSame(rGeom1, rGeom2));

    KRATOS_CHECK_DOUBLE_EQUAL(rGeom1.DomainSize(), rGeom2.DomainSize());

    if (rGeom1.PointsNumber() > 1) {
        KRATOS_CHECK_GREATER(rGeom1.DomainSize(), 0.0);
    } else if (rGeom1.PointsNumber() == 1) {
        KRATOS_CHECK_DOUBLE_EQUAL(rGeom1.DomainSize(), 0.0);
    } else {
        KRATOS_ERROR << "Geometry with no points found! " << rGeom1 << std::endl;
    }

    for (std::size_t i=0; i<rGeom1.PointsNumber(); ++i) {
        CheckEntitiesAreEqual(rGeom1[i], rGeom2[i]);
    }

    KRATOS_CATCH("")
}

template<class TContainerType>
void CheckEntitiesAreEqual(
    const TContainerType& rEntities1,
    const TContainerType& rEntities2)
{
    KRATOS_TRY

    // basic checks
    KRATOS_CHECK_EQUAL(rEntities1.size(), rEntities2.size());

    // check entities
    for (std::size_t i=0; i<rEntities1.size(); ++i) {
        CheckEntitiesAreEqual(*(rEntities1.begin()+i), *(rEntities2.begin()+i));
    }

    KRATOS_CATCH("")
}

void CheckGeometriesAreEqual(
    const ModelPart& rModelPart1,
    const ModelPart& rModelPart2)
{
    KRATOS_TRY

    // basic checks
    KRATOS_CHECK_EQUAL(rModelPart1.NumberOfGeometries(), rModelPart2.NumberOfGeometries());

    auto check_geoms = [](const ModelPart& rModelPart1, const ModelPart& rModelPart2){
        for (const auto& r_geom_1 : rModelPart1.Geometries()) {
            for (const auto& r_geom_2 : rModelPart2.Geometries()) {
                if (have_same_nodes(r_geom_1, r_geom_2)) {
                    CheckGeometriesAreEqual(r_geom_1, r_geom_2);
                    goto here;
                }
            }
            KRATOS_ERROR << "no match found for geometry " << r_geom_1 << std::endl;
            here:;
        }
    };

    check_geoms(rModelPart1, rModelPart2);
    check_geoms(rModelPart2, rModelPart1);


    KRATOS_CATCH("")
}

enum class QuantityType
{
    LENGTH = 1,
    AREA = 2,
    VOLUME = 3,
    DOMAIN_SIZE,
};


double ComputeGeometricalQuantity(
    const ModelPart& rModelPart,
    const QuantityType Quantity)
{
    std::function<double(const GeometryType&)> access_function;

    switch (Quantity)
    {
        case QuantityType::LENGTH:
            access_function = [](const auto& rGeom){return rGeom.Length();};
        case QuantityType::AREA:
            access_function = [](const auto& rGeom){return rGeom.Area();};
        case QuantityType::VOLUME:
            access_function = [](const auto& rGeom){return rGeom.Volume();};
        default:
            access_function = [](const auto& rGeom){return rGeom.DomainSize();};
    }

    double total_quantity = 0.0;
    for (const auto& r_geom : rModelPart.Geometries()) {
        if (Quantity == QuantityType::DOMAIN_SIZE ||
            r_geom.LocalSpaceDimension() == static_cast<std::size_t>(Quantity)) {
            const double quantity = access_function(r_geom);
            if (r_geom.LocalSpaceDimension() > 0) {
                KRATOS_CHECK_GREATER(quantity, 0.0);
            }
            total_quantity += quantity;
        }
    }
    return total_quantity;
}


} // helpers namespace

void MedTestingUtilities::CheckModelPartsAreEqual(
    const ModelPart& rModelPart1,
    const ModelPart& rModelPart2,
    const bool CheckSubModelParts)
{
    KRATOS_TRY

    // check nodes
    CheckEntitiesAreEqual(rModelPart1.Nodes(), rModelPart2.Nodes());

    // no need to further check an empty ModelPart
    if (rModelPart1.NumberOfNodes() == 0) {
        return;
    }

    // make sure the Ids of the Nodes start with one
    auto min_id = [](const ModelPart& rModelPart){
        return block_for_each<MinReduction<std::size_t>>(rModelPart.Nodes(), [](const auto& rNode){
            return rNode.Id();
        });
    };

    KRATOS_CHECK_EQUAL(min_id(rModelPart1), 1);
    KRATOS_CHECK_EQUAL(min_id(rModelPart2), 1);

    // check geometries
    CheckGeometriesAreEqual(rModelPart1, rModelPart2);

    if (!CheckSubModelParts) {
        return;
    }

    KRATOS_CHECK_EQUAL(rModelPart1.NumberOfSubModelParts(), rModelPart2.NumberOfSubModelParts());

    const auto& r_smp2_names = rModelPart2.GetSubModelPartNames();

    for (const auto& r_smp_name : rModelPart1.GetSubModelPartNames()) {
        KRATOS_CHECK(contains(r_smp2_names, r_smp_name));
        CheckModelPartsAreEqual(rModelPart1.GetSubModelPart(r_smp_name), rModelPart2.GetSubModelPart(r_smp_name));
    }

    KRATOS_CATCH("")
}


void MedTestingUtilities::AddGeometriesFromElements(
    ModelPart& rModelPart)
{
    KRATOS_TRY

    for (const auto& r_elem : rModelPart.Elements()) {
        rModelPart.AddGeometry(r_elem.pGetGeometry());
    }

    KRATOS_CATCH("")
}


double MedTestingUtilities::ComputeLength(const ModelPart& rModelPart)
{
    return ComputeGeometricalQuantity(rModelPart, QuantityType::LENGTH);
}

double MedTestingUtilities::ComputeArea(const ModelPart& rModelPart)
{
    return ComputeGeometricalQuantity(rModelPart, QuantityType::AREA);
}

double MedTestingUtilities::ComputeVolume(const ModelPart& rModelPart)
{
    return ComputeGeometricalQuantity(rModelPart, QuantityType::VOLUME);
}

double MedTestingUtilities::ComputeDomainSize(const ModelPart& rModelPart)
{
    return ComputeGeometricalQuantity(rModelPart, QuantityType::DOMAIN_SIZE);
}

} // namespace Kratos
