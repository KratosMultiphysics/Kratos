//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Ricky Aristio
//

// System includes

// External includes

// Project includes
#include "modified_shape_functions/triangle_3d_3_modified_shape_functions.h"
#include "utilities/math_utils.h"

namespace Kratos
{

/// Triangle3D3ModifiedShapeFunctions implementation
/// Default constructor
Triangle3D3ModifiedShapeFunctions::Triangle3D3ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    ModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTriangleSplitter(Kratos::make_shared<DivideTriangle3D3<Node>>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTriangleSplitter->GenerateDivision();
    mpTriangleSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Triangle3D3ModifiedShapeFunctions::~Triangle3D3ModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Triangle3D3ModifiedShapeFunctions::Info() const {
    return "Triangle2D3N modified shape functions computation class.";
};

/// Print information about this object.
void Triangle3D3ModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Triangle2D3N modified shape functions computation class.";
};

/// Print object's data.
void Triangle3D3ModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Triangle2D3N modified shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    std::ostringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};

// Returns a pointer to the splitting utility
const DivideGeometry<Node>::Pointer Triangle3D3ModifiedShapeFunctions::pGetSplittingUtil() const {
    return mpTriangleSplitter;
};

void Triangle3D3ModifiedShapeFunctions::SetCondensationMatrix(Matrix& rIntPointCondMatrix)
{
    ModifiedShapeFunctions::SetCondensationMatrix(
        rIntPointCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

void Triangle3D3ModifiedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    ModifiedShapeFunctions::SetCondensationMatrix(
        rPosSideCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

void Triangle3D3ModifiedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    ModifiedShapeFunctions::SetCondensationMatrix(
        rNegSideCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

void Triangle3D3ModifiedShapeFunctions::ComputeFaceNormalOnOneSide(
    AreaNormalsContainerType& rInterfaceAreaNormalValues,
    const std::vector<IndexedPointGeometryPointerType>& rInterfacesVector,
    const IntegrationMethodType IntegrationMethod)
{
    // The base ModifiedShapeFunctions::ComputeFaceNormalOnOneSide calls
    // IndexedPointGeometryType::Normal(), which for a genuinely 3D-embedded (curved) triangle
    // is a Line3D2
    const GeometryPointerType p_parent_geom = this->GetInputGeometry();
    const array_1d<double, 3> aux_local_coords = ZeroVector(3);
    array_1d<double, 3> parent_unit_normal = p_parent_geom->Normal(aux_local_coords);
    const double parent_normal_norm = norm_2(parent_unit_normal);
    KRATOS_ERROR_IF(parent_normal_norm < std::numeric_limits<double>::epsilon())
        << "Degenerate (zero-area) parent triangle: cannot compute a cut-segment normal." << std::endl;
    parent_unit_normal /= parent_normal_norm;

    const unsigned int n_interfaces = rInterfacesVector.size();
    const unsigned int n_int_pts_per_interface =
        n_interfaces > 0 ? (*rInterfacesVector[0]).IntegrationPointsNumber(IntegrationMethod) : 0;
    const unsigned int n_total_int_pts = n_interfaces * n_int_pts_per_interface;

    rInterfaceAreaNormalValues.clear();
    rInterfaceAreaNormalValues.reserve(n_total_int_pts);

    for (unsigned int i_interface = 0; i_interface < n_interfaces; ++i_interface) {
        const IndexedPointGeometryType& r_interface_geom = *rInterfacesVector[i_interface];
        const unsigned int n_int_pts = r_interface_geom.IntegrationPointsNumber(IntegrationMethod);

        const array_1d<double, 3> tangent = r_interface_geom.GetPoint(1) - r_interface_geom.GetPoint(0);
        array_1d<double, 3> area_normal;
        MathUtils<double>::CrossProduct(area_normal, tangent, parent_unit_normal);

        for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss) {
            rInterfaceAreaNormalValues.push_back(area_normal);
        }
    }
}

}; //namespace Kratos
