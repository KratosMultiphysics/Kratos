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
//

// System includes

// External includes

// Project includes
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"

namespace Kratos
{

/// Triangle2D3AusasModifiedShapeFunctions implementation
/// Default constructor
Triangle2D3AusasModifiedShapeFunctions::Triangle2D3AusasModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    AusasModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTriangleSplitter(Kratos::make_shared<DivideTriangle2D3<Node>>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTriangleSplitter->GenerateDivision();
    mpTriangleSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Triangle2D3AusasModifiedShapeFunctions::~Triangle2D3AusasModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Triangle2D3AusasModifiedShapeFunctions::Info() const {
    return "Triangle2D3N Ausas modified shape functions computation class.";
};

/// Print information about this object.
void Triangle2D3AusasModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Triangle2D3N Ausas modified shape functions computation class.";
};

/// Print object's data.
void Triangle2D3AusasModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Triangle2D3N Ausas modified shape functions computation class:\n";
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
const DivideGeometry<Node>::Pointer Triangle2D3AusasModifiedShapeFunctions::pGetSplittingUtil() const {
    return mpTriangleSplitter;
};

void Triangle2D3AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(
        rPosSideCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

void Triangle2D3AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(
        rNegSideCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

}; //namespace Kratos
