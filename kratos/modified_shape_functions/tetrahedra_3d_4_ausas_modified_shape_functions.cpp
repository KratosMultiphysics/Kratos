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
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

namespace Kratos
{

/// Tetrahedra3D4AusasModifiedShapeFunctions implementation
/// Default constructor
Tetrahedra3D4AusasModifiedShapeFunctions::Tetrahedra3D4AusasModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    AusasModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTetrahedraSplitter(Kratos::make_shared<DivideTetrahedra3D4>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTetrahedraSplitter->GenerateDivision();
    mpTetrahedraSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Tetrahedra3D4AusasModifiedShapeFunctions::~Tetrahedra3D4AusasModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Tetrahedra3D4AusasModifiedShapeFunctions::Info() const {
    return "Tetrahedra3D4N Ausas modified shape functions computation class.";
};

/// Print information about this object.
void Tetrahedra3D4AusasModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Tetrahedra3D4N Ausas modified shape functions computation class.";
};

/// Print object's data.
void Tetrahedra3D4AusasModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Tetrahedra3D4N Ausas modified shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    std::stringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};


// Returns a pointer to the splitting utility
const DivideGeometry::Pointer Tetrahedra3D4AusasModifiedShapeFunctions::pGetSplittingUtil() const {
    return mpTetrahedraSplitter;
};

void Tetrahedra3D4AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(
        rPosSideCondMatrix,
        mpTetrahedraSplitter->mEdgeNodeI,
        mpTetrahedraSplitter->mEdgeNodeJ,
        mpTetrahedraSplitter->mSplitEdges);
}

void Tetrahedra3D4AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(
        rNegSideCondMatrix,
        mpTetrahedraSplitter->mEdgeNodeI,
        mpTetrahedraSplitter->mEdgeNodeJ,
        mpTetrahedraSplitter->mSplitEdges);
}

}; //namespace Kratos
