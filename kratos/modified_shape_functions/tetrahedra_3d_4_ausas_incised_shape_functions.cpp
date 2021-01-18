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
#include "modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions.h"

namespace Kratos
{

/// Tetrahedra3D4AusasIncisedShapeFunctions implementation
/// Default constructor
Tetrahedra3D4AusasIncisedShapeFunctions::Tetrahedra3D4AusasIncisedShapeFunctions(const GeometryPointerType pInputGeometry,
        const Vector& rExtrapolatedNodalDistances, const Vector& rExtrapolatedEdgeRatios) :
    Tetrahedra3D4AusasModifiedShapeFunctions(pInputGeometry, rExtrapolatedNodalDistances), mExtraEdgeRatios(rExtrapolatedEdgeRatios) {};

/// Destructor
Tetrahedra3D4AusasIncisedShapeFunctions::~Tetrahedra3D4AusasIncisedShapeFunctions() {};

/// Turn back information as a string.
std::string Tetrahedra3D4AusasIncisedShapeFunctions::Info() const {
    return "Tetrahedra3D4N Ausas incised shape functions computation class.";
};

/// Print information about this object.
void Tetrahedra3D4AusasIncisedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Tetrahedra3D4N Ausas incised shape functions computation class.";
};

/// Print object's data.
void Tetrahedra3D4AusasIncisedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Tetrahedra3D4N Ausas incised shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    std::stringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tExtrapolated distance values: " << distances_buffer.str();
};

// Returns the nodal distances vector.
const Vector& Tetrahedra3D4AusasIncisedShapeFunctions::GetExtrapolatedEdgeRatios() const {
    return mExtraEdgeRatios;
};

}; //namespace Kratos
